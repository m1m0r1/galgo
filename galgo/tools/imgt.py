from __future__ import absolute_import
from __future__ import print_function
from argtools import command, argument
from builtins import filter, zip, range
from six import string_types
from itertools import islice
from collections import namedtuple, defaultdict
import logging
import re
import os
import time
from functools import wraps
from .. import sh
from ..utils import format as fmt
from ..utils import cached_property, collect_while, skip_until, isblank, make_not, blank_split, fill_text, Counter
from ..bioseq import dna_revcomp, dna_translate, Fasta, mask_unknown_seq
from ..bioseq import STOP_AA, iter_dna_translate
from ..blast import BlastTabFile
from ..imgt import IMGTNucFile, IMGTGenFile, IMGTAlnData, merge_gen_nuc_alns
import pandas as pd
import numpy as np
from joblib import Parallel, delayed
from tqdm import tqdm


class IMGTNuc:
    def __init__(self, name, aln):
        """
        name: A*01:01:01:02N
        aln: ATGGCCGTCATGGCGCCCCGAACCCTCCTC...
        """
        raw_sep_seq = aln.replace('.', '')  # remove deletion char
        raw_exon_seqs = raw_sep_seq.split('|')  # exons with '*'
        exon_frames = get_exon_frames(raw_exon_seqs)

        self.seq = raw_sep_seq.replace('|', '')
        self.exons = raw_exon_seqs
        self.exon_frames = exon_frames

@command.add_sub
@argument('nuc')
def imgt_nuc2exon(args):
    f = IMGTNucFile(args.nuc)
    print ('subtype', 'number', 'is_typed', 'lphase', 'rphase', 'length', 'seq', sep='\t')
    for subtype in f.subtypes:
        nuc = IMGTNuc(subtype, '|'.join(self.data[subtype]))
        for number, (seq, (lphase, rphase)) in enumerate(zip(nuc.exons, nuc.exon_frames), 1):
            is_typed = int('*' not in seq)
            print (subtype, number, is_typed, lphase, rphase, len(seq), seq, sep='\t')


@command.add_sub
@argument('nuc')
@argument('--bar', default='|')
@argument('-f', '--format', choices=['tab'], default='tab')
def imgt_nuc_view(args):
    f = IMGTNucFile(args.nuc)
    if args.format == 'tab':
        f.data.save()


@command.add_sub
@argument('gen')
@argument('--bar', default='|')
@argument('-f', '--format', choices=['tab'], default='tab')
def imgt_gen_view(args):
    f = IMGTGenFile(args.gen)
    if args.format == 'tab':
        f.data.save()


@command.add_sub
@argument('gen')
@argument('nuc')
@argument('-n', '--gene-name', help='gene name restriction')
def imgt_merge_gen_nuc(args):
    gen = IMGTGenFile(args.gen)
    nuc = IMGTNucFile(args.nuc)
    merged = merge_gen_nuc_alns(gen_data=gen.data, nuc_data=nuc.data, gene_name=args.gene_name)
    merged.save()


def msa_aln_impute(aln_data):
    """

    """
    #errors = {}  # {(X1, X2): error}
    #matches = {} # {(X1, X2): match}
    rates = {}  # {(X1, X2): match_rates}
    arr = aln_data.as_array()
    knowns = arr != '*'  # O(NM)
    known_indexes = [tuple(np.where(knowns[:, j])[0]) for j in range(knowns.shape[1])]  # O(NM)

    logging.info('Calculating matches and errors for all pairs')  # O(N^2 M)
    for i in tqdm(range(len(arr))):
        for j in range(i+1, len(arr)):
            ij_known = knowns[i] & knowns[j]
            i_arr = arr[i, ij_known]
            j_arr = arr[j, ij_known]
            match = (i_arr == j_arr).sum()
            #error = i_arr.shape[0] - match
            #matches[i, j] = matches[j, i] = match
            #errors[i, j] = errors[j, i] = error
            rates[i, j] = rates[j ,i] = 1. * match / i_arr.shape[0]

    logging.info('List patterns for closest subtypes')  # O(NM) + O(NP) ~ O(N^2 P)
    closests = {}  # {(k, ...): {u: k1}}  # known set => unknown index => closets known index
    for ks in set(known_indexes):
        kset = set(ks)
        closests[ks] = {u: max(ks, key=lambda k: 1. * rates[u, k]) for u in range(arr.shape[0]) if u not in kset}

    logging.info('Imputing bases from close subtypes')  # O(NM)
    for i, j in tqdm(zip(*np.where(~knowns))):
        ks = known_indexes[j]
        k = closests[ks][i]
        arr[i, j] = arr[k, j]

    assert (arr == '*').sum() == 0   # all the bases are imputed
    imputed = IMGTAlnData.from_array(aln_data.subtypes, aln_data.part_lens, arr)
    return imputed


@command.add_sub
@argument('msa_aln')
@argument('-f', '--format', choices=['msa_aln', 'fasta'])
def imgt_aln_impute(args):
    if args.format == 'fasta':
        data = IMGTAlnData.load_fasta(args.msa_aln)
    else:
        data = IMGTAlnData.load(args.msa_aln)
    data = msa_aln_impute(data)
    if args.format == 'fasta':
        data.save_as_fasta()
    else:
        data.save()


@command.add_sub
@argument('msa_aln')
def imgt_msa2aln(args):
    """
    This cannot reestimate exon-intron boundary '|'
    """
    data = IMGTAlnData.load(args.msa_aln)
    data = msa_aln_impute(data)
    data.save()

@command.add_sub
@argument('msa_aln')
@argument('-r', '--remove-padding', action='store_true')
def imgt_aln2fa(args):
    """
    padding: '.' -> '-'
    missing: '*' -> 'N'
    """
    data = IMGTAlnData.load(args.msa_aln)
    data.save_as_fasta(remove_padding=args.remove_padding)

#@command.add_sub
@argument('prot')
def imgt_prot(args):
    pass


def mask_left(seq, length, mask_char='N', alphabets='ACGT'):
    """
    >>> mask_left('NNAA-NATA', 3)
    'NNNN-NNTA'
    >>> mask_left('NNAA-NATA', 4)
    'NNNN-NNNA'
    >>> mask_left('NNAA-NATA', 1, alphabets='T')
    'NNAA-NANA'
    """
    def emitter():
        count = length
        for b in seq:
            if count and b in alphabets:
                yield mask_char
                count -= 1
            else:
                yield b

    it = emitter()
    return ''.join(it)


def mask_right(seq, length, mask_char='N', alphabets='ACGT'):
    """
    >>> mask_right('NNAA-NATA', 3)
    'NNAA-NNNN'
    >>> mask_right('NNAA-NATA', 4)
    'NNAN-NNNN'
    >>> mask_right('NNAA-NATA', 1, alphabets='T')
    'NNAA-NANA'
    """
    mseq = mask_left(reversed(seq), length, mask_char=mask_char, alphabets=alphabets)
    return ''.join(reversed(mseq))


@command.add_sub
@argument('gen_msa')
@argument('-l', '--exon-layout', required=True)
@argument('-c', '--cds-ranges', help='required in rule `kir_cds`')
@argument('-r', '--rule', required=True, choices=['hla_cds', 'kir_cds', 'kir_exon'])  # TODO HLA 6G (it requires correspondence table)
def ipd_gen_msa_mask(args):
    must_be_unique = True
    will_mask_utr = False
    if args.rule == 'hla_cds':
        def get_name(name):
            n1, n2 = name.split('*')
            n2e = ':'.join(n2.split(':')[:3])
            return n1 + '*' + n2e
    elif args.rule == 'kir_cds':
        will_mask_utr = True
        _cds_range_tab = pd.read_table(args.cds_ranges)
        # Required field for _cds_range_tab:
        # - al_name
        # - cds_start
        # - cds_end
        # - nuc_len
        _prot_name_clips = {row.prot_al_name: (row.cds_start, row.nuc_len - row.cds_end) for row in _cds_range_tab.itertuples()}
        def get_name(name):
            n1, n2 = name.split('*')
            n2e = n2[:3 + 2]
            return n1 + '*' + n2e
        def find_match_name(name):
            if name in _prot_name_clips:   # exact match
                return name
            for name1 in _prot_name_clips:
                if name.startswith(name1):  # prefix matched
                    return name1
        def mask_utr_seq(name, mseq):
            name = find_match_name(name)
            if name is None:  # find match here
                return
            lclip, rclip = _prot_name_clips[name]
            mseq = mask_left(mseq, lclip, mask_char=mask_char, alphabets='ACGT')
            mseq = mask_right(mseq, rclip, mask_char=mask_char, alphabets='ACGT')
            return mseq
    elif args.rule == 'kir_exon':
        must_be_unique = False
        def get_name(name):
            return name

    mask_char = 'N'
    tab = pd.read_table(args.exon_layout)
    regions = list(zip(tab['start'], tab['end']))
    def mask_seq(name, seq):
        parts = []
        s0 = 0
        for s, e in regions:
            parts.append(mask_char * (s - s0))
            parts.append(seq[s:e])
            s0 = e
        parts.append(mask_char * (len(seq) - s0))
        #parts.append(seq[s0:])  # falsy
        mseq = ''.join(parts)
        mseq = mask_unknown_seq(mseq, mask_char=mask_char)  # check this
        # TODO Here, further masking from 5' and 3' sides for cds masking
        if will_mask_utr:
            mseq1 = mask_utr_seq(name, mseq)
            if mseq1 is None:
                logging.warning('Cannot find valid CDS coordinate for %s', name)
            else:
                mseq = mseq1
        return mseq

    def seqs_are_matched(seq1, seq2, mask_char='N'):
        if len(seq1) != len(seq2):
            return False
        return all(b1 == mask_char or b2 == mask_char or b1 == b2 for b1, b2 in zip(seq1, seq2))

    def get_better_seq(seq1, seq2, mask_char='N'):  # first is primary
        n1 = seq1.count(mask_char)
        n2 = seq2.count(mask_char)
        return seq1 if n1 <= n2 else seq2

    def seq_is_superior(seq1, seq2, mask_char='N'):
        if len(seq1) < len(seq2):
            return False
        return all(b2 == mask_char or b1 == b2 for b1, b2 in zip(seq1, seq2))

    with open(args.gen_msa) as fp:
        fasta = Fasta(fp)
        masked_seqs = [(get_name(contig.name), mask_seq(contig.name, contig.seq)) for contig in fasta]
        seq_name_sets = defaultdict(set)
        name_seqs = {}
        emit_names = []
        name_seqs = {}
        for name1, seq in masked_seqs:
            if seq is None:
                logging.warning('Cannot find valid sequences for %s', name1)
                continue
            if seq not in seq_name_sets:
                emit_names.append(name1)  # use OrderedDict ?
                #name_seqs[name1] = seq
            seq_name_sets[seq].add(name1)
            if name1 not in name_seqs:  # seqs for the same name
                name_seqs[name1] = seq
            else:
                assert seqs_are_matched(name_seqs[name1], seq), 'The seqs for name {0} are incompatible (This is always required)'.format(name1)
                name_seqs[name1] = get_better_seq(name_seqs[name1], seq)
            if must_be_unique:
                assert len(seq_name_sets[seq]) == 1, 'seq name should be unique: {0}'.format(seq_name_sets[seq])
            if name1 in name_seqs:
                matched = seqs_are_matched(name_seqs[name1], seq)
                if not matched:
                    logging.warning((name1, [(i, b1, b2) for i, (b1, b2) in enumerate(zip(name_seqs[name1], seq)) if b1 != b2]))
                assert matched
                org_seq = name_seqs[name1]
                if seq_is_superior(org_seq, seq):
                    name_seqs[name1] = org_seq
                elif seq_is_superior(seq, org_seq):
                    name_seqs[name1] = seq
                else:
                    logging.warning('Incompatible seqs were found for %s', name1)
                    raise Exception('unknown case')
            #if len(name_seq_sets[name1]) > 1:
            #    seqs = list(name_seq_sets[name1])
            #    seq1 = seqs[0]
            #    seq2 = seqs[1]
            #    logging.warning([(i, b1, b2) for i, (b1, b2) in enumerate(zip(seq1, seq2)) if b1 != b2])
            #assert len(name_seq_sets[name1]) == 1, 'seq should be unique: {0}'.format(name1)
        #seq_name_count = Counter(n for nset in seq_name_sets.values() for n in nset)
        #assert all(c == 1 for c in seq_name_count.values()), 'All seqs and their names have to be one-to-one: {0}'.format(seq_name_count)

        for name in emit_names:
            seq = name_seqs[name]
            concat_name = name = ','.join(sorted(seq_name_sets[seq]))
            print ('>' + concat_name)
            print (seq)


def get_cds_range(nuc_seq, prot_seq):
    """
    >>> nuc_seq = 'ACTGGTTGCTAT'
    >>> prot_seq = 'TGCY'
    >>> l = get_cds_range(nuc_seq, prot_seq)
    >>> (l['cds_start'], l['cds_end'], l['cds_len'], l['nuc_len'])
    (0, 12, 12, 12)

    >>> nuc_seq = 'CTACTGGTTGCTATACT'
    >>> prot_seq = 'TGCY'
    >>> l = get_cds_range(nuc_seq, prot_seq)
    >>> (l['cds_start'], l['cds_end'], l['cds_len'], l['nuc_len'])
    (2, 14, 12, 17)

    >>> nuc_seq = 'ACTACTGGTTGCTATACT'
    >>> prot_seq = 'TGCY'
    >>> l = get_cds_range(nuc_seq, prot_seq)
    >>> (l['cds_start'], l['cds_end'], l['cds_len'], l['nuc_len'])
    (3, 15, 12, 18)

    >>> nuc_seq = 'ACTGGGACTGGTTGCTATACT'
    >>> prot_seq = 'TGCY'
    >>> l = get_cds_range(nuc_seq, prot_seq)
    >>> (l['cds_start'], l['cds_end'], l['cds_len'], l['nuc_len'])
    (6, 18, 12, 21)
    """
    prot_seq = prot_seq.replace('X', STOP_AA)
    assert STOP_AA not in prot_seq[:-1], 'Stop codon is only supposted to be at the end of prot_seq'
    cds_len = len(prot_seq.replace('$', '')) * 3
    min_len = len(prot_seq) * 3
    nuc_len = len(nuc_seq)
    assert nuc_len >= min_len, 'length of nuc_seq must be over that of amino acid sequnece'

    def matched(nuc_offset):
        iseq = islice(nuc_seq, nuc_offset, None)
        iaas = iter_dna_translate(iseq)
        for a1, a2 in zip(iaas, prot_seq):
            #logging.info((a1, a2))
            if a1 != a2:
                return False
        return True

    def find_offset():
        for s in range(nuc_len - min_len + 1):
        #for s in range(nuc_len):
            if matched(s):
                return s
        f1 = dna_translate(nuc_seq)
        f2 = dna_translate(nuc_seq[1:])
        f3 = dna_translate(nuc_seq[2:])
        logging.warning('aa: %s...', prot_seq)
        logging.warning('F1: %s...', f1)
        logging.warning('F2: %s...', f2)
        logging.warning('F3: %s...', f3)
        raise Exception('Cannot find matches')

    offset = find_offset()

    return {
            'cds_start': offset,
            'cds_end': offset + cds_len,
            'cds_len': cds_len,
            'nuc_len': nuc_len,
    }


@command.add_sub
@argument('nuc_fa')
@argument('prot_fa')
def ipd_kir_cds_range(args):
    """
    Note that X is stop codon in the prot_fa format

    Emits:
        - name
        - al_name
        - cds_start
        - cds_end
        - nuc_len
        - cds_len
    """
    #al_names = {}
    #if args.allele_list:   # alist tend to be incomplete
    #    al_names = {'IPD:' + name: 'KIR' + short_name for name, short_name in
    #            (line.rstrip().split(' ') for line in open(args.allele_list) if line.rstrip())}

    with open(args.nuc_fa) as fp1, \
         open(args.prot_fa) as fp2:
        nuc_fa = Fasta(fp1)
        prot_fa = Fasta(fp2)
        #assert set(nuc_fa.names) == set(prot_fa.names)
        nuc_names = nuc_fa.names
        prot_names = prot_fa.names
        nuc_al_names = {name: nuc_fa.get(name).name_line.split(' ')[1] for name in nuc_names}  # adhoc
        prot_al_names = {name: prot_fa.get(name).name_line.split(' ')[1] for name in prot_names}  # adhoc
        nuc_seqs = {name: nuc_fa.get(name).seq for name in nuc_names}
        prot_seqs = {name: prot_fa.get(name).seq for name in prot_names}
        print ('name', 'nuc_al_name', 'prot_al_name', 'cds_start', 'cds_end', 'nuc_len', 'cds_len', sep='\t')
        for name in nuc_names:
            nuc_al_name = nuc_al_names[name]
            if name not in prot_seqs:
                logging.warning('CDS range is not found for %s (%s)', name, nuc_al_name)
                continue
            prot_al_name = prot_al_names[name]
            #al_name = al_names.get(name, '.')
            nuc_seq = nuc_seqs[name]
            prot_seq = prot_seqs[name]
            try:
                cdsr = get_cds_range(nuc_seq, prot_seq)
                print (
                        name,
                        nuc_al_name,
                        prot_al_name,
                        cdsr['cds_start'],
                        cdsr['cds_end'],
                        cdsr['nuc_len'],
                        cdsr['cds_len'],
                        sep='\t')
            except Exception as e:
                logging.error('Cannot align %s (%s, %s). Length: nuc_len: %s, prot_len: %s', name, nuc_al_name, prot_al_name, len(nuc_seq), len(prot_seq))
                logging.error(e)


@command.add_sub
@argument('gen_msa')
@argument('nuc_msa')
@argument('--show-exon-layout', action='store_true')  # TODO rename to nuc-layout
def ipd_kir_align_msas(args):
    """
    """
    with open(args.gen_msa) as fp1, \
         open(args.nuc_msa) as fp2:
        gen_msa = Fasta(fp1)
        nuc_msa = Fasta(fp2)
        assert set(gen_msa.names) == set(nuc_msa.names)
        names = gen_msa.names
        gen_seqs = [gen_msa.get(name).seq for name in names]
        nuc_seqs = [nuc_msa.get(name).seq for name in names]
        msa = KIRAlignGenNuc(gen_seqs, nuc_seqs)
        positions = list(msa.iter_positions())
        for rec in positions:
            logging.info('gen: %s, nuc: %s, aln length: %s', rec['gen_start'], rec['nuc_start'], rec['length'])
        if args.show_exon_layout:
            print ('exon', 'start', 'end', 'length', sep='\t')
            for n, rec in enumerate(positions, 1):
                print (n, rec['gen_start'], rec['gen_start'] + rec['length'], rec['length'], sep='\t')
            return

        for gen_contig in gen_msa.contigs:
            name = gen_contig.name
            nuc_contig = nuc_msa.get(name)
            gen_seq = gen_contig.seq #get_clip_mask_seq(mask_char='N')
            nuc_seq = nuc_contig.seq #get_clip_mask_seq(mask_char='N')
            if gen_seq.replace('-', '') == nuc_seq.replace('-', ''):
                gen_seq = gen_contig.get_clip_mask_seq(mask_char='N')  # fill totally missing part
                logging.info('gen_seq of %s seems to be nuc_seq', name)
                parts = []
                prev_gen_start = 0
                for rec in positions:
                    gen_start = rec['gen_start']
                    length = rec['length']
                    parts.append('N' * (gen_start - prev_gen_start))
                    part_seq = gen_seq[gen_start:gen_start + length]
                    prev_gen_start = gen_start + length
                    parts.append(part_seq)
                parts.append('N' * (len(gen_seq) - prev_gen_start))
                seq = ''.join(parts)
            else:
                seq = gen_contig.get_clip_mask_seq(mask_char='N')
            print ('>', name, sep='')
            print (seq)


@command.add_sub
@argument('gen_msa')
@argument('nuc_msa')
def ipd_kir_profile_msas(args):
    """
    """
    with open(args.gen_msa) as fp1, \
         open(args.nuc_msa) as fp2:
        gen_msa = Fasta(fp1)
        nuc_msa = Fasta(fp2)

        common_set = set(gen_msa.names) & set(nuc_msa.names)
        if set(gen_msa.names) != set(nuc_msa.names):
            only_gens = set(gen_msa.names) - set(nuc_msa.names)
            only_nucs = set(nuc_msa.names) - set(gen_msa.names)
            logging.warning('%s names were only found in gen_msa: (%s)', len(only_gens), list(sorted(only_gens)))
            logging.warning('%s names were only found in nuc_msa: (%s)', len(only_nucs), list(sorted(only_nucs)))

        names = gen_msa.names
        gen_seqs = [gen_msa.get(name).seq for name in names if name in common_set]
        nuc_seqs = [nuc_msa.get(name).seq for name in names if name in common_set]
        msa = KIRAlignGenNuc(gen_seqs, nuc_seqs)
        positions = list(msa.iter_positions())
        for rec in positions:
            logging.info('gen: %s, nuc: %s, aln length: %s', rec['gen_start'], rec['nuc_start'], rec['length'])

        gen_exon_num = len(positions)

        print ('name', 'gene', 'al_name', 'prot_al_name', 'syn_al_name', 'has_genome', 'gen_len', 'nuc_len',
                'gen_exon_num',
                'exon_num',
                'min_exon',
                'max_exon',
                'exon_lens',
                sep='\t')
        for gen_contig in gen_msa.contigs:
            name = gen_contig.name
            nuc_contig = nuc_msa.get(name)
            gen_seq = gen_contig.seq #get_clip_mask_seq(mask_char='N')
            nuc_seq = nuc_contig.seq #get_clip_mask_seq(mask_char='N')

            gene = name.split('*')[0]
            al_name = name
            prot_name = gene + '*' + name.split('*')[1][:3]
            syn_name = gene + '*' + name.split('*')[1][:5]  # synonymous name
            gen_len = 0
            nuc_len = len(nuc_seq.replace('-', ''))

            exon_seqs = []
            prev_gen_start = 0
            for rec in positions:
                nuc_start = rec['nuc_start']
                length = rec['length']
                exon_seq = nuc_seq[nuc_start:nuc_start + length].replace('-', '')
                exon_seqs.append(exon_seq)

            exon_lens = [len(s) for s in exon_seqs]
            exon_exists = [i for i, s in enumerate(exon_seqs, 1) if s != '']
            min_exon = min(exon_exists)
            max_exon = max(exon_exists)
            exon_num = len(exon_exists)

            if gen_seq.replace('-', '') == nuc_seq.replace('-', ''):
                has_genome = False
                gen_seq = gen_contig.get_clip_mask_seq(mask_char='N')  # fill totally missing part
                logging.info('gen_seq of %s seems to be nuc_seq', name)
            else:
                has_genome = True
                gen_len = len(gen_seq.replace('N', '').replace('-', ''))
            print (name, gene, al_name, prot_name, syn_name, int(has_genome), gen_len, nuc_len,
                    gen_exon_num,
                    exon_num,
                    min_exon,
                    max_exon,
                    ','.join(str(l) for l in exon_lens),
                    sep='\t')


class KIRAlignGenNuc(object):
    """
    >>> res = list(KIRAlignGenNuc(['AAAA----TTTT'], ['A---T']).iter_positions())
    >>> res[0]['nuc_start'], res[0]['gen_start'], res[0]['length']
    (0, 3, 4)
    >>> res[1]['nuc_start'], res[1]['gen_start'], res[1]['length']
    (4, 8, 1)

    >>> res = list(KIRAlignGenNuc(['AAAA----TTTT'], ['A---TTT']).iter_positions())
    >>> res[0]['nuc_start'], res[0]['gen_start'], res[0]['length']
    (0, 3, 4)
    >>> res[1]['nuc_start'], res[1]['gen_start'], res[1]['length']
    (4, 8, 3)
    """

    def __init__(self, gen_seqs, nuc_seqs):
        assert len(gen_seqs) == len(nuc_seqs)
        self._gen_seqs = gen_seqs
        self._nuc_seqs = nuc_seqs

    def iter_positions(self):
        gen_start = 0
        nuc_start = 0
        length = 0
        gen_seqs = self._gen_seqs
        nuc_seqs = self._nuc_seqs

        # init
        gen_len = len(gen_seqs[0])
        nuc_len = len(nuc_seqs[0])
        seq_num = len(gen_seqs)
        L = []  # lower bounds of gen_seqs for nuc_seqs at index position
        U = []  # upper bounds of gen_seqs for nuc_seqs at index position
        def match_all(nuc_pos, gen_pos):
            return all(nuc_seqs[k][nuc_pos] == gen_seqs[k][gen_pos] for k in range(seq_num))

        s = 0
        for i in range(nuc_len):
            for j in range(s, gen_len):
                if match_all(i, j):
                    L.append(j)
                    s = j + 1
                    break
        logging.debug('Lower bounds: %s', L)
        logging.info('Lower bounds: %s, nuc_len: %s', len(L), nuc_len)
        if len(L) < nuc_len:
            raise Exception('nuc_seqs are not contained in gen_seqs')

        e = gen_len - 1  # TODO check again
        for i in range(nuc_len-1, -1, -1):
            for j in range(e, -1, -1):
                if match_all(i, j):
                    U.insert(0, j)
                    e = j - 1
                    break
        logging.debug('Upper bounds: %s', U)
        logging.info('Upper bounds: %s, nuc_len: %s', len(U), nuc_len)
        if len(U) < nuc_len:
            raise Exception('nuc_seqs are not contained in gen_seqs')

        # greedy search matched blocks
        nuc_start = 0
        gen_start = L[0]
        V_prev = set()
        for a in range(nuc_len):
            V = {b for b in range(max(gen_start, L[a]), U[a]+1) if (b-1 in V_prev if V_prev else True) and match_all(a, b)}
            if V:
                V_prev = V
            else:
                length = a - nuc_start
                yield {'nuc_start': nuc_start, 'gen_start': min(V_prev) + 1 - length, 'length': length}
                nuc_start = a
                gen_start = min(V_prev) + 1
                V_prev = {b for b in range(max(gen_start, L[a]), U[a]+1) if match_all(a, b)}

        length = nuc_len - nuc_start
        if V_prev:
            logging.info({'nuc_start': nuc_start, 'gen_start': min(V_prev) + 1 - length, 'length': length})
            yield {'nuc_start': nuc_start, 'gen_start': min(V_prev) + 1 - length, 'length': length}


@command.add_sub
@argument('gen_msa')
@argument('nuc_msa')
def ipd_kir_merge(args):
    """
    """
    with open(args.gen_msa) as fp1, \
         open(args.nuc_msa) as fp2:
        gen_msa = Fasta(fp1)
        nuc_msa = Fasta(fp2)
        assert set(gen_msa.names) == set(nuc_msa.names)

        for gen_contig in gen_msa.contigs:
            name = gen_contig.name
            nuc_contig = nuc_msa.get(name)
            print ('>', gen_contig.name, sep='')
            gen_seq = gen_contig.seq
            nuc_seq = nuc_contig.seq
            if gen_seq.replace('-', '') == nuc_seq.replace('-', ''):
                logging.info('gen_seq of %s seems to be nuc_seq', name)
                seq = gen_seq.replace('-', 'N')
            else:
                seq = gen_contig.get_clip_mask_seq(mask_char='N')
            print (seq)


@command.add_sub
@argument('contigs')
@argument('-d', '--imgt-dir', default=None)
@argument('-g', '--imgt-gen', default=None)
@argument('-l', '--allele-list', default=None)
@argument('-o', '--out-dir', default='imgt_classify_hla')
@argument('-t', '--target-genes', nargs='+', default=['A', 'B', 'C'])
@argument('-f', '--force-rerun', action='store_true')
@argument('--dry-run', action='store_true')
@argument('--EST2GENOME', default='est2genome')
@argument('-j', '--njobs', type=int, default=8)
def imgt_classify_hla(args):
    """
    Currently, classical imgt type only

    Inputs:
        {imgt_dir}/fasta/{gene}_nuc.fasta
        gen.fasta
        required files under {imgt_aln_dir}
        *_nuc.txt is required for identify the exon number contained in cDNA (*_nuc.fasta)
        *_nuc.txt

    Outputs:
        {outdir}/{gene}_nuc.exons.txt
        {outdir}/{gene}_nuc.fa
        {outdir}/{gene}_nuc.fa
        {outdir}/gen.blastdb
        {outdir}/contig.gen.blastn
        {outdir}/contig.gen.blastn.bests
        {outdir}/contigs/{contig_name}.fa
        {outdir}/contigs/{contig_name}.{gene}.est2genome
        {outdir}/contig_summary1.txt
        {outdir}/contig_summary2.txt

    # gene type
    # 4 digit (protein level)
    # 6 digit (exon level)     <= maybe detection of exon 8 will be failed
    # 6G digit (exon 2 and 3) wmda/hla_nom_G.txt
    # 8 digit (outside) if exists
    # variant on exon
    # variant on intron
    """
    contigs_fasta = args.contigs
    gen_fasta_org = args.imgt_gen
    imgt_dir = args.imgt_dir
    out_dir = args.out_dir
    sh.dry_run = args.dry_run

    sh.call(['mkdir', '-p', out_dir])
    gen_blastdb = '{dir}/gen.blastdb'.format(dir=out_dir)
    gen_blastdb_nhr = gen_blastdb + '.nhr'

    if args.allele_list:
        allele_list = args.allele_list
    else:
        allele_list = '{dir}/Allelelist.txt'.format(dir=args.imgt_dir)  # default value

    if args.imgt_gen:
        gen_fasta_org = args.imgt_gen
    else:
        gen_fasta_org = '{dir}/fasta/hla_gen.fasta'.format(dir=args.imgt_dir)  # default value

    gen_fasta = '{dir}/gen.fa'.format(dir=out_dir)
    contig_gen_blast = '{dir}/contigs.gen.blastn'.format(dir=out_dir)
    contig_gen_blast_bests = '{dir}/contigs.gen.blastn.bests'.format(dir=out_dir)
    contig_cdna_dir = '{dir}/contigs'.format(dir=out_dir)
    closest_exons = '{dir}/contig.closest_exons.txt'.format(dir=out_dir)
    closest_genomes = '{dir}/contig.closest_genomes.txt'.format(dir=out_dir)
    contig_summary = '{dir}/contig.summary.txt'.format(dir=out_dir)
    contig_nuc_fa_tmpl = '{dir}/contig.summary.{{gene}}_nuc.fa'.format(dir=out_dir)
    contig_gen_fa_tmpl = '{dir}/contig.summary.{{gene}}_gen.fa'.format(dir=out_dir)
    exon_counts = {'A': 8, 'B': 7, 'C': 8, 'H': 8}   # TODO more systematic way

    if imgt_dir is None:
        logging.error('imgt_dir is required for preprocss')
        raise Exception('No data directory')

    if args.force_rerun or not os.path.exists(gen_fasta) or os.path.getsize(gen_fasta) < 1:
        require(gen_fasta_org)
        galgo(r'''fa_select {fasta} -n $(cat {list} | awk '{{print "HLA:" $1}}' | tr '\n' ' ') > {out_fasta}'''
                .format(fasta=gen_fasta_org, list=allele_list, out_fasta=gen_fasta),
                stdout=2)

    if args.force_rerun or not os.path.exists(gen_blastdb_nhr):
        require(gen_fasta)
        sh.call('makeblastdb -in {fasta} -out {db} -dbtype nucl -parse_seqids'.format(fasta=gen_fasta, db=gen_blastdb))
        #sh.call_if('', 'wc -l {fasta}')

    if args.force_rerun or not os.path.exists(contig_gen_blast):
        require(gen_blastdb_nhr)
        sh.call("""blastn -db {db} -num_threads {njobs} -query {fasta} -evalue 0.1 -outfmt '7 std gaps qlen slen qseq sseq'"""
            .format(db=gen_blastdb, fasta=contigs_fasta, njobs=args.njobs),
            stdout=contig_gen_blast)

    #@after(indexing)
    if args.force_rerun or not os.path.exists(contig_gen_blast_bests):
        require(contig_gen_blast)
        require(allele_list)
        if not args.dry_run:
            _create_best_contigs(contig_gen_blast, allele_list, contig_gen_blast_bests)

    for gene in args.target_genes:
        nuc = '{dir}/{gene}_nuc.fa'.format(dir=out_dir, gene=gene)
        if args.force_rerun or not os.path.exists(nuc) or os.path.getsize(nuc) < 1:
            galgo(r'''fa_select {org_fasta} -n $(cat {list} | awk '{{print "HLA:" $1}}' | tr '\n' ' ') > {out_fasta}'''
                    .format(
                        org_fasta='{imgt_dir}/fasta/{gene}_nuc.fasta'.format(imgt_dir=imgt_dir, gene=gene),
                        list=allele_list,
                        out_fasta=nuc))

    #    {outdir}/contigs/{contig_name}.fa
    #    {outdir}/contigs/{contig_name}.{gene}.est2genome
    if not args.dry_run:
        require(contig_gen_blast_bests)
        sh.call('mkdir -p {dir}'.format(dir=contig_cdna_dir))
        env = {'out_dir': out_dir, 'args': args, 'contig_cdna_dir': contig_cdna_dir}   # workaround
        Parallel(n_jobs=args.njobs, backend='threading')(delayed(run_est2genome)(name, gene, env) for name, gene in extract_contig_genes(contig_gen_blast_bests))

    def run_classify_hla():
        require(contig_gen_blast_bests)
        require(allele_list)
        def find_closest_exons(name, gene):
            est2genome_parsed = '{dir}/{name}.{gene}_nuc.est2genome.txt'.format(dir=contig_cdna_dir, name=name, gene=gene)
            require(est2genome_parsed)
            tab = pd.read_table(est2genome_parsed)
            tab = tab[tab.category == 'span'].sort_values('edit_distance').sort_values('number', ascending=False).reset_index()
            #logging.info(tab)
            min_edit = tab['edit_distance'].loc[0]
            logging.info((name, gene, min_edit, tab.edit_distance.min()))
            return tab[tab.edit_distance == min_edit]

        def find_full_matched_exons(name, gene):  # TODO check if full length nuc was matched (maybe you have to read nuc lengths from fasta)
            est2genome_parsed = '{dir}/{name}.{gene}_nuc.est2genome.txt'.format(dir=contig_cdna_dir, name=name, gene=gene)
            require(est2genome_parsed)
            tab = pd.read_table(est2genome_parsed)
            tab = tab[(tab.category == 'span')
                    & (tab.number == exon_counts[gene])
                    & (tab.regular_intron == exon_counts[gene] - 1)].sort_values('edit_distance').reset_index()

            #logging.info(tab)
            min_edit = tab['edit_distance'].loc[0]
            logging.info((name, gene, min_edit, tab.edit_distance.min()))
            tab = tab[tab.edit_distance == min_edit]

            # TODO should check aa before edit distance?
            tab['qseq_aa'] = tab.qseq.map(lambda x: dna_translate(x.replace('|', '').replace('-', '').upper()))  # db seq
            tab['rseq_aa'] = tab.rseq.map(lambda x: dna_translate(x.replace('|', '').replace('-', '').upper()))  # contig seq
            tab['aa_match'] = tab['qseq_aa'] == tab['rseq_aa']
            return tab


        contig_genes = list(extract_contig_genes(contig_gen_blast_bests))

        def get_matched_exons():
            #tab = pd.concat(find_closest_exons(name, gene) for name, gene in contig_genes)
            tab = pd.concat(find_full_matched_exons(name, gene) for name, gene in contig_genes)
            return tab

        def get_top_genomes():
            tab1 = pd.read_table(contig_gen_blast_bests)
            #tab = pd.concat(find_closest_exons(name, gene) for name, gene in extract_contig_genes(contig_gen_blast_bests))
            def gen_tops():
                for sname, tab in tab1.groupby('qname'):
                    tab = tab.sort_values('edit_distance').reset_index()
                    #logging.info(tab)
                    min_edit = tab['edit_distance'].loc[0]
                    #logging.info(min_edit)
                    yield tab[tab.edit_distance == min_edit]

            tab = pd.concat(gen_tops())
            return tab

        def get_consistent_genomes(qname_subtype_6ds):
            tab1 = pd.read_table(contig_gen_blast_bests)
            #tab = pd.concat(find_closest_exons(name, gene) for name, gene in extract_contig_genes(contig_gen_blast_bests))
            tab1['subtype_6d'] = tab1['subtype'].apply(lambda x: ':'.join(x.split(':')[:3]))

            def gen_tops():
                for qname, tab in tab1.groupby('qname'):
                    subtype_6ds = qname_subtype_6ds.get(qname, [])
                    yield pd.concat(
                            tab[tab['subtype_6d'] == s].sort_values(['edit_distance', 'aln_len'], ascending=[True, False]).nsmallest(1, columns=['edit_distance'])
                            for s in subtype_6ds)

            tab = pd.concat(gen_tops())
            return tab

        tab = get_matched_exons()
        alist = _read_allele_list(allele_list)
        tab = tab.merge(alist, left_on='qname', right_on='id')#, rsuffix='.1')
        tab = tab[tab.columns.difference(['qseq', 'rseq'])] #+ tab[['qseq', 'rseq']]
        tab.sort_values('rname', inplace=True)
        tab.to_csv(closest_exons, sep='\t', index=False)
        logging.info('Create %s', closest_exons)

        # summary contig info
        tab1 = tab[['rname', 'qname', 'subtype', 'edit_distance', 'subtype_6d', 'subtype_4d', 'aa_match']]
        #tab1['subtype_6d'] = tab1['subtype'].apply(lambda x: ':'.join(x.split(':')[:3]))
        #tab1.loc[:, 'subtype_6d'] = tab1['subtype'].apply(lambda x: ':'.join(x.split(':')[:3]))
        #tab1.loc[:, 'subtype_4d'] = tab1['subtype'].apply(lambda x: ':'.join(x.split(':')[:2]))

        contig_info = {}
        # gene
        for name, gene in contig_genes:
            contig_info[name] = {}
            contig_info[name]['gene'] = gene

        # exon
        for name, t in tab1.groupby('rname'):
            contig_info[name]['edit_6d'] = ed = list(t.edit_distance)[0]  # only one for edit 6d
            contig_info[name]['closest_6d'] = closest_6d = list(sorted(t.subtype_6d.unique()))
            contig_info[name]['closest_4d'] = closest_4d = list(sorted(t.subtype_4d.unique()))
            contig_info[name]['aa_match'] = aa_match = [list(t[t.subtype_4d == sub].aa_match)[0] for sub in closest_4d]
            if ed == 0: # exact
                contig_info[name]['exact_6d'] = contig_info[name]['closest_6d'][0]
            if len(aa_match) > 0 and np.any(aa_match): # exact
                idx_match = np.argmax(aa_match)
                contig_info[name]['exact_4d'] = closest_4d[idx_match]

        #tab = get_top_genomes()
        #tab.sort_values('qname', inplace=True)
        #tab.to_csv(closest_genomes, sep='\t', index=False)
        #logging.info('Create %s', closest_genomes)
        contig_closests = dict((name, d['closest_6d']) for name, d in contig_info.items())

        tab = get_consistent_genomes(contig_closests)
        tab.sort_values('qname', inplace=True)
        tab.to_csv(closest_genomes, sep='\t', index=False)
        logging.info('Create %s', closest_genomes)

        #tab1 = tab[['qname', 'subtype', 'edit_distance']]
        tab['exact_6d'] = tab1['subtype'].apply(lambda x: ':'.join(x.split(':')[:3]))

        # genomes
        for name, t in tab.groupby('qname'):
            t = t.sort_values(['edit_distance', 'aln_len'], ascending=[True, False])
            contig_info[name]['edit_8d'] = eds = list(t.edit_distance)
            contig_info[name]['closest_8d'] = list(t.subtype)
            contig_info[name]['ref_left_8d'] = list(t.qlclip)  # reference is query here
            contig_info[name]['ref_right_8d'] = list(t.qrclip)
            contig_info[name]['q2s_edit'] = list(t.q2s_edit)  # reference is query here
            contig_info[name]['s2q_edit'] = list(t.s2q_edit)
            if eds and eds[0] == 0: # exact (not only one, but select one with max aln_len)
                contig_info[name]['exact_8d'] = contig_info[name]['closest_8d'][0]

        # list matched imgt types
        for name, info in contig_info.items():
            if len(info.get('closest_8d', [])):  # if closest subtypes exist, restrict to them
                match_digits = [sub.split(':') for sub in info['closest_8d']]
            else:
                match_digits = [sub.split(':') for sub in info['closest_6d']]
            logging.info('match_digits for %s: %s', name, match_digits)
            is_sub = alist.subtype.map(lambda x: any(x.split(':')[:len(d)] == d for d in match_digits))
            contig_info[name]['subs'] = list(alist[is_sub]['subtype'])
            contig_info[name]['subs_id'] = list(alist[is_sub]['id'])


        # finalize
        def iter_contig_info():   # for pandas
            for name, gene in contig_genes:
                d = contig_info[name]
                yield {
                        'refname': name,
                        'gene': gene,
                        'exact_4d': d.get('exact_4d', ''),
                        'exact_6d': d.get('exact_6d', ''),
                        'exact_8d': d.get('exact_8d', ''),
                        'n_closest_4d': len(d.get('closest_4d', [])),
                        'closest_4d': ','.join(d.get('closest_4d', [])),
                        'n_closest_6d': len(d.get('closest_6d', [])),
                        'closest_6d': ','.join(d.get('closest_6d', [])),
                        'edit_6d': d.get('edit_6d', -1),
                        'n_closest_8d': len(d.get('closest_8d', [])),
                        'closest_8d': ','.join(d.get('closest_8d', [])),
                        'ref_right_8d': ';'.join(map(str, d.get('ref_right_8d', []))),
                        'ref_left_8d': ';'.join(map(str, d.get('ref_left_8d', []))),
                        'edit_8d': ';'.join(map(str, d.get('edit_8d', []))),
                        'edit_8d_q2s': ';'.join(map(str, d.get('q2s_edit', []))),
                        'edit_8d_s2q': ';'.join(map(str, d.get('s2q_edit', []))),
                        'subs': ','.join(d.get('subs', [])),
                        'subs_id': ','.join(d.get('subs_id', [])),
                        }

        tab = pd.DataFrame(iter_contig_info(),
                columns='refname gene exact_4d exact_6d exact_8d n_closest_4d closest_4d n_closest_6d closest_6d edit_6d n_closest_8d closest_8d ref_left_8d ref_right_8d edit_8d edit_8d_q2s edit_8d_s2q subs subs_id'.split(' '))
        tab.to_csv(contig_summary, sep='\t', index=False)

    def run_summary_fasta():
        require(contig_summary)
        require(allele_list)
        tab = pd.read_table(contig_summary)
        imgt6d_refnames = {}  # {imgt_6d: first refnames}
        for row in tab[['refname', 'closest_6d', 'subs']].itertuples():
            if not row.closest_6d:
                continue
            names = row.closest_6d.split(',')
            for name in names:
                imgt6d_refnames.setdefault(name, row.refname)

        alist = _read_allele_list(allele_list)
        subtype_hla_ids = defaultdict(lambda : defaultdict(list))   # e.g. {'A': {'A*01:01:01': ['HLA00001', 'HLA02169', 'HLA14798']}}
        hla_id_set = set([sid for subs in tab.subs_id if isinstance(subs, string_types) for sid in subs.split(',')])  # only create fasta from contigs in subs_id
        for row in alist[['id', 'gene', 'subtype_6d']].itertuples():
            if row.id in hla_id_set: #row.subtype_6d in imgt6d_refnames and row.id in hla_id_set:
                subtype_hla_ids[row.gene][row.subtype_6d].append(row.id)

        # create gen file (contig name is like HLA:HLA00001)
        for gene in args.target_genes:
            fasta_ids = ['HLA:' + id for ids in subtype_hla_ids.get(gene, {}).values() for id in ids]
            out_fasta = contig_gen_fa_tmpl.format(gene=gene)
            if not fasta_ids:
                continue
            galgo('fa_select {fasta} -n {names} > {out_fasta}'.format(names=' '.join(fasta_ids), fasta=gen_fasta, out_fasta=out_fasta))

        # create nuc file (contig name is like HLA00001)
        for gene in args.target_genes:
            sh.call(': > {out_fasta}'.format(out_fasta=contig_nuc_fa_tmpl.format(gene=gene)))

        for subtype_6d, refname in imgt6d_refnames.items():
            gene = subtype_6d[0]  # A, B, C, H
            est_fasta = '{dir}/{name}.{gene}_nuc.est2genome.fa'.format(dir=contig_cdna_dir, name=refname, gene=gene)
            out_fasta = contig_nuc_fa_tmpl.format(gene=gene)
            hla_ids = subtype_hla_ids[gene][subtype_6d]
            galgo('fa_select {fasta} -n {hla_ids} >> {out_fasta}'.format(fasta=est_fasta, hla_ids=' '.join(hla_ids), out_fasta=out_fasta))

    if not args.dry_run:
        run_classify_hla()
        run_summary_fasta()


class IMGTHLANamer:
    """
    >>> namer = IMGTHLANamer()
    >>> namer.gen_name('A')
    'A*x1'
    >>> namer.gen_name('A')
    'A*x2'
    >>> namer.gen_name('A', 'A*01')
    'A*01:x1'
    >>> namer.gen_name('B')
    'B*x1'
    """

    def __init__(self):
        self._offsets = defaultdict(int)  # A*01:01 : 1
        self._prefix = 'x'

    def gen_name(self, gene, subtype=None):
        key = (gene, subtype)
        self._offsets[key] += 1
        last = self._prefix + str(self._offsets[key])
        if subtype is None:
            return gene + '*' + last
        else:
            return subtype + ':' + last


def isnan(obj):
    return isinstance(obj, float) and np.isnan(obj)

@command.add_sub
@argument('contig_summary')
def imgt_name(args):
    tab = pd.read_table(args.contig_summary)
    namer = IMGTHLANamer()
    print ('refname', 'imgt_name', sep='\t')
    for row in tab.itertuples():
        gene = row.gene
        subtype = None
        if not isnan(row.exact_4d):
            subtype = row.exact_4d
        if not isnan(row.exact_6d):
            subtype = row.exact_6d
        if not isnan(row.exact_8d):
            subtype = row.exact_8d
        name = namer.gen_name(gene, subtype)
        print (row.refname, name, sep='\t')


@command.add_sub
@argument('allele_list')
@argument('subtype', help='e.g. A*01:01:01')
def imgt_subtypes(args):
    require(args.allele_list)
    alist = _read_allele_list(args.allele_list)
    digits = args.subtype.split(':')
    for row in alist[['id', 'subtype']].itertuples():
        if row.subtype.split(':')[:len(digits)] == digits:  # this is 
            print (row.id, row.subtype, sep='\t')


def galgo(args, **kwds):
    import sys
    if isinstance(args, string_types):
        args = ' '.join((sys.executable, sys.argv[0], args))
    else:
        args = [sys.executable, sys.argv[0]] + list(args)
    return sh.call(args, **kwds)



def run_est2genome(name, gene, env):
    out_dir = env['out_dir']
    contig_cdna_dir = env['contig_cdna_dir']
    args = env['args']
    contigs_fasta = args.contigs
    nuc = '{dir}/{gene}_nuc.fa'.format(dir=out_dir, gene=gene)
    contig = '{dir}/{name}.fa'.format(dir=contig_cdna_dir, name=name)
    est2genome = '{dir}/{name}.{gene}_nuc.est2genome'.format(dir=contig_cdna_dir, name=name, gene=gene)
    est2genome_parsed = '{dir}/{name}.{gene}_nuc.est2genome.txt'.format(dir=contig_cdna_dir, name=name, gene=gene)
    est2genome_fasta = '{dir}/{name}.{gene}_nuc.est2genome.fa'.format(dir=contig_cdna_dir, name=name, gene=gene)
    logging.info('target: %s', est2genome)
    if args.force_rerun or not os.path.exists(est2genome) or os.path.getsize(est2genome) < 1:
        require(contigs_fasta)
        galgo('fa_select {fasta} -n {name} > {contig}'.format(fasta=contigs_fasta, name=name, contig=contig))
        require(contig)
        sh.call('''
        {EST2GENOME} -estsequence {nuc} -genomesequence {contig} -intronpenalty 40 -splicepenalty 1 -mismatch 1 -align 1 -usesplice -outfile {outfile}
        '''.format(EST2GENOME=args.EST2GENOME, nuc=nuc, contig=contig, outfile=est2genome))

    if args.force_rerun or not os.path.exists(est2genome_parsed) or os.path.getsize(est2genome_parsed) < 1:
        require(est2genome)
        galgo('''est2genome_parse {est2genome} > {out}'''.format(est2genome=est2genome, out=est2genome_parsed))

    if args.force_rerun or not os.path.exists(est2genome_fasta) or os.path.getsize(est2genome_fasta) < 1:
        require(est2genome)
        galgo('''est2genome_parse {est2genome} -f fasta > {out}'''.format(est2genome=est2genome, out=est2genome_fasta))

def extract_contig_genes(blast_bests):
    tab = pd.read_table(blast_bests)
    print (tab[tab['rank'] == 0])
    for row in tab[tab['rank'] == 0][['qname', 'gene']].itertuples():
        yield row.qname, row.gene


def _read_allele_list(allele_list, fasta_prefix='HLA:'):
    alist = pd.read_table(allele_list, sep=' ', names=['id', 'subtype'])
    alist['fasta_id'] = fasta_prefix + alist['id']
    alist['gene'] = alist.subtype.map(lambda x: x.split('*')[0])
    alist['subtype_8d'] = alist.subtype.map(lambda x: ':'.join(x.split(':')[:4]))
    alist['subtype_6d'] = alist.subtype.map(lambda x: ':'.join(x.split(':')[:3]))
    alist['subtype_4d'] = alist.subtype.map(lambda x: ':'.join(x.split(':')[:2]))
    alist['subtype_2d'] = alist.subtype.map(lambda x: ':'.join(x.split(':')[:1]))
    return alist

def _create_best_contigs(blastn, allele_list, bests):
    keys = 'qname rank sname gene subtype aln_len qlen slen identity edit_distance mismatches gaps slclip srclip qlclip qrclip left_overhang right_overhang'.split(' ')
    tab = _get_best_contigs(blastn)
    alist = _read_allele_list(allele_list)
    tab = tab.merge(alist, left_on='sname', right_on='fasta_id')#, rsuffix='.1')
    sseq = np.where(tab.is_reverse, tab.sseq.apply(dna_revcomp), tab.sseq)
    sstart = np.where(tab.is_reverse, tab.send, tab.sstart)
    def encode_edit(edit_info):
        return ','.join(('{0}:{1}>{2}'.format(*e) for e in edit_info))
    tab1 = tab[keys]
    tab1.loc[:, 's2q_edit'] = [encode_edit(calc_edit_info(q, s, start)) for q, s, start in zip(tab.qseq, sseq, sstart)]
    tab1.loc[:, 'q2s_edit'] = [encode_edit(calc_edit_info(s, q, start)) for s, q, start in zip(sseq, tab.qseq, tab.qstart)]

    #tab.to_csv('/dev/stdout', sep='\t', index=False)
    tab1.to_csv(bests, sep='\t', index=False)


def _get_best_contigs(blastn):
    min_sratio = .95
    max_overhang = 20
    #max_edit_distance = 20
    max_edit_distance = 50

    require(blastn)

    with open(blastn) as reader: #, open(summary, 'w+') as writer:
        blast = BlastTabFile(reader)
        def iter_tabs():
            for query, tab in blast.iter_query_tables():
                match_cond = ((tab.sratio >= min_sratio)
                    & (tab.left_overhang <= max_overhang)
                    & (tab.right_overhang <= max_overhang)
                    & (tab.edit_distance <= max_edit_distance))
                tab = tab[match_cond]
                tab = tab.sort_values(['edit_distance', 'aln_len'], ascending=[True, False])
                tab.reset_index(drop=True, inplace=True)
                tab.loc[:, 'rank'] = tab.index
                yield tab

        return pd.concat(iter_tabs()) #ignore_index=True)


def calc_edit_info(query_aln, ref_aln, ref_start=1): #TODO simplify
    """
    #        123 4567   890123456 789012 # offset
    >>> q = 'ACGTCGTTTGG--GACC-GTTTATTCA'
    >>> r = 'AAG-CGTT---ATGTCCAGT-TGTTAA'
    >>> l = list(calc_edit_info(q, r, 3001))
    >>> l[0]
    (3002, 'A', 'C')
    >>> l[1]
    (3004, '', 'T')
    >>> l[2]
    (3008, 'AT', 'TGG')
    >>> l[3]
    (3011, 'T', 'A')
    >>> l[4]
    (3014, 'A', '')
    >>> l[5]
    (3017, '', 'T')
    >>> l[6]
    (3018, 'G', 'A')
    >>> l[7]
    (3021, 'A', 'C')
    """
    assert len(query_aln) == len(ref_aln)
    edits = []
    offset = 0
    apos = ref_start
    for r, q in zip(ref_aln, query_aln):
        if r != q:
            edits.append((apos, offset, r.replace('-', ''), q.replace('-', '')))
        if r == '-':
            offset += 1
        apos += 1

    prev = -float('inf')
    outs = []
    for apos, offset, r, q in edits:
        if apos - prev == 1:
            outs[-1]['ref'] += r
            outs[-1]['query'] += q
        else:
            outs.append({'pos': apos - offset, 'ref': r, 'query': q})
        prev = apos

    return [(out['pos'], out['ref'], out['query']) for out in outs]


def require(file):
    if not os.path.exists(file):
        raise Exception('File does not exist: {0}'.format(file))


@command.add_sub
@argument('est2genome')
@argument('-f', '--format', choices=['tab', 'fasta'], default='tab', help='output format')
def est2genome_parse(args):
    """
    span
        alignment bases
        '|': intron
        # TODO check intron category (regular or irregular)
    """
    with open(args.est2genome) as fp:
        aln_keys = ['edit_distance', 'mismatches', 'gaps', 'qseq', 'rseq']
        if args.format == 'tab':
            print (*('qname qstart qend category number rname rstart rend score match_rate regular_intron'.split(' ') + aln_keys), sep='\t')
            for rec in Est2GenomeRecord.parse(fp):
                aln = [getattr(rec.span, k) for k in aln_keys]
                regular_intron = len([i for i in rec.introns if i.type in '+-'])
                print (rec.qname, rec.span.qstart, rec.span.qend, 'span', len(rec.exons), rec.span.rname, rec.span.rstart, rec.span.rend, rec.span.score, rec.span.match_rate, regular_intron, *aln, sep='\t')
                for num, rec in enumerate(rec.exons, 1):
                    aln = [getattr(rec, k) for k in aln_keys]
                    print (rec.qname, rec.qstart, rec.qend, 'exon', num, rec.rname, rec.rstart, rec.rend, rec.score, rec.match_rate, 0, *aln, sep='\t')
        elif args.format == 'fasta':
            for rec in Est2GenomeRecord.parse(fp):
                print ('>' + rec.qname)
                seqs = []
                last_rend = None
                for num, rec in enumerate(rec.exons, 1):
                    if last_rend is not None:
                        fill_len = rec.rstart - last_rend
                        seqs.append('N' * fill_len)
                    seqs.append(rec.qseq)
                    last_rend = rec.rend
                print (''.join(seqs))



class Est2GenomeRecord:
    #Span = namedtuple('Span', 'score match_rate rstart rend rname qstart qend qname rseq qseq')
    #Exon = namedtuple('Exon', 'score match_rate rstart rend rname qstart qend qname rseq qseq')
    class Feature(namedtuple('Feature', 'score match_rate rstart rend rname qstart qend qname rseq qseq')):
        #__slots__ = ('mismatches', 'edit_distance')
        def __init__(self, *args, **kwds):
            #super(Est2GenomeRecord.Feature, self).__init__(*args, **kwds)
            self.mismatches = self.gaps = self.edit_distance = None
            if self.rseq is not None:
                self.mismatches = self.gaps = 0
                for r, q in zip(self.rseq, self.qseq):
                    if r == '|':
                        continue
                    elif r == '-' or q == '-':
                        self.gaps += 1
                    elif r != q:
                        self.mismatches += 1
                self.edit_distance = self.mismatches + self.gaps

    Intron = namedtuple('Intron', 'type score match_rate rstart rend rname')

    """
    """

    _intron_pat = re.compile(r'[\.a-z]+')   # split at 'CGAGGCCGgtgag.....gccagGGTCTCACACTTGGCAGACGATGTATG'

    @staticmethod
    def _is_alignment_start(line):
        tokens = blank_split(line)
        return len(tokens) == 3 and tokens[1] == 'vs' and tokens[2][-1] == ':'

    @staticmethod
    def _is_new_block(line):
        return line.startswith('Note')

    @classmethod
    def parse(cls, fp):
        try:
            it = fp
            #it = dropwhile(lambda x: not _is_new_block(x), it)  # find a line starts with 'Note'
            line, it = skip_until(cls._is_new_block, it)  # find new block
            # body
            while 1:
                splice_lines, it = collect_while(make_not(isblank), it)  # collect exon, intron records
                next(it)  # blank
                span_line = next(it)  # blank
                next(it)  # blank
                segment_lines, it = collect_while(make_not(isblank), it)

                splicings = [blank_split(line.rstrip()) for line in splice_lines]
                segments = [blank_split(line.rstrip()) for line in segment_lines]
                span_row = blank_split(span_line)

                rname = span_row[5]
                qname = span_row[8]
                logging.debug('%s vs %s', rname, qname)

                # find alignment or next block
                line, it = skip_until(lambda x: cls._is_alignment_start(x) or cls._is_new_block(x), it)
                rseqs = qseqs = None
                if cls._is_alignment_start(line):
                    next(it)  # skip blank line
                    rbuf = []
                    qbuf = []
                    while 1:  # reading alignment
                        tokens = blank_split(next(it))  # ref line
                        if not tokens or tokens[0] != rname:  # end of alignment block
                            break
                        rbuf.append(tokens[2])
                        next(it)  # match line
                        tokens = blank_split(next(it))  # query line
                        qbuf.append(tokens[2])
                        next(it)  # blank line

                    rseqs = cls._intron_pat.split(''.join(rbuf))
                    qseqs = cls._intron_pat.split(''.join(qbuf))

                yield cls(span_row, splicings, segments, rseqs, qseqs)
                skip_until(cls._is_new_block, it)
        except StopIteration as e:
            pass

    def __init__(self, span_row, splice_rows, segment_rows, rseqs, qseqs):
        row = span_row
        if rseqs is not None:
            rseq = '|'.join(rseqs)
            qseq = '|'.join(qseqs)
        else:
            rseq = qseq = None
        self.span = self.Feature(int(row[1]), float(row[2]), int(row[3])-1, int(row[4]), row[5], int(row[6])-1, int(row[7]), row[8], rseq, qseq)
        self.exons = []
        self.introns = []
        exon_i = 0
        for row in splice_rows:
            if row[0] == 'Exon':
                if rseqs is not None:
                    rseq = rseqs[exon_i]
                    qseq = qseqs[exon_i]
                else:
                    rseq = qseq = None
                rec = self.Feature(int(row[1]), float(row[2]), int(row[3])-1, int(row[4]), row[5], int(row[6])-1, int(row[7]), row[8], rseq, qseq)
                self.exons.append(rec)
                exon_i += 1
            elif row[0][1:] == 'Intron':
                intron_type = row[0][0] # (+, -, ?)
                rec = self.Intron(intron_type, int(row[1]), float(row[2]), int(row[3])-1, int(row[4]), row[5])
                self.introns.append(rec)
            else:
                logging.error(row)
                raise NotImplementedError

        self.qname = self.span.qname
        self.rname = self.span.rname
        #self.segment_rows = segment_rows  # TODO



def read_vbseq(result, mean_rlen):
    """
    ID      LENGTH  Z       FPKM    THETA
    HLA:HLA00001    1047    0.00    0.0     0.000000e+00
    HLA:HLA00002    1047    0.00    0.0     0.000000e+00
    HLA:HLA00398    1047    0.00    0.0     0.000000e+00
    HLA:HLA00644    1047    0.00    0.0     0.000000e+00
    HLA:HLA00003    1047    45.74   7877.1781230    1.456456e-06
    """
    tab = pd.read_table(result)
    tab['mean_depth'] = 1. * mean_rlen * tab['Z'] / tab['LENGTH']
    return tab


@command.add_sub
@argument('vbseq')
@argument('-a', '--allele-list', required=True, help='IMGT HLA Allelelist')
@argument('--fasta-prefix', default='HLA:', help='HLA prefix of Allelelist')
@argument('-r', '--mean-rlen', required=True, type=float, help='mean single read length')
@argument('--is-paired', action='store_true', help='if set, twice the mean rlen for depth calculation')
@argument('-D', '--depth-threshold', type=float, default=5., help='threshold depth')
@argument('-g', '--genes', nargs='+')
@argument('-s', '--sampleid', default='NA')
def vbseq_hla_gt_org(args):
    """ Backward compatibility for the original definition
    """
    require(args.allele_list)
    logging.info('a prefix of refname is: %s', args.fasta_prefix)
    alist = _read_allele_list(args.allele_list, fasta_prefix=args.fasta_prefix)
    default_cn = 2  # fixed

    mean_rlen = args.mean_rlen
    if args.is_paired:
        mean_rlen *= 2
    tab = read_vbseq(args.vbseq, mean_rlen)
    tab = tab.merge(alist, left_on='ID', right_on='fasta_id')#, rsuffix='.1')
    if args.genes:
        gene_set = set(args.genes)
        tab = tab[tab.gene.map(lambda x: x in gene_set)]

    del tab['fasta_id']
    del tab['id']
    gene_depths = {}
    for gene, gtab in tab.groupby('gene'):
        gene_depths[gene] = gtab.mean_depth.sum()

    keys = '8d 6d 4d 2d'.split(' ')
    depths = {key: {} for key in keys}
    for key in keys:
        for gkey, gtab in tab.groupby('subtype_' + key):
            depths[key][gkey] = gtab.mean_depth.sum()

    tab['gene_cn'] = tab.gene.map(lambda x: default_cn)
    tab['gene_depth'] = tab.gene.map(lambda x: gene_depths[x])
    tab['ratio_in_gene'] = tab.mean_depth / tab.gene_depth
    for key in keys:
        d = tab['subtype_' + key].map(lambda x: depths[key][x])
        tab['ratio_' + key] = d / tab.gene_depth

    # add genotype
    def add_gt(tab, id_key='ID', ratio_key='ratio_in_gene', gt_key='gt', sel_key='sel'):
        def get_2copy(gtab):
            gtab = gtab.drop_duplicates(id_key, keep='first')
            depths = gtab[ratio_key] * gtab['gene_depth']
            depth_sorted = depths.sort_values(ascending=False)
            depth_tops = depth_sorted[depth_sorted >= args.depth_threshold]  # filter
            #print (depth_tops.index)
            #print (gtab.index)
            tops = gtab.loc[depth_tops.index]
            logging.debug(tops)
            if len(tops) >= 2:
                if depth_tops.iloc[0] >= depth_tops.iloc[1] * 2:
                    return ('hom', [tops.iloc[0][id_key]])
                else:
                    return ('het', [tops.iloc[0][id_key], tops.iloc[1][id_key]])
            elif len(tops) == 1:
                if depth_tops.iloc[0] >= args.depth_threshold * 2:
                    return ('hom', [tops.iloc[0][id_key]])
                else:
                    return ('het', [tops.iloc[0][id_key]])  # partly missing
            else:
                return ('unknown', [])

        gene_gts = {}
        selected_ids = set()
        for gene, gtab in tab.groupby('gene'):
            gt, ids = get_2copy(gtab)
            gene_gts[gene] = gt
            selected_ids.update(set(ids))

        tab[gt_key] = tab.gene.map(lambda x: gene_gts[x])
        tab[sel_key] = tab[id_key].map(lambda x: int(x in selected_ids))

    add_gt(tab)
    for key in keys:
        add_gt(tab, id_key='subtype_' + key, ratio_key='ratio_' + key, gt_key='gt_' + key, sel_key='sel_' + key)

    #print (tab)
    cond_show = (tab['mean_depth'] > 0)
    tab = tab[cond_show]
    tab = tab.sort_values('gene')
    tab = convert_hla_gt_format(tab, sampleid=args.sampleid)
    tab.to_csv('/dev/stdout', index=False, sep='\t')


@command.add_sub
@argument('vbseq')
@argument('-a', '--allele-list', required=True, help='IMGT HLA Allelelist')
@argument('--fasta-prefix', default='HLA:', help='HLA prefix of Allelelist')
@argument('-r', '--mean-rlen', required=True, type=float, help='mean single read length')
@argument('--is-paired', action='store_true', help='if set, twice the mean rlen for depth calculation')
@argument('--min-depth', type=float, default=0, help='min read depth (only used for output)')
@argument('--min-ratio', type=float, default=0, help='min read ratio (only used for output)')
@argument('--min-het-ratio', type=float, default=.25)
@argument('--min-hom-ratio', type=float, default=.75)
@argument('--hom-x', type=float, default=2.5)
@argument('-g', '--genes', nargs='+')
@argument('-c', '--copy-numbers', nargs='+', type=int)
@argument('-f', '--format', choices=['v0', 'v1'], default='v1', help='version of output format')
@argument('-s', '--sampleid', default='NA')
def vbseq_hla_gt(args):
    """
    """
    require(args.allele_list)
    logging.info('a prefix of refname is: %s', args.fasta_prefix)
    alist = _read_allele_list(args.allele_list, fasta_prefix=args.fasta_prefix)
    default_cn = 2
    gene_cns = {}
    if args.copy_numbers:
        assert len(args.genes) == len(args.copy_numbers)
        gene_cns = dict(zip(args.genes, args.copy_numbers))

    mean_rlen = args.mean_rlen
    if args.is_paired:
        mean_rlen *= 2
    tab = read_vbseq(args.vbseq, mean_rlen)
    tab = tab.merge(alist, left_on='ID', right_on='fasta_id')#, rsuffix='.1')
    if args.genes:
        gene_set = set(args.genes)
        tab = tab[tab.gene.map(lambda x: x in gene_set)]

    del tab['fasta_id']
    del tab['id']
    gene_depths = {}
    for gene, gtab in tab.groupby('gene'):
        gene_depths[gene] = gtab.mean_depth.sum()

    keys = '8d 6d 4d 2d'.split(' ')
    depths = {key: {} for key in keys}
    for key in keys:
        for gkey, gtab in tab.groupby('subtype_' + key):
            depths[key][gkey] = gtab.mean_depth.sum()

    tab['gene_cn'] = tab.gene.map(lambda x: gene_cns.get(x, default_cn))
    tab['gene_depth'] = tab.gene.map(lambda x: gene_depths[x])
    tab['ratio_in_gene'] = tab.mean_depth / tab.gene_depth
    for key in keys:
        d = tab['subtype_' + key].map(lambda x: depths[key][x])
        tab['ratio_' + key] = d / tab.gene_depth

    # add genotype
    def add_gt(tab, id_key='ID', ratio_key='ratio_in_gene', gt_key='gt', sel_key='sel'):
        def get_2copy(gtab):
            gtab = gtab.drop_duplicates(id_key, keep='first')
            tops = gtab[gtab[ratio_key] >= args.min_het_ratio].sort_values(ratio_key, ascending=False)[:2]
            if len(tops) >= 2:
                if tops.iloc[0][ratio_key] >= args.min_hom_ratio:
                    if tops.iloc[0][ratio_key] > tops.iloc[1][ratio_key] * args.hom_x:
                        return ('hom', [tops.iloc[0][id_key]])
                    else:
                        return ('het', [tops.iloc[0][id_key], tops.iloc[1][id_key]])
                else:
                    return ('het', [tops.iloc[0][id_key], tops.iloc[1][id_key]])
            elif len(tops) == 1:
                if tops.iloc[0][ratio_key] >= args.min_hom_ratio:
                    return ('hom', [tops.iloc[0][id_key]])
                else:
                    return ('het', [tops.iloc[0][id_key]])
            else:
                return ('unknown', [])

        def get_1copy(gtab):
            gtab = gtab.drop_duplicates(id_key, keep='first')
            tops = gtab[gtab[ratio_key] >= args.min_hom_ratio].sort_values(ratio_key, ascending=False)[:2]
            if len(tops) >= 1:
                return ('hom', [tops.iloc[0][id_key]])
            else:
                return ('unknown', [])

        gene_gts = {}
        selected_ids = set()
        for gene, gtab in tab.groupby('gene'):
            cn = gtab.iloc[0]['gene_cn']
            if cn == 2:
                gt, ids = get_2copy(gtab)
            elif cn == 1:
                gt, ids = get_1copy(gtab)
            elif cn == 0:
                gt, ids = 'hom', []
            else:
                raise NotImplemented
            gene_gts[gene] = gt
            selected_ids.update(set(ids))

        tab[gt_key] = tab.gene.map(lambda x: gene_gts[x])
        tab[sel_key] = tab[id_key].map(lambda x: int(x in selected_ids))

    add_gt(tab)
    for key in keys:
        add_gt(tab, id_key='subtype_' + key, ratio_key='ratio_' + key, gt_key='gt_' + key, sel_key='sel_' + key)

    #print (tab)
    cond_show = (tab['ratio_in_gene'] >= args.min_ratio) & (tab['mean_depth'] >= args.min_depth)
    tab = tab[cond_show]
    tab = tab.sort_values('gene')
    if args.format == 'v0':
        pass
    if args.format == 'v1':
        tab = convert_hla_gt_format(tab, sampleid=args.sampleid)
    tab.to_csv('/dev/stdout', index=False, sep='\t')


def convert_hla_gt_format(org_tab, sampleid='sampleid'):
    """ # TODO
    emit untyped record?
    emit copy number 0 record?
    """
    digit = 0
    tabs = []
    cols = 'sampleid gene digit subtype gene_cn allele_cn allele sel gene_depth ratio refnames'.split(' ')
    digits = [2, 4, 6, 8, 0]

    for gene, tab1 in org_tab.groupby('gene'):
        for digit in digits:
            if digit == 0:
                id_key = 'ID'
                ratio_key = 'ratio_in_gene'
                gt_key = 'gt'
                sel_key = 'sel'
            else:
                id_key = 'subtype_{0}d'.format(digit)
                ratio_key = 'ratio_{0}d'.format(digit)
                gt_key = 'gt_{0}d'.format(digit)
                sel_key = 'sel_{0}d'.format(digit)

            rtab = tab1.groupby(id_key).agg({'ID': lambda x: ','.join(x)}).rename(columns={'ID': 'refnames'}).reset_index()
            tab2 = tab1.drop_duplicates(id_key, keep='first')
            tab2 = tab2.merge(rtab)
            tab2 = tab2.sort_values(ratio_key, ascending=False).reset_index()
            recs = {col: [] for col in cols}
            for an, rec in enumerate(tab2.itertuples(), 1):
                recs['sampleid'].append(sampleid)
                recs['gene'].append(gene)
                recs['digit'].append(digit)
                recs['gene_cn'].append(rec.gene_cn)
                recs['gene_depth'].append(rec.gene_depth)
                recs['allele'].append(an)
                gt = getattr(rec, gt_key)
                if rec.gene_cn == 0:
                    acn = 0
                elif rec.gene_cn == 1 and gt == 'hom':
                    acn = 1
                elif rec.gene_cn == 2 and gt == 'hom':
                    acn = 2
                elif rec.gene_cn == 2 and gt == 'het':
                    acn = 1
                else:
                    acn = 0
                sel = getattr(rec, sel_key)

                recs['allele_cn'].append(acn * sel)
                recs['subtype'].append(getattr(rec, id_key))
                recs['sel'].append(sel)
                recs['ratio'].append(getattr(rec, ratio_key))
                recs['refnames'].append(rec.refnames)
            tab_out = pd.DataFrame(recs, columns=cols)
            tabs.append(tab_out)

    tab = pd.concat(tabs)[cols]
    return tab


@command.add_sub
@argument('hla_gt_table')
@argument('-s', '--sampleid')
@argument('--min-depth', type=float, default=0, help='min read depth (only used for output)')
@argument('--min-ratio', type=float, default=0, help='min read ratio (only used for output)')
def convert_hla_gt(args):
    """ Update format
    """
    tab = pd.read_table(args.hla_gt_table)
    tab = convert_hla_gt_format(tab, sampleid=args.sampleid)
    cond_show = (tab['ratio_in_gene'] >= args.min_ratio) & (tab['mean_depth'] >= args.min_depth)
    tab = tab[cond_show]
    tab.to_csv('/dev/stdout', index=False, sep='\t')


# pipeline experimental

class ValidationError(Exception):
    pass

class File(object):
    def __init__(self, name, check_exists=True, min_size=None):
        self.name = name
        self._validators = []  # {'fn', 'name'}
        if check_exists:
            self.add_validator(lambda self: self.exists(), name='exists')
        if min_size:
            self.add_validator(lambda self: self.get_size() >= min_size, name=fmt('check_size >= {min_size}'))

    def get_ctime(self):
        """ The last time the file's inode was changed
        Returns None if not exist
        """
        return self.exists() and os.path.getctime(self.name)

    def get_mtime(self):
        """ The last time the file's contents were changed
        Returns None if not exist
        """
        return self.exists() and os.path.getmtime(self.name)

    def get_size(self):
        return os.path.getsize(self.name)

    def exists(self):
        return os.path.exists(self.name)

    def add_validator(self, fn, name=None):
        v = {'fn': fn, 'name': name or fn.__name__}
        self._validators.append(v)

    def validate(self):
        """
        Raises: ValidationError if invalid
        """
        logging.info('Validation: %s', self.name)
        for v in self._validators:
            result = v['fn'](self)
            logging.info('- %s => %s', v['name'], result)
            if not result:
                raise ValidationError

    def get_line_count(self):  # TODO wc -l ?
        c = 0
        with open(self.name) as fp:
            for _ in fp:
                c += 1
        return c

    def __str__(self):
        return self.name  # as this is a file

def listize(item, list_types=(list, tuple)):
    """
    >>> listize(None)
    []
    >>> listize([1])
    [1]
    >>> listize((1, 2))
    (1, 2)
    >>> listize('abc')
    ['abc']
    """
    if item is None:
        return []
    if isinstance(item, list_types):
        return item
    return [item]


class Task(object):
    def __init__(self, task, name=None, **opts):
        if isinstance(task, string_types):
            self.name = name or 'sh'
            self._desc = task
            cmd = task
            task = lambda self: sh.call(cmd)
        else:
            self.name = name or task.__name__
            self._desc = task.__doc__
        assert callable(task)
        self._targets = []
        self._requires = []
        self._task = task

        targets = listize(opts.get('target'))
        for t in targets:
            self._append_target(t)

        requires = listize(opts.get('require'))  # TODO flatten ?
        for r in requires:
            self._append_require(r)

    def get_targets(self):
        return self._targets

    def get_requires(self):
        return self._requires

    def _append_target(self, target):
        if isinstance(target, string_types):
            target = File(target)
        assert isinstance(target, File)
        self._targets.append(target)

    def _append_require(self, require):  # TODO accept callable
        if isinstance(require, Task):
            for t in require._targets:
                self._append_require(t)
            return
        if isinstance(require, string_types):
            require = File(require)
        assert isinstance(require, File), fmt('{require}')
        self._requires.append(require)

    def execute(self):
        return self._task(self)

    def is_ready(self):
        """ Test if all requirements were satisfied
        """
        for item in self._requires:
            assert isinstance(item, File)  # TODO accept callable
            try:
                item.validate()
            except ValidationError:
                logging.info('* %s is not ready *', item)
                return False
        return True

    def is_fresh(self):
        """ Test if all targets have newer or equal mtime than that of requires
        """
        if not self._targets:
            return False   # cannot say fresh if target does not exist
        if not self._requires:
            return True    # always fresh
        req_max_mtime = max([x.get_mtime() or float('inf') for x in self._requires])
        logging.info('Latest req_max_mtime %s', req_max_mtime)
        return all(req_max_mtime <= (item.get_mtime() or 0) for item in self._targets)

    def is_done(self):
        """ Test if all targets were done (freshness of targets won't be tested)
        """
        for item in self._targets:
            assert isinstance(item, File)
            try:
                item.validate()
            except ValidationError:
                logging.info('* %s is not done *', item)
                return False
        return True

    def __str__(self):  # TODO
        head = self._desc and self._desc.split('\n', 1)[0]
        suffix = ' $ {}'.format(head) if head else ''
        return self.name + suffix


class OrderedPipeline(object):
    def __init__(self):
        self._tasks = []

    def append(self, task=None, name=None, target=None, require=None, **opts):
        """ Returns Task
        """
        if task is None:  # decorator
            return lambda task: self.append(task=task, name=name, target=target, require=require, **opts)
        task = Task(task, name=name, target=target, require=require, **opts)
        self._tasks.append(task)
        return task

    def _iter_tasks(self, steps=None):
        if steps is None:
            for step, task in enumerate(self._tasks, 1):
                step

    def run(self, steps=None, force_run=False, dry_run=False):
        if dry_run:
            logging.info('* dry_run mode *')
        if force_run:
            logging.info('* force_run mode *')
        for step, task in enumerate(self._tasks, 1):  # TODO iter_tasks
            logging.info('[Step %s] %s', step, task)
            for i, t in enumerate(task.get_targets(), 1):
                logging.info('Target  %s: %s (mtime: %s)', i, t, t.get_mtime())
            for i, r in enumerate(task.get_requires(), 1):
                logging.info('Require %s: %s (mtime: %s)', i, r, r.get_mtime())
            if not force_run and task.is_fresh() and task.is_done():
                logging.info('Target is fresh and done')
                continue
            logging.info('* Target should be executed *')
            if not task.is_ready():
                if not dry_run:
                    raise Exception(fmt('Task {task} is not ready'))
            if not dry_run:
                logging.info('Running task: %s', task)
                task.execute()

    @classmethod
    def as_recipe(cls, setup_fn):
        assert callable(setup_fn)
        self = cls()

        #TODO make option names changable
        @argument('--status', action='store_true')
        @argument('--dry-run', action='store_true')
        @argument('--force-run', action='store_true')
        @argument('--steps', help='e.g. $step_name 1 1- -4 1,2,5 1-5')
        @wraps(setup_fn, assigned=('__module__', '__name__'))
        def fn(args):
            """
            Pipeline operations
            --steps {step_names}
                $step_name
                $step_name1,$step_name2,$step_name3
                $step_name1-
                -$step_name2
                $step_name1-$step_name2

            default is automatic determination
            """
            setup_fn(self, args)
            self.run(steps=args.steps, force_run=args.force_run, dry_run=args.dry_run)

        # Concatenate docstring
        fn.__doc__ = '\n'.join((setup_fn.__doc__, fn.__doc__))
        return fn


as_recipe = OrderedPipeline.as_recipe

@command.add_sub
@argument('bam')
@argument('-a', '--allele-list', required=True, help='IMGT HLA Allelelist')
@argument('-o', '--outdir', default='bam2hla_results')
@as_recipe
def bam2hla(pipe, args):
    """
    bam file (assumed multimapped bam)

    Outputs:
        {outdir}/bam_depth_nm.txt
        {outdir}/bam_depth_nm_imgt.txt
    """
    fn_bam_depth_nm = File(fmt('{args.outdir}/bam_depth_nm.txt'), min_size=100)
    fn_bam_depth_nm_rc = File(fmt('{args.outdir}/bam_depth_nm_rc.txt'), min_size=100)
    fn_bam_depth_nm_imgt = File(fmt('{args.outdir}/bam_depth_nm_imgt.txt'), min_size=100)

    fn_bam_depth_nm_rc.add_validator(lambda self: self.get_line_count() == fn_bam_depth_nm.get_line_count(), name='has same lines')
    fn_bam_depth_nm_imgt.add_validator(lambda self: self.get_line_count() == fn_bam_depth_nm.get_line_count(), name='has same lines')

    pipe.append(fmt('mkdir -p {args.outdir}'), name='mkdir', target=fmt('{args.outdir}'))   # TODO mkdir shortcut

    @pipe.append(target=fn_bam_depth_nm, require=args.bam)
    def bam_depth_with_nm(self):
        args
        output = self.get_targets()[0].name # local variable is required
        galgo(fmt('bam_depth_with_nm {args.bam} --max-nm 10 --summary > {output}'))

    @pipe.append(target=fn_bam_depth_nm_rc, require=args.bam)
    def bam_depth_with_nm_rc(self):
        args
        output = self.get_targets()[0].name # local variable is required
        galgo(fmt('bam_depth_with_nm {args.bam} --max-nm 10 --read-count > {output}'))

    @pipe.append(target=fn_bam_depth_nm_imgt, require=[fn_bam_depth_nm, fn_bam_depth_nm_rc, args.allele_list])
    def join_table(self):
        input1 = self.get_requires()[0].name
        input2 = self.get_requires()[1].name
        output = self.get_targets()[0].name

        tab1 = pd.read_table(input1)
        tab2 = pd.read_table(input2).drop(['length'], axis=1)
        alist = _read_allele_list(args.allele_list)
        merged = tab1.merge(tab2, on='contig', suffixes=('', '.rc'))
        merged = alist.merge(merged, left_on='fasta_id', right_on='contig')
        merged.to_csv(output, index=False, sep='\t')
