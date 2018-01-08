#from __future__ import absolute_import
from __future__ import print_function
import logging
from .utils import collect_while, skip_until, isblank, make_not, blank_split
import sys
from six import string_types
import numpy as np
import pandas as pd
from .bioseq import Fasta


def decode_seq(seq, ref, same='-'):
    """
    >>> decode_seq('-G-.-A.-*-', 'ATACTTTG..')
    'AGA.TA.G*.'
    """
    return ''.join(r if a == same else a for a, r in zip(seq, ref))


def get_exon_frames(exon_seqs):
    """
    >>> get_exon_frames(['***GC*', 'GA****'])
    [(0, 0), (0, 0)]
    >>> get_exon_frames(['AT', 'GCTTTA', 'AT'])
    [(0, 2), (2, 2), (2, 1)]
    >>> get_exon_frames(['ATT'])
    [(0, 0)]
    """
    frames = []
    right = 0
    for seq in exon_seqs:
        left = right
        right = (left + len(seq)) % 3
        frames.append((left, right))
    return frames


class IMGTAlnData(object):
    """
    Similar to OrderedDict

    * : missing
    . : deletion
    """

    def __init__(self):
        self.subtypes = []
        self._sub_alns = {}
        self._part_lens = None  # e.g. (10, 25, 10, 5, 5)

    @property
    def part_lens(self):
        return self._part_lens

    def append(self, subtype, alns):
        """
        Args:
            subtype:
            alns: [aln]
        """
        assert subtype not in self._sub_alns
        alns = list(alns)  # copy here
        part_lens = tuple(map(len, alns))
        if not self.subtypes:
            self._part_lens = part_lens
        #logging.info(self._part_lens)
        #logging.info(part_lens)
        # if alns are truncated at last element, pad it
        if self._part_lens[-1] > part_lens[-1]:
            logging.warning('The last element of part_lens for subtype: %s is truncated', subtype)
            pad_len = self._part_lens[-1] - part_lens[-1]
            alns[-1] = alns[-1] + '*' * pad_len
            part_lens = tuple(map(len, alns))  # reculc
            logging.warning('%s missing bases were padded for subtype: %s', pad_len, subtype)

        try:
            assert part_lens == self._part_lens
        except AssertionError as e:
            logging.error('%s has inconsistent part_lens: %s => %s', subtype, self._part_lens, part_lens)
            logging.error('alns: %s', alns)
            logging.error('First subtype: %s', self.subtypes[0])
            logging.error('alns: %s', self._sub_alns[self.subtypes[0]])
            raise
        self.subtypes.append(subtype)
        self._sub_alns[subtype] = alns

    def get_extended(self, part_lens):
        assert len(self._part_lens) <= len(part_lens), 'Metric is not suitable for extension'
        assert self._part_lens <= part_lens, 'Metric is not suitable for extension'
        new = self.__class__()
        pad_len = sum(part_lens) - sum(self._part_lens)
        for sub in self.subtypes:
            alns = self._sub_alns[sub]
            aln = ''.join(alns) + '*' * pad_len
            alns = []
            s = 0
            for l in part_lens:
                alns.append(aln[s:s+l])
                s += l
            new.append(sub, alns)
        return new

    def __contains__(self, subtype):
        return subtype in self._sub_alns

    def __getitem__(self, subtype):
        return self._sub_alns[subtype]

    def __len__(self):
        return len(self.subtypes)

    def __iter__(self):
        return iter(self.subtypes)

    def items(self):
        for sub in self.subtypes:
            yield (sub, self._sub_alns[sub])

    def as_array(self):
        return np.array([list(''.join(alns)) for sub, alns in self.items()])

    @classmethod
    def load_fasta(cls, fname):
        with open(fname) as fp:
            fa = Fasta(fp)
            self = cls()
            for contig in fa:
                subtype = contig.name
                seq = contig.seq.replace('-', '.').replace('N', '*')
                self.append(subtype, [seq])
            return self

    @classmethod
    def load(cls, fname, sep='\t'):
        tab = pd.read_table(fname, sep=sep)
        self = cls()
        for rec in tab.itertuples():
            self.append(rec.subtype, rec.alignment.split('|'))
        return self

    @classmethod
    def from_array(cls, subtypes, part_lens, arr):
        assert arr.shape[0] == len(subtypes), ''
        assert arr.shape[1] == sum(part_lens), ''
        indexes = []
        s = 0
        for l in part_lens:
            indexes.append((s, s + l))
            s += l
        self = cls()
        for i, sub in enumerate(subtypes):
            alns = [''.join(arr[i, s:e]) for s, e in indexes]
            self.append(sub, alns)
        return self

    def save(self, file=sys.stdout, sep='\t'):
        if isinstance(file, string_types):
            with open(file) as fp:
                return self.save(file=fp, sep=sep)
        print ('subtype', 'alignment', sep=sep, file=file)
        for sub in self.subtypes:
            alns = self._sub_alns[sub]
            aln = '|'.join(alns)
            print (sub, aln, sep=sep, file=file)

    def save_as_fasta(self, file=sys.stdout, remove_padding=False):
        if isinstance(file, string_types):
            with open(file) as fp:
                return self.save_as_fasta(file=fp, remove_padding=remove_padding)

        for sub in self.subtypes:
            print ('>' + sub, file=file)
            alns = self._sub_alns[sub]
            aln = ''.join(alns)
            aln = aln.replace('*', 'N').replace('.', '-')
            if remove_padding:
                aln = aln.replace('-', '')
            print (aln, file=file)


class IMGTAlnFileBase(object):
    """
    * : missing
    . : deletion
    | : partition
    """

    def _parse_lines(self, aln_file):
        with open(aln_file) as fp:
            try:
                it = fp
                #headers, it = collect_while(make_not(isblank), it)
                isnot_blank = make_not(isblank)
                headers, it = collect_while(lambda l: isnot_blank(l) and not l.strip().startswith(self._head_token), it)
                bodies = []
                logging.info('Read headers')
                # body
                while 1:
                    data = {}
                    _, it = collect_while(isblank, it)  # skip blank lines

                    line = next(it).rstrip()
                    logging.debug(line)
                    tokens = blank_split(line)
                    if tokens[0] != self._head_token:
                        raise StopIteration
                    data['dna_poss'] = list(map(int, tokens[1:]))

                    line = next(it).rstrip()
                    logging.debug(line)
                    if line.strip().startswith('AA'):
                        tokens = blank_split(line)
                        assert tokens[0] == 'AA' and tokens[1] == 'codon'
                        data['aa_poss'] = list(map(int, tokens[2:]))
                        line = next(it).rstrip()
                        logging.debug(line)

                    data['offsets'] = [i for i, s in enumerate(line) if s == '|']

                    data['lines'], it = collect_while(make_not(isblank), it)
                    logging.debug('%s', data['lines'][0].rstrip())
                    if len(data['lines']) > 1:
                        logging.debug('%s', data['lines'][1].rstrip())
                    logging.debug('Read %s subtypes', len(data['lines']))
                    bodies.append(data)
            except StopIteration as e:
                pass
            logging.info('Read bodies')
        return headers, bodies

    def __init__(self, aln_file):
        headers, bodies = self._parse_lines(aln_file)
        # headers
        self.name = headers[0].rstrip()
        self.gene = self.name.split(' ')[0]  # HLA-A
        self.version = None
        if len(headers) > 2:
            self.version = headers[1].rstrip().split(':')[1].strip()  # 3.24.0.1
        # bodies
        offsets = []
        subtypes = subs = [blank_split(line)[0] for line in bodies[0]['lines']]
        self.data = IMGTAlnData()
        sub_alns = dict((subtype, []) for subtype in subs)  # {subtype: [exon alignment]}
        for body in bodies:
            offsets.extend(body['offsets'])

            # the first line of each body block is the refenence alignment
            line = body['lines'][0]
            tokens = blank_split(line)
            subtype = tokens[0]
            ref_aln = ''.join(tokens[1:])
            sub_alns[subtype].append(ref_aln)
            for line in body['lines'][1:]:
                tokens = blank_split(line)
                subtype = tokens[0]
                raw_aln = ''.join(tokens[1:])
                aln = decode_seq(raw_aln, ref_aln)
                sub_alns[subtype].append(aln)

        for subtype in subtypes:
            try:
                aln = ''.join(sub_alns[subtype])
                alns = aln.split('|')
                self.data.append(subtype, alns)
            except AssertionError:
                logging.warning('subtype: %s was ignored', alns)

        logging.info('Load %s subtypes', len(subs))
        logging.info('Alignment block length: %s', len(ref_aln))
        self._sub_alns = sub_alns

    @property
    def subtypes(self):
        return self.data.subtypes

    def __iter__(self):
        """
        Yields:
            (subtype, [aln_str])
        """
        for sub in self.subtypes:
            yield (sub, self.data[sub])

    def __len__(self):
        return len(self.subtypes)


class IMGTNucFile(IMGTAlnFileBase):
    '''
    Example:
        <header>

        HLA-A Nucleotide Sequence Alignments
        IMGT/HLA Release: 3.24.0.1
        Sequences Aligned: 2016 May 04
        Steven GE Marsh, Anthony Nolan Research Institute.
        Please see http://hla.alleles.org/terms.html for terms of use.
    '''
    _head_token = 'cDNA'


class IMGTGenFile(IMGTAlnFileBase):
    _head_token = 'gDNA'


def _gen_swap_exons(gen_alns, nuc_alns):
    """
    >>> _gen_swap_exons(['AAA', 'TTT', 'GGG'], ['TT..T..'])
    ['AAA', 'TT..T..', 'GGG']
    """
    assert len(gen_alns) == len(nuc_alns) * 2 + 1
    return [aln if i%2==0 else nuc_alns[(i-1)//2] for i, aln in enumerate(gen_alns)]


def _nuc_pad_missing(nuc_alns, gen_lens):
    """
    >>> _nuc_pad_missing(['AAAAAAA'], gen_lens=[3, 5, 4])
    ['***', 'AAAAAAA', '****']
    """
    #logging.info((len(nuc_alns), len(gen_lens)))
    assert len(gen_lens) == len(nuc_alns) * 2 + 1
    alns = []
    for i, aln in enumerate(nuc_alns):
        alns.append('*' * gen_lens[i * 2])
        alns.append(aln)
    alns.append('*' * gen_lens[-1])
    return alns


def merge_gen_nuc_alns(gen_data, nuc_data, gene_name=None):
    """
    Args:
        gen_data: IMGTAlnData
        nuc_data: IMGTAlnData
    Returns:
        IMGTAlnData
    """
    merged = IMGTAlnData()
    gen_lens = gen_data.part_lens
    if nuc_data.part_lens != gen_lens[1::2]:
        logging.warning('gen_data and nuc_data are inconsistence and need padding (%s, %s)', gen_lens, nuc_data.part_lens)
        nuc_data = nuc_data.get_extended(gen_lens[1::2])  # exons
    logging.info(gen_data.part_lens)
    logging.info(nuc_data.part_lens)
    logging.info('Gene name restriction: %s', gene_name)

    for sub, nuc_alns in nuc_data.items():
        if gene_name is not None:
            sub_gene = sub.split('*')[0]
            if sub_gene != gene_name:
                logging.info('skip gene %s', sub_gene)
                continue
        if sub in gen_data:
            gen_alns = gen_data[sub]
            alns = _gen_swap_exons(gen_alns, nuc_alns)
            merged.append(sub, alns)
        else:
            alns = _nuc_pad_missing(nuc_alns, gen_lens)
            merged.append(sub, alns)

    return merged
