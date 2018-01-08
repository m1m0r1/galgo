from __future__ import absolute_import
from __future__ import print_function
from builtins import map
import sys
import logging
from argtools import command, argument
from itertools import groupby
from collections import defaultdict, namedtuple
from .. import sh
import re
import numpy as np
import pandas as pd
import pysam
from ..utils import memoize, iter_tabs, Counter
from ..io import open
import tqdm


@command.add_sub
@argument('vcf', help='VCF or gzipped VCF')
def vcf_header(args):
    vcffile = '/dev/stdin' if args.vcf == '-' else args.vcf
    reader = pysam.VariantFile(vcffile, ignore_truncation=True)
    h = reader.header
    print ('category', 'key', 'number', 'type', 'description', sep='\t')
    for (key, rec) in h.filters.items():
        print ('FILTER', rec.id, '', '', rec.description, sep='\t')

    for (key, rec) in h.info.items():
        print ('INFO', rec.id, rec.number, rec.type, rec.description, sep='\t')

    for (key, rec) in h.formats.items():
        print ('FORMAT', rec.id, rec.number, rec.type, rec.description, sep='\t')

    for (key, rec) in h.alts.items():
        print ('ALT', key, '', '', rec.get('Description', '.'), sep='\t')

    for (key, rec) in h.metadata.items():
        print ('METADATA', key, '', '', rec, sep='\t')


@command.add_sub
@argument('vcf', help='VCF or gzipped VCF')
def vcf_samples(args):
    vcffile = '/dev/stdin' if args.vcf == '-' else args.vcf
    reader = pysam.VariantFile(vcffile, ignore_truncation=True)
    for sample in reader.header.samples:
        print (sample, sep='\t')


def _get_sample_flags(reader, args):
    sample_set = set()
    if args.samples:
        sample_set = set(args.samples.split(','))
    if args.samples_file:
        sample_set = set([line.rstrip() for line in open(args.samples_file)])

    if sample_set:
        sample_flags = np.array([s in sample_set for s in reader.header.samples], dtype=bool)
    else:
        sample_flags = np.ones(len(reader.header.samples), dtype=bool)

    return sample_flags


class VCFRefFix(object):  # for VCF converted from Plink
    """ TODO
    - plink2chr
    Chrom name mapping
    Omni manifest file is required.
    """
    def __init__(self, reader, ref_fasta):
        self._reader = reader
        it = reader
        fasta = pysam.FastaFile(ref_fasta)
        for chrom, recs in groupby(it, lambda rec: rec.chrom):
            chrom_seq = fasta.fetch(reference=chrom)
            fixed_count = 0
            for rec in recs:
                try:
                    is_fixed = self._fix_rec(rec, chrom_seq)
                    fixed_count += is_fixed
                    yield rec
                except Exception as e:
                    logging.error('Cannot handle indels %s: %s', rec, e)
                    yield rec
            logging.error('Cannot handle %s: %s', rec, e)

    def _fix_rec(self, rec, chrom_seq):
        pos0 = rec.pos - 1
        assert len(rec.ref) == 1 and all(len(a) for a in rec.alts)
        true_ref = chrom_seq[pos0:pos0 + refl]
        if true_ref == rec.ref:  # No need to fix
            return 0
        # You need to fix ref and alt


@command.add_sub
@argument('vcf', help='VCF or gzipped VCF')
@argument('-f', '--ref-fasta', required=True, help='')
@argument('-r', '--region')
def vcf_ref_check(args):
    vcffile = '/dev/stdin' if args.vcf == '-' else args.vcf
    reader = pysam.VariantFile(vcffile, ignore_truncation=True)
    if args.region:
        it = reader.fetch(region=args.region)
    else:
        it = reader
    fasta = pysam.FastaFile(args.ref_fasta)
    for chrom, recs in groupby(it, lambda rec: rec.chrom):
        chrom_seq = fasta.fetch(reference=chrom)
        for rec in tqdm.tqdm(recs):
            start = rec.pos-1
            rlen = len(rec.ref)
            end = start + rlen
            seq = chrom_seq[start:end]
            if rec.ref != seq:
                logging.warning('REF is inconsistent: %s != %s at %s:%s', rec.ref, seq, chrom, rec.pos)


class VCFChromMap(object):
    def __init__(self, reader, chrom_map):
        self._chrom_map = chrom_map
        self._reader = reader
        header_rows = str(reader.header).split('\n')

        ##contig=<ID=1,length=249212726>  # length is falsy for plink VCF
        pat = re.compile('##contig=<ID=([\w]+),.*')
        def replacer(match):
            name = match.group(1)
            line = match.group(0)
            name1 = self._chrom_map.get(name, name)
            return line[:match.start(1)] + name1 + line[match.end(1):]

        for i, line in enumerate(header_rows):
            if line.startswith('##contig'):
                line = pat.sub(replacer, line)
                header_rows[i] = line

        self.header = '\n'.join(header_rows) # TODO currently no way to rewrite header by pysam api

    def __iter__(self):
        for chrom, recs in groupby(self._reader, lambda rec: rec.chrom):
            chrom1 = self._chrom_map.get(chrom, chrom)
            for rec in recs:
                # TODO cannot replace rec.chrom that is not defined in the header which is immutable
                # We need more flexible library
                rec = str(rec).split('\t')
                rec[0] = chrom1
                yield '\t'.join(rec)


def get_chrom_map(mode=None, maps=None):
    chrom_map = {}
    n_autosomal = 22
    num_chr_pairs = [(str(n), 'chr' + str(n)) for n in range(1, n_autosomal+1)] + [('X', 'chrX'), ('Y', 'chrY'), ('XY', 'chrXY'), ('MT', 'chrM')]
    plink_num_pairs = [(str(n), str(n)) for n in range(1, n_autosomal+1)] + [('23', 'X'), ('24', 'Y'), ('25', 'XY'), ('26', 'MT')]  # 25 is PAR
    num_chrs = dict(num_chr_pairs)
    if mode is None:
        pass
    elif mode == 'chr2num':
        chrom_map.update((c, n) for n, c in num_chr_pairs)
    elif mode == 'num2chr':
        chrom_map.update(num_chrs)
    elif mode == 'plink2num':
        chrom_map.update(plink_num_pairs)
    elif mode == 'plink2chr':
        chrom_map.update((p, num_chrs[n]) for p, n in plink_num_pairs)
    else:
        logging.error('Unknown mode: %s', mode)
        raise NotImplementedError()
    if maps:
        chrom_map.update(maps)
    return chrom_map


# TODO move to bed module
@command.add_sub
@argument('bed', help='Bed or gzipped VCF')
@argument('--header', action='store_true')
@argument('-c', '--chrom-mode', choices=['chr2num', 'num2chr', 'plink2num', 'plink2chr'])
@argument('-m', '--chrom-maps', nargs='+', help='before,after')
def bed_chrom_map(args):
    """
    the bed file should be separated by tab
    """
    infile = '/dev/stdin' if args.bed == '-' else args.bed
    chrom_map = get_chrom_map(mode=args.chrom_mode, maps=args.chrom_maps)
    with open(infile) as fp:
        if args.header:
            line = next(fp)
            print (line, end='')

        for line in fp:
            if line.startswith('#'):
                print (line, end='')
            else:
                rec = line.rstrip('\n\r').split('\t')
                chrom = rec[0]
                chrom = chrom_map.get(chrom, chrom)
                print (chrom, *rec[1:], sep='\t')


@command.add_sub
@argument('vcf', help='VCF or gzipped VCF')
@argument('-c', '--chrom-mode', choices=['chr2num', 'num2chr', 'plink2num', 'plink2chr'])
@argument('-m', '--chrom-maps', nargs='+', help='before,after')
def vcf_chrom_map(args):
    vcffile = '/dev/stdin' if args.vcf == '-' else args.vcf
    chrom_map = get_chrom_map(mode=args.chrom_mode, maps=args.chrom_maps)
    reader = pysam.VariantFile(vcffile, ignore_truncation=True)
    reader = VCFChromMap(reader, chrom_map)

    print (reader.header, end='')
    for rec in reader:
        print (rec, end='')


@memoize
def get_gt_hom_idxs(alt_num):
    """
    >>> get_gt_hom_idxs(0)  # 0/0
    [0]
    >>> get_gt_hom_idxs(1)  # 0/0, 0/1, 1/1
    [0, 2]
    >>> get_gt_hom_idxs(2)  # 0/0, 0/1, 1/1, 0/2, 1/2, 2/2
    [0, 2, 5]
    >>> get_gt_hom_idxs(3)
    [0, 2, 5, 9]
    """
    last = -1
    hom_idxs = []
    for a in range(alt_num + 1):
        last = last + (a + 1)
        hom_idxs.append(last)
    return hom_idxs


@memoize
def get_default_gl(alt_num):
    """
    >>> get_default_gl(1)
    '-0.301,0,-0.301'
    >>> get_default_gl(2)
    '-0.301,0,-0.301,0,0,-0.301'
    """
    het_gl = '0'
    hom_gl = '{0:.3f}'.format(- np.log10(2))
    gt_num = (alt_num + 1) * (alt_num + 2) // 2
    gl = [het_gl] * gt_num
    hom_idxs = get_gt_hom_idxs(alt_num)
    for idx in hom_idxs:
        gl[idx] = hom_gl
    return ','.join(gl)


@command.add_sub
@argument('vcf', help='VCF or gzipped VCF')
@argument.exclusive(
    argument('-s', '--samples', help='comma separated sample list'),
    argument('-S', '--samples-file', help='sample list file'),
)
def vcf_reset_gl(args):
    """ whitening format tag for specified samples

    # TODO also reset PL
    """
    vcffile = '/dev/stdin' if args.vcf == '-' else args.vcf
    import pysam
    reader = pysam.VariantFile(vcffile)
    sample_flags = _get_sample_flags(reader, args)
    print (reader.header, end='')
    sample_offset = 9

    # hom: log10 p**2 = 2 log10 p
    # het: log10 2p**2 = 2 log10 p + log10 2

    def replace_fmt(alt_num, gl_index, fmt_info):
        splt = fmt_info.split(':')
        splt[gl_index] = get_default_gl(gt_num)
        return ':'.join(splt)

    for rec in reader:
        srec = str(rec).split('\t')
        fmts = rec.format.keys()
        alt_num = len(rec.alts)
        gl_index = fmts.index('GL')

        fmt_rec = [(replace_fmt(alt_num, gl_index, s) if flg else s) for flg, s in zip(sample_flags, srec[sample_offset:])]
        print (*(srec[:sample_offset] + fmt_rec), sep='\t', end='')


@command.add_sub
@argument('vcf', help='VCF or gzipped VCF')
@argument.exclusive(
    argument('-s', '--samples', help='comma separated sample list'),
    argument('-S', '--samples-file', help='sample list file'),
)
def vcf_untype(args):
    """ whitening format tag for specified samples

    # TODO also reset PL
    """
    vcffile = '/dev/stdin' if args.vcf == '-' else args.vcf
    import pysam
    reader = pysam.VariantFile(vcffile)
    sample_flags = _get_sample_flags(reader, args)
    print (reader.header, end='')
    sample_offset = 9

    import re
    gt_sep = re.compile('/|')
    def replace_fmt(gt_index, fmt_info):
        splt = fmt_info.split(':')
        gt = splt[gt_index]
        ngt = len(gt_sep.split(gt))
        splt[gt_index] = '/'.join(('.',) * ngt)
        return ':'.join(splt)

    for rec in reader:
        srec = str(rec).split('\t')
        fmts = rec.format.keys()
        gt_index = fmts.index('GT')

        fmt_rec = [(replace_fmt(gt_index, s) if flg else s) for flg, s in zip(sample_flags, srec[sample_offset:])]
        print (*(srec[:sample_offset] + fmt_rec), sep='\t', end='')


def _load_pos_sets(poslist):
    chrom_pos_sets = defaultdict(set)
    with open(poslist) as fp:
        for row in iter_tabs(fp):
            chrom = row[0]
            pos = int(row[1])
            chrom_pos_sets[chrom].add(pos)
    return chrom_pos_sets


@command.add_sub
@argument('vcf', help='VCF or gzipped VCF')
@argument('poslist', help='chrom:pos')
@argument('-r', '--region', help='VCF or gzipped VCF')
def vcf_pick(args):
    """ poslist: chrom:pos
    """
    vcffile = '/dev/stdin' if args.vcf == '-' else args.vcf
    import pysam
    reader = pysam.VariantFile(vcffile)

    # iteration
    if args.region:
        it = reader.fetch(region=args.region)
    else:
        it = reader

    logging.info('Loading poslist: %s', args.poslist)
    chrom_pos_sets = _load_pos_sets(args.poslist)
    for chrom in sorted(chrom_pos_sets):
        logging.info('Loaded %s pos from chrom: %s', len(chrom_pos_sets[chrom]), chrom)

    print (reader.header, end='')
    for chrom, recs in groupby(it, lambda x: x.contig):
        pos_sets = chrom_pos_sets[chrom]
        for rec in recs:
            if rec.pos in pos_sets:
                print (str(rec), end='')


@command.add_sub
@argument('vcf', help='VCF or gzipped VCF')
@argument('-r', '--region')
@argument('-c', '--check-variant', action='store_true', help='consider duplicated if check ref and alts are concordant')
def vcf_list_dup(args):
    """
    Emit:
        chrom, pos, id, ref, alt
    """
    vcffile = '/dev/stdin' if args.vcf == '-' else args.vcf
    import pysam
    reader = pysam.VariantFile(vcffile)

    # iteration
    if args.region:
        it = reader.fetch(region=args.region)
    else:
        it = reader

    if args.check_variant:
        def get_id(rec):
            return ':'.join((rec.chrom, str(rec.pos), rec.ref, ','.join(rec.alts)))
    else:
        def get_id(rec):
            return ':'.join((rec.chrom, str(rec.pos)))

    print ('chrom', 'pos', 'id', 'ref', 'alt', sep='\t')
    checked = set()  # pos ':' ref ':' ','.join(alts)
    for chrom, recs in groupby(it, lambda x: x.contig):
        last_pos = 0
        for rec in recs:
            rec_id = get_id(rec)
            if rec_id in checked:
                print (rec.chrom, rec.pos, rec.id, rec.ref, ','.join(rec.alts), sep='\t')
                continue
            checked.add(rec_id)


@memoize
def mendel_violate(child, father, mother):
    """
    >>> mendel_violate((0, 0), (0, 0), (0, 0))
    0
    >>> mendel_violate((0,), (0,), (0, 0))
    0
    >>> mendel_violate((1,), (0,), (0, 0))
    1
    >>> mendel_violate((1,), (0,), (1, 2))
    0
    >>> mendel_violate((2, 0), (0,), (1, 2))
    0
    >>> mendel_violate((1, 2), (0,), (1, 2))
    1
    >>> mendel_violate((0, 1), (1, 1), (0, 0))
    0
    >>> mendel_violate((0, 1), (1, 1), (0, 1))
    0
    >>> mendel_violate((1, 1), (1, 1), (0, 1))
    0
    >>> mendel_violate((1,), (1,), ())
    0
    >>> mendel_violate((1,), (0,), ())
    1
    >>> mendel_violate((1,), (None,), ())
    0
    >>> mendel_violate((0, 1), (None,), (None, None))
    0
    >>> mendel_violate((0, 0), (None,), (1, 1))
    1
    >>> mendel_violate((None, None), (1,), (0, 1))
    0
    >>> mendel_violate((0, 1), (None, None), (0, 0))
    0
    >>> mendel_violate((0, 0), (0, 0), (None,))
    0
    """
    assert len(child) <= 2
    assert len(mother) <= 2
    assert len(father) <= 2
    assert child == (None,) * len(child) or None not in child
    assert mother == (None,) * len(mother) or None not in mother
    assert father == (None,) * len(father) or None not in father
    if len(child) == 0:
        return 0
    if len(child) == 1:
        if None in child:
            return 0
        if None in mother:
            return 0
        if None in father:
            return 0
        if child[0] in mother or child[0] in father:
            return 0
    # child is diploid
    if None in child:
        return 0
    if None in mother and None in father:
        return 0
    if None in mother:
        if child[0] in father or child[1] in father:
            return 0
    if None in father:
        if child[0] in mother or child[1] in mother:
            return 0
    if child[0] in mother and child[1] in father:
        return 0
    if child[0] in father and child[1] in mother:
        return 0
    return 1


class PedFile(object):
    """
    Family ID
    Individual ID
    Paternal ID
    Maternal ID
    Sex (1=male; 2=female; other=unknown)
    Phenotype ...
    """
    def __init__(self, pedfile, header=None, sep=' '):
        columns = ['family', 'sample_id', 'father', 'mother', 'sex']
        self._ped_tab = pd.read_table(pedfile, header=header, sep=sep,
                #dtype={'family': str, 'sample_id': str, 'father': str, 'mother': str, 'sex': int})
                dtype={'x01': str, 'x02': str, 'x03': str, 'x04': str, 'x05': int})
        self._ped_tab.columns = columns + list(self._ped_tab.columns[len(columns):])
        self._recs = {t.sample_id: t for t in self._ped_tab.itertuples()}

    def __len__(self):
        return len(self._ped_tab)

    def get_sex(self, sample):
        return self._recs[sample].sex

    def iter_trios(self):
       """
       Returns (sample_id, father, mother)
       """
       for row in self._ped_tab.itertuples():
           if row.father == '0':
               continue
           if row.mother == '0':
               continue
           yield (row.sample_id, row.father, row.mother)

    def iter_male_duos(self):
        return self.iter_duos(father=True, mother=False, sex=1)

    def iter_duos(self, father=True, mother=True, sex=None):
        """
        Returns (sample_id, parent_id)
        """
        for row in self._ped_tab.itertuples():
            if sex is None or row.sex == sex:
                if father and row.father != '0':
                    yield (row.sample_id, row.father)
                if mother and row.mother != '0':
                    yield (row.sample_id, row.mother)



def _check_is_snp(ref, alts):
    return len(ref) == 1 and all(len(alt) == 1 for alt in alts)

class VCFSampleData(object):
    #logging.info((dir(sample_data),
    #    sample_data['GT'],
    #    str(sample_data['GT']),
    #    str(sample_data),
    #    sample_data.items(),
    #    sample_data.index,
    #    sample_data.alleles,
    #    sample_data.allele_indices,
    #    sample_data.phased))

    def __init__(self, d):
        self._d = d
        self.gt = self._d['GT']

    #@property
    #def gt_filled(self):
    #    return tuple('.' if index is None else index for index in self.gt)

    @property
    def gt_str(self):
        sep = '|' if self._d.phased else '/'
        return sep.join('.' if index is None else str(index) for index in self.gt)


MENDEL_CHECK_MODES = ('trio', 'male_duo')

class VCFMendelCheck(object):
    var_key = ['chrom', 'pos', 'end', 'ref', 'alt', 'is_snp', 'al_count']
    pattern_key = ['invalid', 'gt_c', 'gt_f', 'gt_m', 'sex']
    sample_key = ['trio']

    Var = namedtuple('Var', var_key)
    Pattern = namedtuple('Pattern', pattern_key)

    def __init__(self, vcf, pedfile, regions=None, mode='trio'):
        # VCF
        self._ped = ped = pedfile
        trios = list(ped.iter_trios())
        male_duos = list(ped.iter_male_duos())
        self.mode = mode
        assert mode in MENDEL_CHECK_MODES
        if regions:
            reader = pysam.VariantFile(vcf)
            def gen():
                for region in regions:
                    for rec in reader.fetch(region=region):
                        yield rec
            self._it = gen()
        else:
            fp = open(vcf)
            reader = pysam.VariantFile(fp)
            self._it = reader
        self.sample_idxs = {sample: i for i, sample in enumerate(reader.header.samples)}
        sample_set = set(reader.header.samples)
        #logging.info(sample_set)
        #logging.info(trios)
        #logging.info(male_duos)
        self.trio_samples = [trio for trio in trios if all(sid in sample_set for sid in trio)]
        self.male_duo_samples = [male_duo for male_duo in male_duos if all(sid in sample_set for sid in male_duo)]
        logging.debug(sample_set)

    def __iter__(self):
        sample_idxs = self.sample_idxs
        Var = self.Var
        Pattern = self.Pattern
        ped = self._ped
        for row in self._it:
            info = row.info
            ref = row.ref
            alts = row.alts
            end = row.pos + len(ref) - 1
            is_snp = int(_check_is_snp(ref, alts))
            al_count = 1 + len(alts)
            var = Var(row.chrom, row.pos, info.get('END', end), ref, ','.join(alts), is_snp, al_count)
            if self.mode == 'trio':
                for c, f, m in self.trio_samples:
                    sex = ped.get_sex(c)
                    c_idx = sample_idxs[c]
                    f_idx = sample_idxs[f]
                    m_idx = sample_idxs[m]
                    #logging.info(dir(row.format['GT']))
                    dat_c = VCFSampleData(row.samples[c_idx])
                    dat_f = VCFSampleData(row.samples[f_idx])
                    dat_m = VCFSampleData(row.samples[m_idx])
                    invalid = mendel_violate(dat_c.gt, dat_f.gt, dat_m.gt)
                    pattern = Pattern(invalid, dat_c.gt_str, dat_f.gt_str, dat_m.gt_str, sex)
                    yield var, pattern, (c, f, m)
            elif self.mode == 'male_duo':
                for c, f in self.male_duo_samples:
                    c_idx = sample_idxs[c]
                    f_idx = sample_idxs[f]
                    #logging.info(dir(row.format['GT']))
                    dat_c = VCFSampleData(row.samples[c_idx])
                    dat_f = VCFSampleData(row.samples[f_idx])
                    invalid = mendel_violate(dat_c.gt, dat_f.gt, ())
                    pattern = Pattern(invalid, dat_c.gt_str, dat_f.gt_str, '.', 1)
                    yield var, pattern, (c, f)


import codecs
def unescaped_str(s):
    return codecs.decode(str(s), 'unicode_escape')

@command.add_sub
@argument('vcf', help='VCF or gzipped VCF')
@argument('-p', '--pedigree', help='ped file')
@argument('-s', '--ped-sep', default=' ', type=unescaped_str)
@argument('-r', '--regions', nargs='+')
@argument('--stats', action='store_true', help='Only emit statistics for each variant')
@argument('-m', '--mode', choices=MENDEL_CHECK_MODES, default='trio', help='Use male_duo for check chrY consistency')
def vcf_mendel_check(args):
    """
    """
    pedfile = PedFile(args.pedigree, sep=args.ped_sep)
    checker = VCFMendelCheck(args.vcf, pedfile, regions=args.regions, mode=args.mode)
    if args.mode == 'trio':
        logging.info('Trios in PED: %s (%s)', len(list(pedfile.iter_trios())), args.pedigree)
        logging.info('Trios in the VCF: %s (%s)', len(checker.trio_samples), args.vcf)
    elif args.mode == 'male_duo':
        logging.info('Male duos in PED: %s (%s)', len(list(pedfile.iter_male_duos())), args.pedigree)
        logging.info('Male duos in the VCF: %s (%s)', len(checker.male_duo_samples), args.vcf)
    else:
        raise NotImplementedError

    mode = args.mode
    if args.stats:
        keys = checker.var_key + checker.pattern_key + ['count']
        print (*keys, sep='\t')
        for v, recs in groupby(checker, (lambda v, p, t: v)):
            patterns = Counter(p for _, p, _ in recs)
            for pat, count in patterns.items():
                print (*(v + pat + (count,)), sep='\t')
    else:
        keys = checker.var_key + checker.pattern_key + [mode]
        print (*keys, sep='\t')
        for v, pat, members in checker:
            mem = ','.join(members)
            print (*(v + pat + (mem,)), sep='\t')



@command.add_sub
@argument('mendel_tab')
def vcf_mendel_summary(args):
    """
    """
    group_key = ('chrom', 'is_snp', 'al_count', 'sex', 'gt_f', 'gt_m')
    level = list(range(len(group_key)))
    tab = pd.read_table(open(args.mendel_tab),
            dtype={'gt_c': 'str', 'gt_f': 'str', 'gt_m': 'str'})
    gtab = tab.groupby(group_key).apply(lambda tab: pd.DataFrame({
        'count': pd.Series([len(tab)]),
        'invalid': pd.Series([tab.invalid.sum()]),
        'rate': pd.Series([1. * tab.invalid.sum() / len(tab)]),
        },
        columns=['count', 'invalid', 'rate'])
    ).reset_index(level=level)
    gtab.to_csv(sys.stdout, sep='\t', index=False)
