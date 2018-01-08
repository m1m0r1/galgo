from __future__ import print_function
from argtools import command, argument
from collections import defaultdict
import logging
import os
import tempfile
import pysam
import numpy as np
import pandas as pd
from statsmodels import robust
from ..utils import iter_tabs
from ..samutil import Read
from .. import plot
from .. import sh


class SamMethods(object):
    def __init__(self, sam):
        self._sam = sam
        self.filename = sam.filename
        index1 = '{0}.bai'.format(sam.filename)
        index2 = '{0}.bai'.format(sam.filename.rsplit('.', 2)[0])
        #if os.path.exists(index1):  # no use
        #    self.indexpath = index1
        #elif os.path.exists(index2):
        #    self.indexpath = index2
        #else:
        #    self.indexpath = None

    def _gen_it(self, chroms=None, start=0):
        if chroms is None and not start:
            for rec in self._sam:
                yield rec
        else:
            if not chroms:
                chroms = self._sam.references
            for chrom in chroms:
                sam_it = self._sam.fetch(reference=chrom, start=start)
                for rec in sam_it:
                    yield rec

    def _get_stat(self, vals, stat_fn='mean'):
        if isinstance(stat_fn, (tuple, list)):
            return [self._get_stat(vals, fn) for fn in stat_fn]

        if stat_fn is None or stat_fn == 'mean':
            return np.mean(vals)
        if stat_fn == 'median':
            return np.median(vals)
        if stat_fn == 'mad':
            return robust.mad(vals)
        logging.warning(stat_fn)
        raise NotImplementedError

    def estimate_rlen(self, chroms=None, rsample=1000, use_alen=False, skip_flag=0x900, start=0, stat_fn=None):
        if use_alen:
            fn = lambda rec: rec.alen
        else:
            fn = lambda rec: len(rec.seq)
        vals = self._estimate_fn(fn=fn, chroms=chroms, rsample=rsample, skip_flag=skip_flag, start=start)
        return self._get_stat(vals, stat_fn=stat_fn)

    def estimate_inslen(self, chroms=None, rsample=1000, use_alen=False, skip_flag=0x904 + 0x008, start=0, stat_fn=None):
        fn = lambda rec: abs(rec.template_length)
        test_fn = lambda rec: rec.reference_id == rec.next_reference_id
        vals = self._estimate_fn(fn=fn, chroms=chroms, rsample=rsample, skip_flag=skip_flag, start=start,
                test_fn=test_fn)
        return self._get_stat(vals, stat_fn=stat_fn)

    def _estimate_fn(self, fn, chroms=None, rsample=1000, skip_flag=0x900, start=None, test_fn=None):
        """ Rough estimattion
        """
        vals = []
        sam_it = self._gen_it(chroms=chroms, start=start)
        try:
            i = rsample
            while i:
                rec = next(sam_it)
                if test_fn and not test_fn(rec):
                    continue
                if not (rec.flag & skip_flag):
                    vals.append(fn(rec))
                    i -= 1

        except StopIteration:
            pass
        self._sam.reset()
        return vals

    def get_reference_length(self, chroms=None):
        if chroms is None:
            return sum(self._sam.lengths) - nlen
        chrom_set = set(chroms)
        tlen = sum(length for ref, length in zip(self._sam.references, self._sam.lengths) if ref in chrom_set)
        return tlen

    def _get_index_stats(self):
        with tempfile.NamedTemporaryFile() as fp:
            logging.info('Emit index_stats to %s', fp.name)
            sh.call(['samtools', 'idxstats', self.filename], stdout=fp.name)
            tab = pd.read_table(fp, header=None, names=['reference', 'length', 'mapped', 'unmapped'])
        return tab

    def get_mapped(self, chroms=None):
        if chroms is None:
            return self._sam.mapped
        stats = self._get_index_stats()
        chrom_set = set(chroms)
        maps = [mapped for chrom, mapped in zip(stats.reference, stats.mapped) if chrom in chrom_set]
        return sum(maps)


class Nbed(object):
    def __init__(self, nbed):
        self._tab = pd.read_table(nbed, header=None, names=['contig', 'start', 'end'])

    def get_length(self, contig):
        return (self._tab.end - self._tab.start)[self._tab.contig == contig].sum()


# deprecated
@command.add_sub
@argument('bamfile')
@argument('-r', '--rsample', type=int, default=1000)
@argument('-s', '--start', type=int, default=0)
def bam_mean_rlen(args):
    sam = pysam.Samfile(args.bamfile, 'rb')
    samm = SamMethods(sam)
    rlen = samm.estimate_rlen(rsample=args.rsample, start=args.start, stat_fn='mean')
    print (rlen)


@command.add_sub
@argument('bamfile')
@argument('-r', '--rsample', type=int, default=1000)
@argument('-s', '--start', type=int, default=0)
@argument('-f', '--stat-fn', choices=['mean', 'median', 'mad'])
def bam_rlen(args):
    sam = pysam.Samfile(args.bamfile, 'rb')
    samm = SamMethods(sam)
    rlen = samm.estimate_rlen(rsample=args.rsample, start=args.start, stat_fn=args.stat_fn)
    print (rlen)


@command.add_sub
@argument('bamfile')
@argument('-r', '--rsample', type=int, default=1000)
@argument('-s', '--start', type=int, default=0)
@argument('-f', '--stat-fn', choices=['mean', 'median', 'mad'])
def bam_inslen(args):
    sam = pysam.Samfile(args.bamfile, 'rb')
    samm = SamMethods(sam)
    inslen = samm.estimate_inslen(rsample=args.rsample, start=args.start, stat_fn=args.stat_fn)
    print (inslen)


@command.add_sub
@argument('bamfile')
def bam_mean_cov(args):
    """ Emit mean coverage (#reads) per base
    """
    sam = pysam.Samfile(args.bamfile, 'rb')
    totlen = sum(sam.lengths)    # reference length (contains N)
    cov = 1. * sam.mapped / totlen # (mapped may contains secondary or supplementary alignments)
    print (cov)


@command.add_sub
@argument('bamfile')
@argument('-r', '--rsample', type=int, default=1000)
@argument('--pacbio', action='store_true', help='use suplementary alignemnt, use only aligned length')
@argument('--nbed', help='N bed regions')
@argument.exclusive(
    argument('--chroms', nargs='+'),
    argument('--autosomal', action='store_true'),
)
def bam_mean_depth(args):
    """ Roughly estimate bam depth

    - estimate mean read length by scanning top {rsample} reads
    - divide the total length by total reference lengths

    # TODO autosomal option
    """
    sam = pysam.Samfile(args.bamfile, 'rb')
    samm = SamMethods(sam)
    if args.autosomal:
        if 'chr1' in sam.references:
            chroms = ['chr' + str(i + 1) for i in range(22)]
        else:
            chroms = [str(i + 1) for i in range(22)]
    else:
        chroms = args.chroms

    if args.pacbio:
        logging.info('PacBio mode')
        rlen = samm.estimate_rlen(rsample=args.rsample, chroms=chroms, use_alen=True, skip_flag=0x104)
    else:
        rlen = samm.estimate_rlen(rsample=args.rsample, chroms=chroms)
    totlen = samm.get_reference_length(chroms=chroms)

    # negating N bases if available
    nlen = 0
    if args.nbed:
        nbed = Nbed(args.nbed)
        if chroms is None:
            nlen = sum(nbed.get_length(ref) for ref in sam.references)
        else:
            chrom_set = set(chroms)
            nlen = sum(nbed.get_length(ref) for ref in sam.references if ref in chrom_set)

    logging.info('Negating N region from total length: %s', nlen)
    mapped = samm.get_mapped(chroms=chroms)
    d = 1. * mapped * rlen / (totlen - nlen) # (mapped may contains secondary or supplementary alignments)
    logging.info('chroms: %s', chroms)
    logging.info('Sampled reads: %s', args.rsample)
    logging.info('Estimated mean read length: %s', rlen)
    logging.info('Total mapped fragments: %s', mapped)
    logging.info('Total refenrece length: %s', totlen)
    print (d)


@command.add_sub
@argument('bamfile')
@argument('-r', '--region')
@argument('-p', '--plot')
def bam_stats_len(args):
    """ Default is length show distribution
    """
    skip_flags = 0x800 | 0x100  # ignore secondary or suppl. alignment
    sam = pysam.Samfile(args.bamfile, mode='rb', check_sq=False)  # check_sq=False is required for PacBio bam
    if args.region:
        samit = sam.fetch(region=args.region)
    else:
        samit = sam

    lengths = []
    for rec in samit:
        logging.debug(rec)
        if rec.flag & skip_flags:
            continue
        length = rec.infer_query_length(always=False)
        if length is None:
            length = rec.query_length
        lengths.append(length)

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    lens = pd.Series(lengths)
    print (lens.describe())
    if args.plot:
        logging.info('Plot to %s', args.plot)
        plot.init()
        lens.hist(bins=100)
        plt.savefig(args.plot)


@command.add_sub
@argument('bamfile')
@argument('-g', '--group', help='tsv of refname group_name (if missing, refname is emited)')
@argument('-c', '--columns', nargs='+', default=['ngroup', 'groups'], help='additional columns (e.g. qlen, NM, AS, ...)')
def bam_refgroup(args):
    """ Default is length show distribution
    """
    skip_flags = 0x004 # ignore unmapped
    sam = pysam.Samfile(args.bamfile, mode='rb', check_sq=False)  # check_sq=False is required for PacBio bam
    it = sam

    reads = {}  # {readname '/' first_or_not}
    if args.group:
        qname_groups = {qname: grp for qname, grp in iter_tabs(open(args.group), sep='\t')}
        groups = {sam.get_tid(name): qname_groups[name] for name in sam.references}
    else:
        groups = {sam.get_tid(name): name for name in sam.references}
    lengths = {sam.get_tid(name): length for name, length in zip(sam.references, sam.lengths)}

    qname_grps = defaultdict(list)
    qname_vals = defaultdict(lambda : defaultdict(list))

    keys = args.columns or []
    getters1 = {
            'group': lambda read, rec: groups[rec.tid],
            'start': lambda read, rec: read.start,
            'end': lambda read, rec: read.end,
            'tid': lambda read, rec: rec.tid,
            'refname': lambda read, rec: rec.reference_name,
            'lclip': lambda read, rec: read.lclip,
            'rclip': lambda read, rec: read.rclip,
            'qlen': lambda read, rec: read.qlen,
            'NM': lambda read, rec: rec.get_tag('NM'),
            'AS': lambda read, rec: rec.get_tag('AS'),
    }
    keys1 = list(filter(lambda k: k in getters1, keys))
    if 'group' not in keys1:
        keys1.append('group')  # required

    for rec in it:
        if rec.flag & skip_flags:
            continue
        read = Read(rec)
        if read.both_clipped:
            continue
        if read.lclip > 0:
            if not read.start > 0:   # not an end clip
                continue
        if read.rclip > 0:
            if not read.end < lengths[rec.tid]:   # not an end clip
                continue
        #if read.has_lclip and read.has_rclip:  # both clip
        #    continue
        grp = groups[rec.tid]
        if rec.is_paired:
            name = '{}/{}'.format(rec.qname, '21'[rec.is_read1])
        else:
            name = rec.qname
        for key in keys1:
            qname_vals[name][key].append(getters1[key](read, rec))

    getters2 = {
        'ngroup': lambda vals: len(vals['group']),
        'uniq': lambda vals: ','.join(list(set(vals['group']))),
        'nuniq': lambda vals: len(set(vals['group'])),
    }
    print ('qname', *keys, sep='\t')
    for qname, vals in qname_vals.iteritems():
        val_map = qname_vals[qname]
        vals = [','.join(map(str, val_map[key]))
            if key in getters1
            else getters2[key](val_map)
            for key in keys]
        print (qname, *vals, sep='\t')


