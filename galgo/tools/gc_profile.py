#!/usr/bin/env python
from __future__ import print_function
from argtools import command, argument
from itertools import groupby, chain
from collections import namedtuple, defaultdict
import logging
from operator import attrgetter, itemgetter
import re
import sys
import os
import gzip
import numpy as np
from ..utils import iter_tabs, with_header, safeint, safefloat, safediv, weighted_choice, Counter
from ..interval import Interval, iter_bin_interval_flanks, align_intervals, ChromCompare
from pysam import Fastafile, Tabixfile



def validated_bin_intervals(intervals, bin_size):
    """ Confirm chrom ordered, sorted in start position, and equally binned.
    """
    chrom_set_done = set()
    chrom = None
    start = -bin_size

    for chrom, ivs in groupby(intervals, itemgetter(0)):  # group by chrom
        assert chrom not in chrom_set_done, '{0} was appeared again!'.format(chrom)
        chrom_set_done.add(chrom)
        start = -bin_size
        for iv in ivs:
            assert start < iv.start, 'Not sorted in start position!'
            assert iv.end - iv.start == bin_size, 'Interval {0} != bin size: {1}'.format(iv.end - iv.start, bin_size)
            assert iv.start % bin_size == 0, 'Start position {0} should be divided by bin size: {1}'.format(iv.start, bin_size)
            start = iv.start
            yield iv


@command.add_sub
@argument('nuc', help='output of fa_nuc')
@argument('-s', '--bin-size', type=int, default=50)
@argument('-l', '--left-flank', type=int, default=200, help='size of the left side of the window start (a multiple of the window size)')
@argument('-r', '--right-flank', type=int, default=200, help='size or the right side of the window start (a multiple of the window size)')
@argument('--no-validate', dest='validate', action='store_false', default=True)
def add_gc_flanks(args):
    """
    Required columns:

    """
    logging.info('Validation: %s', args.validate)
    logging.info('Window size: %s bp', args.bin_size)
    logging.info('Left flank: %s bp', args.left_flank)
    logging.info('Right flank: %s bp', args.right_flank)

    assert args.left_flank % args.bin_size == 0, 'Size of the left flank should be a multiple of the bin size'
    assert args.right_flank % args.bin_size == 0, 'Size of the right flank should be a multiple of the bin size'

    lbins = args.left_flank / args.bin_size
    rbins = args.right_flank / args.bin_size

    lines = open(args.nuc)
    nucs = with_header(iter_tabs(lines), use_header=True)
    header = next(nucs)

    def iter_ivs(nucs):
        for chrom, nucs in groupby(nucs, lambda x: x['chrom']):
            for nuc_rec in nucs:
                #cov_rec['GC'] = int(rec['G']) + int(rec['C'])
                iv = Interval(chrom, int(nuc_rec['start']), int(nuc_rec['end']), data=nuc_rec)
                yield iv

    intervals = iter_ivs(nucs)

    if args.validate:
        intervals = validated_bin_intervals(intervals, args.bin_size)

    iv_flanks = iter_bin_interval_flanks(intervals, bin_size=50, lflank=args.left_flank, rflank=args.right_flank)

    # summing up
    bin_size = args.bin_size
    flank_size = args.left_flank + args.right_flank
    header = header + ('flank_gc', 'flank_nonn', 'flank_gc_ratio')
    print (*header, sep='\t')
    for chrom, recs in groupby(iv_flanks, lambda x: x.iv.contig):
        G = C = N = 0
        for rec in recs:
            # calculation by difference

            for iv in rec.adds:
                G += int(iv.data['Gs'])
                C += int(iv.data['Cs'])
                N += int(iv.data['Ns'])
            for iv in rec.removes:
                G -= int(iv.data['Gs'])
                C -= int(iv.data['Cs'])
                N -= int(iv.data['Ns'])

            iv = rec.iv
            flank_gc =  G + C
            flank_nonn = len(rec.flanks) * bin_size - N
            iv.data['flank_gc'] = flank_gc
            iv.data['flank_nonn'] = flank_nonn
            iv.data['flank_gc_ratio'] = 1. * flank_gc / flank_nonn if flank_nonn else 0.
            print (*iv.data.values(), sep='\t')


@command.add_sub
@argument('fasta')
@argument('-s', '--bin-size', type=int, default=50)
@argument('-N', '--omit-N', action='store_true', default=False)
@argument('-H', '--no-header', dest='header', action='store_false', default=True)
@argument('--chroms', nargs='*', default=None)
def fa_nuc(args):
    """
    """
    header = 'chrom start end As Cs Gs Ts Ns'.split(' ')
    bin_size = args.bin_size
    if args.chroms:
        chroms = list(line.rstrip('\r\n').split('\t')[0] for line in open(args.chrom_list))
    else:
        chroms = None
    if args.header:
        print (*header, sep='\t')
    for chrom, start, end, seq in iter_fa_seq(args.fasta, bin_size=bin_size, chroms=chroms):
        counter = Counter(seq)
        if len(seq) < bin_size:
            counter['N'] += bin_size - len(seq)

        if args.omit_N and counter['N'] == bin_size:
            continue

        print (chrom, start, end, counter['A'], counter['C'], counter['G'], counter['T'], counter['N'], sep='\t')


def iter_fa_seq(fasta, bin_size=50, chroms=None):
    fasta = Fastafile(fasta)
    if chroms is None:
        chroms = fasta.references
    for chrom in chroms:
        start = 0
        while 1:
            seq = fasta.fetch(chrom, start, start + bin_size).upper()
            if not seq:
                break
            yield chrom, start, start + bin_size, seq
            start += bin_size


@command.add_sub
@argument('bed1')
@argument('bed2')
@argument('-H', '--no-header', dest='use_header', action='store_false', default=True)
@argument('-l', '--chrom-list', required=True, help='list of chromosome to use in order')
def bed_join(args):
    if args.use_header:
        def get_bed_ivs(bed):
            bed = iter_tabs(bed)
            bed = with_header(bed, use_header=True)
            header = next(bed)
            ivs = (Interval(rec['chrom'], int(rec['start']), int(rec['end']), data=rec) for rec in bed
                    if rec['chrom'] in chrom_set)
            return header, ivs

        def get_row(iv):
            iv1 = iv.data[0]
            iv2 = iv.data[1]
            return (iv.contig, iv.start, iv.end,) + tuple(iv1[k] for k in data_key1) + tuple(iv2[k] for k in data_key2)

    else:
        def get_bed_ivs(bed):
            bed = iter_tabs(bed)
            return None, (Interval(rec[0], int(rec[1]), int(rec[2]), data=rec) for rec in bed
                            if rec[0] in chrom_set)

        def get_row(iv):
            iv1 = iv.data[0]
            iv2 = iv.data[1]
            return (iv.contig, iv.start, iv.end,) + tuple(iv1[3:]) + tuple(iv2[3:])

    header1, ivs1 = get_bed_ivs(open(args.bed1))
    header2, ivs2 = get_bed_ivs(open(args.bed2))

    if args.use_header:
        data_key1 = header1[3:]
        data_key2 = header2[3:]
        # header
        print (*(('chrom', 'start', 'end') + data_key1 + data_key2), sep='\t')

    chroms = list(line.rstrip('\r\n').split('\t')[0] for line in open(args.chrom_list))
    chrom_set = set(chroms)
    comp = ChromCompare(chroms)

    for iv in align_intervals(ivs1, ivs2, comp):
        iv1 = iv.data[0]
        iv2 = iv.data[1]
        if iv1 is None or iv2 is None:
            logging.info('%s was skipped', iv)
            continue
        row = get_row(iv)
        print (*row, sep='\t')


_re_autosomal = re.compile('^(chr[0-9]+|[0-9]+)$')

def is_autosomal(chrom):
    """
    >>> is_autosomal('chr1')
    True
    >>> is_autosomal('19')
    True
    >>> is_autosomal('22')
    True
    >>> is_autosomal('X')
    False
    >>> is_autosomal('chrY')
    False
    >>> is_autosomal('chr6_random')
    False
    """
    return bool(_re_autosomal.match(chrom))


@command.add_sub
@argument('gc_cov')
@argument('--gc-col', default='flank_gc_ratio')
@argument('--low-mapq-col', default='mapq1')
@argument('--gc-bin-size', type=float, default=.02)
@argument('-b', '--low-mapq-bound', type=float, default=.04)
@argument('--chroms', nargs='*', default=None)
def gc_cov_count(args):
    """
    Default chroms are autosomals
    {
        gc_bin_size,
        max_multimaps,
        max_multimaps,
        counts: {gc_bin: {cov: count}},
    }
    """
    logging.info('Chromosomes: %s', args.chroms)
    if args.chroms:
        chrom_set = set(args.chroms)
        def is_target_chrom(chrom):
            return chrom in chrom_set
    else:
        is_target_chrom = is_autosomal

    gc_col = args.gc_col
    low_mapq_col = args.low_mapq_col
    gc_bin_size = args.gc_bin_size
    bound = args.low_mapq_bound

    with open(args.gc_cov) as fp:
        gc_covs = iter_tabs(fp)
        gc_covs = with_header(gc_covs)

        gc_cov_counters = defaultdict(Counter)  # {gc_bin: {cov: count}}
        for chrom, gc_covs in groupby(gc_covs, lambda x: x['chrom']):
            if not is_target_chrom(chrom):
                continue

            logging.info('Processing chrom: %s', chrom)
            for gc_cov in gc_covs:
                cov = int(gc_cov['count'])
                low_cov = int(gc_cov[low_mapq_col])
                gc_ratio = float(gc_cov[gc_col])
                if cov > 0 and 1. * low_cov / cov >= args.low_mapq_bound:
                    continue  # skip low mapq record

                gc_bin = int(gc_ratio / gc_bin_size)
                gc_cov_counters[gc_bin][cov] += 1


    gc_bin_max = int(round(1. / gc_bin_size))
    print ('gc_bin', 'gc_ratio', 'cov', 'count', sep='\t')
    for gc_bin in xrange(gc_bin_max):
        for cov, count in sorted(gc_cov_counters[gc_bin].items()):
            print (gc_bin, '{0:.2f}'.format(gc_bin * gc_bin_size), cov, count, sep='\t')


_gc_stats_header = ['gc_bin', 'gc_ratio', 'nbins_total', 'cov_mean_total', 'nbins_total_trim', 'cov_mean_total_trim', 'nbins', 'cov_mean', 'cov_std', 'nbins_trim', 'cov_mean_trim']
def make_gc_stats(gc_cov_count, out=sys.stdout, min_nbins=100):
    """  # TODO use panads
    gc_bin gc_ratio gc_ratio
    """
    gc_cov_counts = defaultdict(Counter)
    cov_counts = Counter()
    gc_ratios = {}
    q1 = .1
    q2 = .9

    with open(gc_cov_count) as fp:
        for rec in with_header(iter_tabs(fp)):
            gc_bin = int(rec['gc_bin'])
            gc_ratio = float(rec['gc_ratio'])
            cov = int(rec['cov'])
            count = int(rec['count'])
            cov_counts[cov] += count
            gc_cov_counts[gc_bin][cov] += count
            gc_ratios[gc_bin] = gc_ratio

    print (*_gc_stats_header, file=out, sep='\t')

    covs_total = np.array(sorted(cov_counts.keys()))
    weights_total = np.array([cov_counts[cov] for cov in covs_total]).astype('float')
    cov_mean_total = (covs_total * weights_total).sum() / weights_total.sum()

    nbins_total = int(weights_total.sum())
    weights_total_trim = trim_weights(weights_total, q1, q2)
    nbins_total_trim = weights_total_trim.sum().astype('int')
    cov_mean_total_trim = trimmed_mean(covs_total, weights_total, q1, q2)

    logging.info('nbins_total: %s', nbins_total)
    logging.info('cov_mean_total: %s', cov_mean_total)
    logging.info('nbins_total_trim: %s', nbins_total_trim)
    logging.info('cov_mean_total_trim: %s', cov_mean_total_trim)

    for gc_bin in sorted(gc_cov_counts):
        gc_ratio = gc_ratios[gc_bin]
        covs = np.array(list(gc_cov_counts[gc_bin].keys())) * 1.
        bins = list(gc_cov_counts[gc_bin].values())
        weights = np.array(bins) * 1.
        cov_mean = np.average(covs, weights=weights)
        cov_std = np.sqrt((np.average(covs**2, weights=weights) - cov_mean**2).sum())
        #remains = (covs >= cov_mean_total * trim_rate_min) & (covs <= cov_mean_total * trim_rate_max)
        weights_trim = trim_weights(weights, q1, q2)
        nbins_trim = int(weights_trim.sum().round())
        cov_mean_trim = trimmed_mean(covs, weights, q1, q2) if nbins_trim >= min_nbins else float('nan')

        nbins = sum(bins)
        row = [gc_bin, '{0:.2f}'.format(gc_ratio),
            nbins_total, '{0:.4f}'.format(cov_mean_total),
            nbins_total_trim, '{0:.4f}'.format(cov_mean_total_trim),
            nbins, '{0:.4f}'.format(cov_mean), '{0:.4f}'.format(cov_std),
            nbins_trim, '{0:.4f}'.format(cov_mean_trim)]
        print (*row, file=out, sep='\t')


def load_gc_stats(stats_file):
    gc_stats = {}
    with open(stats_file) as fp:
        for row in with_header(iter_tabs(fp)):
            gc_bin = int(row['gc_bin'])
            gc_stats[gc_bin] = row
            for k in 'gc_bin nbins_total nbins_total_trim nbins nbins_trim'.split(' '):
                gc_stats[gc_bin][k] = safeint(row[k], default=row[k])
            for k in 'gc_ratio cov_mean_total cov_mean_total_trim cov_mean cov_std cov_mean_trim'.split(' '):
                gc_stats[gc_bin][k] = safefloat(row[k], default=row[k])

    return gc_stats


def trim_weights(weight, q1=0.1, q2=0.9):
    """ Trim weight outside of quantiles
    Args:
        weight: numpy array or pandas.Series object

    # print (hist_quantile(df['x'], df['w'], 0.5))
    >>> weight = np.array([5, 10, 20, 30, 20, 10, 5])
    >>> list(trim_weights(weight, 0.025, 0.975))
    [2.5, 10.0, 20.0, 30.0, 20.0, 10.0, 2.5]
    >>> list(trim_weights(weight, 0.05, 0.95))
    [0.0, 10.0, 20.0, 30.0, 20.0, 10.0, 0.0]
    >>> list(trim_weights(weight, 0.1, 0.9))
    [0.0, 5.0, 20.0, 30.0, 20.0, 5.0, 0.0]
    >>> list(trim_weights(weight, 0.2, 0.8))
    [0.0, 0.0, 15.0, 30.0, 15.0, 0.0, 0.0]
    >>> list(trim_weights(weight, 0.25, 0.75))
    [0.0, 0.0, 10.0, 30.0, 10.0, 0.0, 0.0]
    >>> list(trim_weights(weight, 0.4, 0.6))
    [0.0, 0.0, 0.0, 20.0, 0.0, 0.0, 0.0]
    >>> list(trim_weights(weight, 0.025, 0.03))
    [0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    >>> list(trim_weights(weight, 0.97, 0.975))
    [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5]
    """
    assert q1 < q2
    assert np.all(weight > 0)
    wcum = np.array([0] + list(weight.cumsum()))
    tot = weight.sum()

    idx1 = wcum[wcum <= tot*q1].shape[0] - 1
    wsub1 = tot * q1- wcum[idx1]
    idx2 = wcum[wcum <= tot*q2].shape[0] - 1
    wsub2 = wcum[idx2] - tot * q2

    new_weight = weight.copy().astype('float')

    new_weight[:idx1] = 0
    new_weight[idx2:] = 0
    new_weight[idx1] -= wsub1
    new_weight[idx2] -= wsub2
    return new_weight

def trimmed_mean(vals, weights, q1, q2):
    w = trim_weights(weights, q1, q2)
    s = w.sum()
    return (1. * (vals * w).sum() / s) if s else float('nan')


@command.add_sub
@argument('gc_cov_count')
@argument('-o', '--output')
@argument('--gc-bin-size', type=float, default=None)
@argument('--plot-box', action='store_true')
def plot_gc_cov_count(args):
    #logging.info('Creating stats file: %s', stats_file)
    #with open(stats_file, 'w+') as out:
    #    make_gc_stats(args.gc_cov_count, out=out)
    #stats_file = '{0}.stats'.format(args.gc_cov_count)
    #stats = load_gc_stats(stats_file)
    import pandas as pd
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.cm as cm
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    try:
        import seaborn
    except ImportError:
        pass

    outname = args.output if args.output else '{0}.pdf'.format(args.gc_cov_count)
    logging.info('Make plot to %s', outname)

    def get_quantile_coverage(tab, q):
        tab = pd.DataFrame({'cov': list(tab['cov']), 'count': list(tab['count'])})   # workaround
        cov_idxs = tab['cov'].argsort()

        cum = tab['count'][cov_idxs].cumsum()
        tot = tab['count'].sum()

        last_idx = cum[1. * cum / tot <= q].last_valid_index()
        return tab['cov'][last_idx] if last_idx is not None else float('nan')

    # setup data
    tab = pd.read_table(args.gc_cov_count)

    q_cov = get_quantile_coverage(tab, .9)
    logging.info('The %s quantile of the coverage: %s', .9, q_cov)
    max_cov = safeint(q_cov * 2., default=10)

    def gen_data_or_sample(values, weights, size=100):
        if len(values) > size:
            return values[weighted_choice(weights, size)]
        else:
            return np.repeat(values, weights)
    # plotting
    pp = PdfPages(outname)

    gc_ratios = []
    cov_1s = []
    cov_3s = []
    cov_25s = []
    cov_75s = []
    cov_samples = []
    cov_trims = []
    for gc_bin, tab1 in tab.groupby('gc_bin'):
        gc_ratio = tab1['gc_ratio'][tab1.index[0]]
        gc_ratios.append(gc_ratio)
        cov_sample = gen_data_or_sample(tab1['cov'].values, tab1['count'].values)
        covs = np.sort(tab1['cov'].values)
        counts = tab1['count'].values[tab1['cov'].values.argsort()]
        #logging.info('%s', (covs, counts))
        cov_samples.append(cov_sample)
        cov_trims.append(trimmed_mean(covs, counts, .1, .9))
        #cov_25s.append(trimmed_mean(covs, counts, .1, .4))
        #cov_75s.append(trimmed_mean(covs, counts, .6, .9))
        cov_1s.append(get_quantile_coverage(tab1, .25))
        cov_3s.append(get_quantile_coverage(tab1, .75))

    #covs = pd.DataFrame({'gc_ratio': gc_ratios, 'cov_1': cov_1s, 'cov_med': cov_meds, 'cov_3': cov_3s})
    fig, ax = plt.subplots()
    ax.set_xlim((0, 1.))
    ax.set_ylim(0, max_cov)
    plt.plot(gc_ratios, cov_trims)
    plt.plot(gc_ratios, cov_1s, ls='dashed')
    plt.plot(gc_ratios, cov_3s, ls='dashed')
    pp.savefig(fig)

    if args.plot_box:
        fig, ax = plt.subplots()
        ax.set_xlim((0, 1.))
        ax.set_ylim(0, max_cov)
        ax.set_xticklabels(ax.xaxis.get_majorticklabels(), rotation=90, size=7)
        plt.boxplot(cov_samples, labels=gc_ratios)

        #covs.plot(kind='line', x='gc_ratio', ax=ax)
        pp.savefig(fig)

    #fig = plt.figure()
    fig, ax = plt.subplots()
    ax.set_xlim((0, max_cov))
    ax.set_ylim((0, 1.))
    x = np.arange(0, 1, 0.025)
    for gc_bin, tab1 in tab.groupby('gc_bin'):
        gc_ratio = tab1['gc_ratio'][tab1.index[0]]
        color = cm.rainbow(gc_ratio, 0.6)
        #logging.info('Color for gc %s : %s', gc_ratio, color)
        s = tab1['count'].sum()
        ax.plot(tab1['cov'].values, 1. * tab1['count'].values / s, c=color, label=gc_ratio)
    ax.legend(tab['gc_ratio'].unique(), loc='upper right', ncol=3, prop={'size': 9}) #bbox_to_anchor=(1.05, 1), ncol=3)
    pp.savefig(fig)

    #fig, ax = plt.subplots()
    #ax.set_xlim((0, 1.))
    #ax.set_ylabel('region count')
    ##ax.set_ylim((0, max_cov))
    ##tab.plot(kind='hexbin', x='gc_ratio', y='cov', C='count', ax=ax)
    #tab.groupby('gc_ratio')['count'].sum().plot(kind='area', ax=ax)
    #tab.plot(kind='scatter', x='gc_ratio', y='cov', ax=ax)
    #pp.savefig(fig)

    pp.close()


@command.add_sub
@argument('gc_cov_count')
@argument('--min-nbins', type=int, default=100)
def gc_cov_stats(args):
    make_gc_stats(args.gc_cov_count, min_nbins=args.min_nbins)


# deprecated
@command.add_sub
@argument('gc_cov')
@argument('-c', '--gc-cov-count', required=True)
@argument('-b', '--low-mapq-bound', type=float, default=.04)
@argument('--gc-bin-size', type=float, default=None)
@argument('--gc-col', default='flank_gc_ratio')
@argument('--min-nbins', type=int, default=100)
def gc_cov_normalize(args):
    stats_file = '{0}.stats'.format(args.gc_cov_count)
    #if not os.path.exists(stats_file) or os.path.getsize(stats_file) == 0 or os.path.getatime(stats_file) < os.path.getatime(args.gc_cov_count):
    logging.info('Creating stats file: %s', stats_file)
    if not os.path.exists(stats_file) or os.path.getatime(stats_file) < os.path.getatime(args.gc_cov_count):
        with open(stats_file, 'w+') as out:
            make_gc_stats(args.gc_cov_count, min_nbins=args.min_nbins, out=out)
    stats = load_gc_stats(stats_file)
    gc_means = dict((gc_bin, row['cov_mean_trim']) for gc_bin, row in stats.items())
    total_mean = list(stats.values())[0]['cov_mean_total_trim']
    gc_col = args.gc_col

    # estimate gc_bin_size
    if args.gc_bin_size is None:
        gc_ratio_lis = np.array(list(sorted(v['gc_ratio'] for v in stats.values())))
        gc_bin_size = np.min(np.diff(gc_ratio_lis))  # calculate from the difference of adjacent bin
    else:
        gc_bin_size = args.gc_bin_size

    logging.info('gc_bin_size: %s', gc_bin_size)

    #for gc_bin in sorted(stats):
    #    print (stats[gc_bin])

    with open(args.gc_cov) as fp:
        it = with_header(iter_tabs(fp), use_header=True)
        keys = next(it)
        header = keys + ('gc_bin', 'cov_normal', 'cov_gc_normal')
        print (*header, sep='\t')

        for row in it:
            vals = [row[k] for k in keys]
            gc_ratio = float(row[gc_col])
            gc_bin = int(round(gc_ratio / gc_bin_size))
            gc_mean = gc_means.get(gc_bin, float('nan'))
            cov_normal = safediv(float(row['count']), total_mean)
            cov_gc_normal = safediv(float(row['count']), gc_mean)

            vals.append(gc_bin)
            vals.append('{0:.4f}'.format(cov_normal))
            vals.append('{0:.4f}'.format(cov_gc_normal))
            print (*vals, sep='\t')


@command.add_sub
@argument('region', default='')
@argument('cov', default='coverage file tabixed')
@argument('-c', '--cov-col', default='cov_gc_normal', help='column name of coverage info')
@argument('-b', '--low-mapq-bound', type=float, default=.04)
@argument('--low-mapq-col', default='mapq1')
@argument('-r', '--use-region-header', action='store_true', default=False)
@argument('-s', '--stats-file', default=None, help='used for normalization')
@argument('--gc-col', default='flank_gc_ratio', help='used for normalizatoin')
@argument('--max-coeff', type=float, default=3.0)
@argument('--min-coeff', type=float, default=0.3)
def cov_intersect(args):
    """

    Output bed file:
      - 1 - region.ncol : region file record
      - extra fields:
        mean coverage
        valid length
        filtered length (due to excess lowq reads or unknown coverage)
        total length
    """
    cov_header = gzip.open(args.cov).readline().rstrip('\r\n').split('\t')
    idx_cov = cov_header.index(args.cov_col)  # column index of coverage data
    idx_count = cov_header.index('count')
    idx_low_mapq = cov_header.index(args.low_mapq_col)
    low_mapq_bound = args.low_mapq_bound
    #header_ext = ['mean_cov', 'region_len', 'valid_len', 'valid_ratio', 'raw_mean_count', 'raw_mean_low_mapq', 'raw_low_mapq_ratio']
    header_ext = ['mean_cov', 'region_len', 'valid_len', 'valid_ratio', 'valid_count', 'valid_exp_count', 'count', 'low_mapq_ratio']
    NEARLY_ZERO = 1e-10

    if args.stats_file:
        # for normalization
        idx_gc = cov_header.index(args.gc_col)
        max_coeff = args.max_coeff
        min_coeff = args.min_coeff

        stats = load_gc_stats(args.stats_file)
        total_mean = list(stats.values())[0]['cov_mean_total_trim']
        gc_covs = {}
        for gc_bin, row in stats.items():
            gc_ratio = row['gc_ratio']
            coeff = 1. * safediv(total_mean, row['cov_mean_trim'], error=float('nan'))   # strength of correction
            gc_covs[gc_bin] = row['cov_mean_trim'] if (min_coeff <= coeff <= max_coeff) else float('nan')
            logging.info('Reference coverage for bin %s, gc_ratio %s: %s (coeff: %s)', gc_bin, gc_ratio, gc_covs[gc_bin], coeff)

        # estimate gc_bin_size
        gc_ratio_lis = np.array(list(sorted(v['gc_ratio'] for v in stats.values())))
        gc_bin_size = np.min(np.diff(gc_ratio_lis))  # calculate from the difference of adjacent bin
        logging.info('gc_bin_size: %s', gc_bin_size)
        def get_cov(rec):
            gc_ratio = float(rec[idx_gc])
            gc_bin = int(round(gc_ratio / gc_bin_size))
            cov = gc_covs.get(gc_bin, float('nan'))
            return safediv(float(rec[idx_count]), cov)

        def get_exp_count(rec):
            gc_ratio = float(rec[idx_gc])
            gc_bin = int(round(gc_ratio / gc_bin_size))
            return gc_covs.get(gc_bin, float('nan'))
    else:
        def get_cov(rec):
            return float(rec[idx_cov])

        def get_exp_count(rec):
            return float('nan')

    with Tabixfile(args.cov) as cov_fp:
        with open(args.region) as reg_fp:
            regs = iter_tabs(reg_fp)

            first = next(regs, None)  # peaking first line of the region file
            if first is None:
                return
            if args.use_region_header:
                header = list(first) + header_ext
            else:
                header = ['chrom', 'start', 'end'] \
                    + ['regcol{0}'.format(i+1) for i in xrange(len(first) - 3)] + header_ext
                regs = chain([first], regs)

            print (*header, sep='\t')
            for reg_row in regs:
                chrom = reg_row[0]
                start = int(reg_row[1])
                end = int(reg_row[2])
                region_len = end - start

                # TODO unit test
                recs = list(iter_tabs(cov_fp.fetch(chrom, start, end)))
                if recs:
                    bin1_start = int(recs[0][1])  # first bin
                    bin1_end  = int(recs[0][2])
                    bin2_start = int(recs[-1][1]) # last bin
                    bin2_end  = int(recs[-1][2])
                    bin_size = bin1_end - bin1_start

                    mapqs = [float(rec[idx_low_mapq]) for rec in recs]  # mapqs of each recs
                    bin_len1 = bin_len2 = bin_size
                    if bin1_start < start:
                        if end < bin1_end:   # the bin contains the whole region
                            bin_len1 = end - start
                        else:
                            bin_len1 = bin1_end - start
                    if end < bin2_end:
                        if bin2_start < start:   # the bin contains the whole region
                            bin_len2 = end - start
                        else:
                            bin_len2 = end - bin2_start

                    exp_counts = np.array([get_exp_count(rec) for rec in recs])
                    covs = np.array([get_cov(rec) for rec in recs])
                    counts = np.array([int(cov[idx_count]) for cov in recs])
                    low_mapqs = np.array([int(cov[idx_low_mapq]) for cov in recs])
                    invalids = (np.isnan(covs) | (1. * (low_mapqs / (counts + NEARLY_ZERO)) >= low_mapq_bound))

                    bin_lens = np.array([bin_size] * len(recs))
                    bin_lens[0] = bin_len1 # weight correction for the left end bin
                    bin_lens[-1] = bin_len2 # weight correction for the right end bin

                    valid_len = bin_lens[~ invalids].sum()
                    mean_cov = (covs * bin_lens)[~ invalids].sum() / valid_len if valid_len else float('nan')
                    valid_ratio = 1. * valid_len / region_len

                    valid_count = 1. * (counts * bin_lens / bin_size)[~ invalids].sum()
                    valid_exp_count = 1. * (exp_counts * bin_lens / bin_size)[~ invalids].sum()
                    # counts[~ invalids].sum() * bin_size / valid_len if valid_len else float('nan')
                    count = 1.* counts.sum() #* bin_size / region_len
                    low_mapq_ratio = safediv(1.* low_mapqs.sum(), count) #* bin_size / region_len
                else:
                    mean_cov = float('nan')
                    valid_len = 0
                    valid_ratio = 0.
                    #valid_mean_count = 0.

                    valid_count = 0.
                    valid_exp_count = float('nan')
                    count = 0.
                    low_mapq_ratio = 0.

                out_row = list(reg_row)
                out_row.append('{0:.3f}'.format(mean_cov))
                out_row.append(region_len)
                out_row.append(valid_len)
                out_row.append('{0:.3f}'.format(valid_ratio))
                #out_row.append('{0:.3f}'.format(valid_mean_count))
                out_row.append('{0:.2f}'.format(valid_count))
                out_row.append('{0:.2f}'.format(valid_exp_count))
                out_row.append('{0:.2f}'.format(count))
                out_row.append('{0:.2f}'.format(low_mapq_ratio))
                #out_row.append('{0:.3f}'.format(1. * raw_mean_low_mapq / raw_mean_count if raw_mean_count else 0))
                print (*out_row, sep='\t')



