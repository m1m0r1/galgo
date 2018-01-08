from __future__ import print_function
from argtools import command, argument
import logging
import pysam
from operator import attrgetter
from ..samutil import ReadCountGenerator, sam_intervals
from ..interval import Interval, parse_bed


@command.add_sub
@argument('sam')
@argument.exclusive(
  argument('-r', '--regions', nargs='*', help='region e.g. one based indexes chr1, chr1:1001, or chr1:1001-2000'),
  argument('-b', '--bed', help='region defined bed file'))
@argument('-w', '--window', type=int, default=50, help='width of counter')
@argument('--offset', type=int, default=0, help='offset of regions (or bed file)')
@argument('-H', '--no-header', dest='header', action='store_false', default=True, help='add header')
@argument('--skip-flag', type=lambda x: int(x, 16), default=0x904)
@argument('-m', '--mass', choices=['middle', 'start', 'end', 'overlap'], default='middle', help='mass points of read used to determine the read is inside of the window')
@argument('--no-mate', dest='use_mate', action='store_false', default=True)
def bam_read_counter(args):
    sam = pysam.Samfile(args.sam)

    with sam:
        attrs = ['start', 'end', 'count', 'mapq1', 'clip', 'bclip']
        if args.use_mate:
            mate_info = ['unmapped', 'jumped', 'far', 'overlap', 'penetrate', 'rr', 'ff']
            attrs += ['mate.{0}'.format(m) for m in mate_info]
        attrs += ['covlen', 'covlen_mapq1']

        if args.header:
            print ('chrom', *attrs, sep='\t')

        if args.bed:
            it = parse_bed(args.bed, offset=args.offset)
        else:
            it = sam_intervals(sam, args.regions)
        for (chrom, start, end, _) in it:
            counters = ReadCountGenerator(sam, chrom, start=start, end=end, window=args.window, skip_flag=args.skip_flag, mass=args.mass)

            getter = attrgetter(*attrs)
            try:
                for counter in counters:
                    print (chrom, *getter(counter), sep='\t')
            except Exception as e:
                logging.error(e)

@command.add_sub
@argument('sam')
@argument.exclusive(
  argument('-r', '--regions', nargs='*', help='regions e.g. one based indexes chr1, chr1:1001, or chr1:1001-2000'),
  argument('-b', '--bed', help='region defined bed file'))
@argument('-w', '--window', type=int, default=50, help='width of counter')
@argument('--offset', type=int, default=0, help='offset of regions (or bed file)')
@argument('--skip-flag', type=lambda x: int(x, 16), default=0x904)
@argument('--row', type=int, default=3)
@argument('--col', type=int, default=2)
@argument('--use-name', action='store_true')
@argument('--name-col', type=int, default=4)
@argument('-o', '--output', default='output.pdf')
@argument('--font-size', type=float, default=20.)
@argument('--min-count', type=int, default=0)
def bam_read_count_viz(args):
    sam = pysam.Samfile(args.sam)
    window = args.window
    logging.info('Output : %s', args.output)
    logging.info('Window size : %s', window)

    import pandas as pd
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.cm as cm
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    import seaborn
    def vline(x, *args, **kwds):
        plt.gca().add_line(matplotlib.lines.Line2D([x, x], plt.gca().get_ylim(), *args, **kwds))

    def hline(y, *args, **kwds):
        plt.gca().add_line(matplotlib.lines.Line2D(plt.gca().get_xlim(),[y, y], *args, **kwds))

    def subplot_gen(row, col, on_page_change=None):
        page = 1
        while 1:
            if page > 1 and on_page_change:
                on_page_change()
            plt.figure()
            for i in xrange(row * col):
                yield page, plt.subplot(row, col, i+1)
            page += 1

    def on_page_change():
        pp.savefig(plt.gcf())

    with sam:
        attrs = ['start', 'end', 'count', 'mapq1', 'clip', 'bclip']
        mate_info = ['unmapped', 'jumped', 'far', 'overlap', 'penetrate', 'rr', 'ff']
        attrs += ['mate.{0}'.format(m) for m in mate_info]
        attrs += ['covlen', 'covlen_mapq1']
        getter = attrgetter(*attrs)

        with PdfPages(args.output) as pp:
            zoom = 1
            matplotlib.rcParams['figure.figsize'] = (11.7 * zoom, 8.3 * zoom)
            matplotlib.rcParams['font.size'] = args.font_size
            matplotlib.rcParams['axes.titlesize'] = 'medium'
            matplotlib.rcParams['axes.labelsize'] = 'medium'
            matplotlib.rcParams['axes.labelpad'] = 8.
            matplotlib.rcParams['xtick.labelsize'] = 'small'
            matplotlib.rcParams['ytick.labelsize'] = 'small'
            matplotlib.rcParams['legend.fontsize'] = 'small'

            if args.bed:
                it = parse_bed(args.bed, offset=args.offset)
            else:
                it = sam_intervals(sam, args.regions)
            subplot = subplot_gen(args.row, args.col, on_page_change=on_page_change)
            for (chrom, start, end, data) in it:
                logging.info('Plotting %s:%s-%s', chrom, start+1, end)
                counters = ReadCountGenerator(sam, chrom, start=start, end=end, window=args.window, skip_flag=args.skip_flag, mass='overlap')
                tab = pd.DataFrame([getter(counter) for counter in counters], columns=attrs)
                tab['mid'] = (tab['start'] + tab['end']) / 2.
                tab['density'] = 1. * tab['covlen'] / window   # like depth
                tab['density_mapq1'] = 1. * tab['covlen_mapq1'] / window   # like depth

                density_mean = tab['density'].mean()  # like depth
                logging.debug('\n%s', tab.head())
                logging.debug('\n%s', tab.tail())

                if len(tab) < args.min_count:  # skip plotting
                    logging.info('skipped')
                    continue

                page, ax = next(subplot)
                plt.tight_layout()
                #plt.title('{0}:{1}-{2}', chrom, start, end)
                #ax = plt.gca()
                if len(tab):
                    tab.plot(kind='area', x='mid', y='density', ax=ax)
                    tab.plot(kind='area', x='mid', y='density_mapq1', ax=ax)
                if start is not None and end is not None:
                    ax.set_xlabel('{0}:{1}-{2} ({3})'.format(chrom, start+1, end, (end - start)))
                else:
                    ax.set_xlabel('{0}'.format(chrom))
                hline(density_mean, color='red')
                if args.use_name:
                    title = data[args.name_col - 1]
                    plt.title(title)
            on_page_change()
