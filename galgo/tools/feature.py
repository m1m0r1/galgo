from argtools import command, argument
import logging
import os
import re
import math
from .. import plot
from itertools import groupby
from collections import defaultdict, namedtuple
from ..utils import cached_property, remove_suffix, get_process_pool, bucket, unzip, iter_tabs
from ..interval import IntervalMixin, interval_pileup, get_containers
from .. import samutil
from ..samutil import parse_region
import pysam
from pysam import TabixFile, asGTF


class GenomeBrowser(object):
    def __init__(self):
        self._tracks = []

    def add_gene_track(self, gtf):
        GeneTrack(gene_models)

    def show(self, outfile):
        logging.info('Plot to %s', outfile)
        for track in self._tracks:
            track.show(args.region)


class GeneTrack(object):
    pass


class GTFRecord(namedtuple('GTFRecord', 'seqname,source,feature,start1,end,score,strand,frame,attribute')):
    """

    seqname - name of the chromosome or scaffold; chromosome names can be given with or without the 'chr' prefix. Important note: the seqname must be one used within Ensembl, i.e. a standard chromosome name or an Ensembl identifier such as a scaffold ID, without any additional content such as species or assembly. See the example GFF output below.
    source - name of the program that generated this feature, or the data source (database or project name)
    feature - feature type name, e.g. Gene, Variation, Similarity
    start - Start position of the feature, with sequence numbering starting at 1.
    end - End position of the feature, with sequence numbering starting at 1.
    score - A floating point value.
    strand - defined as + (forward) or - (reverse).
    frame - One of '0', '1' or '2'. '0' indicates that the first base of the feature is the first base of a codon, '1' that the second base is the first base of a codon, and so on..
    attribute - A semicolon-separated list of tag-value pairs, providing additional information about each feature.

    chr22   protein_coding  transcript      42095503        42195458        .       +       .       gene_id "ENSG00000167077"; transcript_id "ENST00000401548"; gene_name "MEI1"; gene
    """
    _pat_attr = re.compile('(\w+) "([^;"]*)";')

    def __new__(cls, seqname, source, feature, start1, end, score, strand, frame, attribute):
        start1 = int(start1)
        end = int(end)
        self = super(GTFRecord, cls).__new__(seqname, source, feature, start1, end, score, strand, frame, attribute)
        for match in self._pat_attr.findall(attribute):
            g = match.groups()
            g[0], g[1]

        return self

    def as_interval(self):
        return Interval(self.seqname, self.start, self.end, data=self)

    @property
    def start(self):
        return self.start1 - 1


class GTFFile(object):
    """
    +----------+----------+-------------------------------+
    |*Column*  |*Name*    |*Content*                      |
    +----------+----------+-------------------------------+
    |1         |contig    |the chromosome name            |
    +----------+----------+-------------------------------+
    |2         |feature   |The feature type               |
    +----------+----------+-------------------------------+
    |3         |source    |The feature source             |
    +----------+----------+-------------------------------+
    |4         |start     |genomic start coordinate       |
    |          |          |(0-based)                      |
    +----------+----------+-------------------------------+
    |5         |end       |genomic end coordinate         |
    |          |          |(0-based)                      |
    +----------+----------+-------------------------------+
    |6         |score     |feature score                  |
    +----------+----------+-------------------------------+
    |7         |strand    |strand                         |
    +----------+----------+-------------------------------+
    |8         |frame     |frame                          |
    +----------+----------+-------------------------------+
    |9         |attributes|the attribute field            |
    +----------+----------+-------------------------------+

    GTF formatted entries also define the following fields that
    are derived from the attributes field:

    +--------------------+------------------------------+
    |*Name*              |*Content*                     |
    +--------------------+------------------------------+
    |gene_id             |the gene identifier           |
    +--------------------+------------------------------+
    |transcript_id       |the transcript identifier     |
    +--------------------+------------------------------+
    """
    def __init__(self, gtf):  # if not seekable, load all info
        self._gtf = gtf
        self._parser = asGTF()
        self.has_tabix = os.path.exists('{0}.tbi'.format(self._gtf))

    @cached_property
    def _data(self):  # load all
        return list(iter(self))

    def __iter__(self):
        with open(self._gtf) as fp:
            for line in fp:
                yield self._parser(line)

    def _fetch(self, region):
        if not self.has_tabix:
            raise Exception('Currently, tabix is required for region query')
        with TabixFile(self._gtf, parser=self._parser) as tabix:
            for row in tabix.fetch(region):
                yield row

    def get_gene_models(self, region=None):
        if region:
            return GeneModels(self._fetch(region), region=region)
        return GeneModels(self._data)


class Transcript(IntervalMixin):
    def __init__(self, gtf_rows):
        self.contig = None
        self.start = None
        self.end = None
        self.exons = []
        self.CDSs = []
        self.UTRs = []
        self.start_codon = None
        self.stop_codon = None
        self.row = None
        self.source = None
        self.id = None
        self.gene_id = None
        self.gene_name = None
        for row in gtf_rows:
            attrs = row.asDict()
            self.id = attrs['transcript_id']
            self.gene_id = attrs['gene_id']
            if row.feature == 'transcript':
                self.row = row
                self.gene_name = attrs['gene_name']
            elif row.feature == 'exon':
                self.exons.append(row)
            elif row.feature == 'CDS':
                self.CDSs.append(row)
            elif row.feature == 'UTR':
                self.UTRs.append(row)
            elif row.feature == 'start_codon':
                self.start_codon = row
            elif row.feature == 'stop_codon':
                self.stop_codon = row

        if self.row:
            self.contig = self.row.contig
            self.start = self.row.start
            self.end = self.row.end
            self.strand = self.row.strand
            self.source = self.row.source

        # sort
        if self.exons:
            self.exons = sorted(self.exons, key=lambda x:int(x.asDict()['exon_number']))
        if self.CDSs:
            self.CDSs = sorted(self.CDSs, key=lambda x:int(x.asDict()['exon_number']))

    def __str__(self):
        return 'Transcript[{self.gene_id}({self.gene_name}):{self.id},{self.source},{self.contig}:{self.start}-{self.end}({length}),nexons:{nexon},CDS:{nCDS}]'.format(
                self=self, nexon=len(self.exons), nCDS=len(self.CDSs), length=self.end - self.start)


class GeneModels(object):
    def __init__(self, gtf_rows, region=None):
        self.region = region
        self.transcripts = trs = []

        gtf_rows = sorted(gtf_rows, key=lambda x: x.asDict()['gene_id'])
        for gene_id, rows in groupby(gtf_rows, lambda x: x.asDict()['gene_id']):
            rows = list(rows)
            tr_rows = defaultdict(list)
            gene_row = None
            for row in rows:
                if row.feature == 'gene':
                    gene_row = row
                    continue
                tr_rows[row.asDict()['transcript_id']].append(row)
            for rows in tr_rows.values():
                tr = Transcript(rows)
                trs.append(tr)

        # Reorder by start asc, length desc
        trs.sort(key=lambda x: (x.start, - x.end))
        for tr in trs:
            logging.info(tr)

    def plot(self, ax=None, text_size=2.5, gene_only=False):
        import matplotlib
        import matplotlib.pyplot as plt
        import matplotlib.patches as patches
        layouts = []
        trs = self.transcripts  # TODO filter some?
        levels = interval_pileup(trs)

        y_margin = 70
        width_CDS = 40
        width_exon = 25
        strand_styles = {'+': '>', '-': '<', '.': '--'}

        containers = get_containers(trs)
        if len(containers) != 1:
            logging.error('Cannot plot genes of multiple or zero contigs: %s', containers)
        container = tuple(containers.values())[0]
        if self.region is None:
        #matplotlib.rcParams['font.size'] = 2.  # TODO local
            xlim = (container.start, container.end)
        else:
            parsed = parse_region(self.region)
            xlim = (parsed[1], parsed[2])

        def plot_transcript(tr, level):
            y = - level * y_margin
            # basic line (TODO arrows for strand)
            ls = strand_styles[tr.strand]  # TODO
            ax.plot((tr.start, tr.end), (y, y), linestyle='-', color='b', linewidth=0.5)
            # exon
            for row in (tr.exons or []):
                ax.add_patch(patches.Rectangle((row.start, y - width_exon / 2.), row.end - row.start, width_exon, facecolor='b', linewidth=0))
            for row in (tr.CDSs or []):
                ax.add_patch(patches.Rectangle((row.start, y - width_CDS / 2.), row.end - row.start, width_CDS, facecolor='b', linewidth=0))
            #
            tmpl = '{tr.gene_name}' if gene_only else '{tr.id}({tr.gene_name})'
            ax.text((tr.start + tr.end)/2., - (level-0.5) * y_margin, tmpl.format(tr=tr), ha='center', size=text_size)

        if ax is None:
            ax = plt.gca()
        ax.axes.get_yaxis().set_visible(False)
        ax.set_xlim(*xlim)
        ax.set_ylim(- (max(levels) + 1) * y_margin, 1 * y_margin)
        for tr, level in zip(trs, levels):
            plot_transcript(tr, level)
        #--------------------
        #   -------------
        #   ------------
        #   ----------
        #browser.add_gene_track(args.gtf)   # check if seekable
        #browser.set_region(args.region)
        #browser.show(outfile)


@command.add_sub
@argument('gtf', help='need tabix')
@argument('region', help='e.g. chr1:42100000-42150000')
@argument('-o', '--output')
def gtf_plot(args):
    outfile = 'gtf_plot.{region}.pdf'.format(region=args.region)
    plot.init()
    import matplotlib.pyplot as plt
    logging.info('plot to %s', outfile)
    plt.figure(figsize=(10., 2.))
    gtf = GTFFile(args.gtf)
    gm = gtf.get_gene_models(region=args.region)
    gm.plot()
    plt.savefig(outfile)


class SAMAnalyzer(object):
    def __init__(self, fname, region):
        self.filename = fname
        self.region = region
        self.chrom, self.start, self.end = parse_region(region)
        sam = pysam.Samfile(fname)
        self.chrom_length = dict(zip(sam.references, sam.lengths))[self.chrom]
        self.start = self.start or 0
        self.end = self.end or self.chrom_length

        self.aligned = []
        self.mapq1s = []
        self.lclips = []
        self.rclips = []
        self.lmiss = []
        self.rmiss = []
        self.linvs = []
        self.rinvs = []
        for rec in sam.fetch(reference=self.chrom, start=self.start, end=self.end):
            read = samutil.Read(rec)
            if read.unmapped:
                continue
            if read.lclip > 0 and read.rclip > 0:  # discard both clipped reads
                continue
            self.aligned.append(read)
            if read.mapq <= 1:
                self.mapq1s.append(read)
            if read.lclip > 0:
                self.lclips.append(read)
            if read.rclip > 0:
                self.rclips.append(read)
            if read.mate_miss:
                if read.is_left:
                    self.rmiss.append(read)
                else:
                    self.lmiss.append(read)
            if read.mate_invert:
                if read.is_left:
                    self.rinvs.append(read)
                else:
                    self.linvs.append(read)

        #self.lmiss_x = [(a.start + a.end)/2. for a in lmiss]
        #self.rmiss_x = [(a.start + a.end)/2. for a in rmiss]
        #self.linvs_x = [(a.start + a.end)/2. for a in linvs]
        #self.rinvs_x = [(a.start + a.end)/2. for a in rinvs]

    def get_width(self):
        return self.end - self.start

    def plot_amount(self, ax=None, bin_size=1000):
        import matplotlib.pyplot as plt
        ax = ax or plt.gca()
        bins = int(math.ceil((self.end - self.start) / bin_size))
        aligned_x = [(a.start + a.end)/2. for a in self.aligned]
        mapq1s_x = [(a.start + a.end)/2. for a in self.mapq1s]
        ax.hist(aligned_x, bins=bins, range=(self.start, self.end), linewidth=0)
        ax.hist(mapq1s_x, bins=bins, range=(self.start, self.end), linewidth=0, color='silver')

    def plot_clips(self, ax=None, bin_size=1000):
        import matplotlib.pyplot as plt
        ax = ax or plt.gca()
        bins = int(math.ceil((self.end - self.start) / bin_size))
        lclips_x = [a.start for a in self.lclips]
        rclips_x = [a.end for a in self.rclips]
        ax.hist(lclips_x, bins=bins, range=(self.start, self.end), linewidth=0, color='red')
        ax.hist(rclips_x, bins=bins, range=(self.start, self.end), linewidth=0, color='green')


def get_sample_name(fname):
    name = os.path.basename(fname)
    name = remove_suffix(name, '.sam')
    name = remove_suffix(name, '.bam')
    name = remove_suffix(name, '.sorted')
    name = remove_suffix(name, '.sort')
    return name


# TODO density or read depth
@command.add_sub
@argument('sam', nargs='+')
@argument.exclusive(
    argument('-r', '--regions', nargs='+'),
    argument('-R', '--region-file'),
    required=True,
)
@argument('--sample', nargs='+')
@argument('--gtf', help='need tabix')
@argument('--show-clips', action='store_true')
@argument('-o', '--output', default='bam_bp_plot.pdf')
@argument('--sam-per-page', default=10, type=int)
@argument('--bin-size', default=1000, type=int)
@argument('--njobs', default=8, type=int)
@argument('--style', default='darkgrid')
@argument('--max-depth', type=int)
@argument('--xlab-rot', type=int)
@argument('--region-propto', action='store_true', help='set column width to that proportional to each region length')
#@argument('--skip-flag', type=lambda x: int(x, 16), default=0x904)
def bam_bp_plot(args):
    """ Show Breakpoints

    [break_point -> read]   # position, right, left, read properties (has_another_break, aligned portion, ...)
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib import ticker
    import seaborn as sns
    #sns.despine()
    #binned_read_stream = read_stream.subscribe(binning)
    #binned_read_stream.subscribe

    samass = []
    names = []
    region_names = []
    if args.regions:
        regions = args.regions
    else:
        regions = []
        for row in iter_tabs(open(args.region_file)):
            regions.append(row[0])
            if len(row) > 1:
                region_names.append(row[1])
        if not region_names:
            region_names = regions

    def task(a):
        logging.info('Analyzing %s', a)
        i, fname = a
        sample = args.sample[i] if args.sample else get_sample_name(fname)
        samas = [SAMAnalyzer(fname, region=region) for region in regions]
        return samas, sample

    #pool = get_process_pool(args.njobs)
    #for sama, sample in pool.imap(task, enumerate(args.sam), 1):
    for samas, sample in map(task, enumerate(args.sam)):
        samass.append(samas)
        names.append(sample)

    region_widths = [sama.get_width() for sama in samass[0]]

    #for i, fname in enumerate(args.sam):

    for region in regions:
        chrom, start, end = parse_region(region)
        if start:
            logging.info('Plotting %s:%s-%s', chrom, start+1, end)
        else:
            logging.info('Plotting %s', chrom)

    zoom = 1
    matplotlib.rcParams['figure.figsize'] = (11.7 * zoom, 8.3 * zoom)
    matplotlib.rcParams['font.size'] = 8.
    matplotlib.rcParams['axes.titlesize'] = 'medium'
    matplotlib.rcParams['axes.labelsize'] = 'medium'
    matplotlib.rcParams['axes.labelpad'] = 8.
    matplotlib.rcParams['xtick.labelsize'] = 'small'
    matplotlib.rcParams['ytick.labelsize'] = 'small'
    matplotlib.rcParams['legend.fontsize'] = 'small'

    #counters = ReadCountGenerator(sam, chrom, start=start, end=end, window=args.window, skip_flag=args.skip_flag, mass='overlap')
    #tab = pd.DataFrame([getter(counter) for counter in counters], columns=attrs)
    #tab['mid'] = (tab['start'] + tab['end']) / 2.
    #tab['density'] = 1. * tab['covlen'] / window   # like depth
    #tab['density_mapq1'] = 1. * tab['covlen_mapq1'] / window   # like depth

    #plt.tight_layout()
    #plt.title('{0}:{1}-{2}', chrom, start, end)
    # tracks

    def plot_page(samass, names):
        track_n = track_per_sample * len(samass)
        ratios = [1] * (track_per_sample * len(samass))

        if args.gtf:
            track_n += 1
            ratios.append(3)
        region_n = len(regions)

        gridspec_kw = {'height_ratios': ratios}
        gridspec_kw['width_ratios'] = [1] * region_n
        if args.region_propto:
            gridspec_kw['width_ratios'] = region_widths
        fig, sub_ax_list = plt.subplots(track_n, region_n, sharex='col', sharey='row', gridspec_kw=gridspec_kw)
        if track_n == 1:
            sub_ax_list = [sub_ax_list]
        fig.subplots_adjust(hspace=0, wspace=0)

        sub_axss = iter(sub_ax_list)
        bin_size = args.bin_size
        yrot = 0
        for samas, name in zip(samass, names):
            sub_axs = next(sub_axss)
            #plt.setp([a.get_yticklabels() for a in sub_axs[1:]], visible=False)  # remove axis lables
            for sama, ax in zip(samas, sub_axs):
                logging.info('Plotting %s (%s)', sama.filename, sama.region)
                sama.plot_amount(ax=ax, bin_size=bin_size)
            ax = sub_axs[0]
            ax.set_ylabel('#Reads\n{0}'.format(name), rotation=yrot, ha='right')
            if args.max_depth:
                ax.set_ylim((0, args.max_depth))
            if args.show_clips:
                sub_axs = next(sub_axss)
                for sama, ax in zip(samas, sub_axs):
                    sama.plot_clips(ax=ax, bin_size=bin_size)
                ax = sub_axs[0]
                ax.set_ylabel('#Clips\n{0}'.format(name), rotation=yrot, ha='right')

        #for i, (ax, x) in enumerate(zip(sub_axs, bp_xs), 1):
        #    #plt.subplot(n, 1, i, sharex=True)
        #    ax.hist(x, bins=bins, range=(start, end), linewidth=0)
            #ax = plt.gca()
            #ax.xaxis.set_major_formatter(ticker.NullFormatter())
            #ax.xaxis.set_major_locator(ticker.NullLocator())
            #ax.xaxis.set_major_locator(ticker.LinearLocator())

        if args.gtf:
            sub_axs = next(sub_axss)
            for region, ax in zip(regions, sub_axs):
                gtf = GTFFile(args.gtf)
                gm = gtf.get_gene_models(region=region)
                gm.plot(ax=ax, text_size=6.5, gene_only=True)
                ax = plt.gca()
                ax.set_ylabel('Genes', rotation=yrot)

        #ax = plt.gca()
        #ax.xaxis.set_major_formatter(ticker.ScalarFormatter())
        #ax.xaxis.set_major_locator(ticker.LinearLocator())
        sub_axs = sub_ax_list[-1]
        for i, (region, ax) in enumerate(zip(regions, sub_axs)):
            chrom, start, end = parse_region(region)
            if region_names:
                region_name = region_names[i]
            else:
                if start is not None and end is not None:
                    region_name = '{0}:{1}-{2} ({3})'.format(chrom, start+1, end, (end - start))
                else:
                    region_name = chrom
            logging.info('Plotting %s', region_name)
            ax.set_xlabel(region_name, rotation=args.xlab_rot)

    track_per_sample = 1
    if args.show_clips:
        track_per_sample += 1
    from matplotlib.backends.backend_pdf import PdfPages
    pp = PdfPages(args.output)

    for grp in bucket(zip(samass, names), args.sam_per_page):
        samass1, names = unzip(grp, 2)
        with sns.axes_style(args.style, rc={'axes.linewidth': 2.}):
            plot_page(samass1, names)
        pp.savefig(plt.gcf())
    pp.close()
