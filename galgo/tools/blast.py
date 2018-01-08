from __future__ import absolute_import
import logging
from argtools import command, argument
from ..blast import BlastTabFile
from ..utils import iter_tabs
import pysam
from galgo import plot
from ..samutil import parse_region
import numpy as np
import io


def _read_genomes(fai_file):
    with io.open(fai_file) as reader:
        return [(rec[0], int(rec[1])) for rec in iter_tabs(reader)]


@command.add_sub
@argument('blast_tab')
@argument('-g', '--genome', required=True, help='genome file (.fai file)')
@argument('-o', '--outfile', default='/dev/stdout')
@argument('--max-evalue', type=float)
def blast2sam(args):
    """
    Assumed format: blastn -task blastn -evalue 0.1 -db {hg19} -query {fastsa} -outfmt '7 std gaps qlen slen qseq sseq'
    blastn -db {db} -num_threads {njobs} -perc_identity 95 -max_target_seqs 30 -best_hit_overhang 0.1 -query {region_fasta} -evalue 0.1 -outfmt '7 std gaps qlen slen qseq sseq'
    """
    header = {
            'HD': {'VN': '1.0'},
    }

    max_evalue = args.max_evalue
    def test(rec):
        if max_evalue is None:
            return True
        elif dict(rec.tags)['ZE'] <= max_evalue:
            return True

    with io.open(args.blast_tab) as reader:
        tab = BlastTabFile(reader)
        genomes = _read_genomes(args.genome)
        SQ = [{'LN': length, 'SN': ref} for ref, length in genomes]
        header['SQ'] = SQ
        reference_ids = dict((ref, tid) for tid, (ref, length) in enumerate(genomes))

        with pysam.AlignmentFile(args.outfile, 'wh', header=header) as out:
            for rec in tab.iter_sam_record(reference_ids):
                if test(rec):
                    out.write(rec)


@command.add_sub
@argument('blast_tabs', nargs='+', help='blastn with -outfmt 7')
@argument('--max-evalue', type=float)
@argument('--max-overhang-ratio', type=float)
@argument('--min-qratio', type=float)
def blast_table(args):
    max_evalue = args.max_evalue
    import sys
    for i, blast_tab in enumerate(args.blast_tabs):
        with io.open(blast_tab) as reader:
            try:
                tab = BlastTabFile(reader)
                for j, (qname, qtab) in enumerate(tab.iter_query_tables()):
                    header = i == 0 and j == 0
                    mode = 'w' if i == 0 and j == 0 else 'wa'
                    if args.max_evalue is not None:
                        qtab = qtab[qtab.evalue <= args.max_evalue]
                    if args.max_overhang_ratio is not None:
                        qtab = qtab[(qtab.left_overhang + qtab.right_overhang) / qtab.qlen <= args.max_overhang_ratio]
                    if args.min_qratio is not None:
                        qtab = qtab[qtab.qratio >= args.min_qratio]
                    qtab.to_csv(sys.stdout, header=header, index=False, mode=mode, sep='\t')
            except Exception as e:
                logging.warning('Error for %s (%s:%s)', blast_tab, type(e), e)


class Measure():
    def __init__(self, a, b):
        self._a = a
        self._b = b
        self._l = b - a

    def __call__(self, value):
        return 1. * (value - self._a) / self._l


# TODO load gtf file for annotation
# TODO add gridline (and ticks for high identity features)
@command.add_sub
@argument('blast_tab', help='blastn with -outfmt 7')
def blast2dotplot(args):
    import sys
    plot.init()
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    output = '{0}.pdf'.format(args.blast_tab)
    match_rate = Measure(.8, 1.)
    with io.open(args.blast_tab) as reader:
        tab = BlastTabFile(reader)
        for j, (qname, qtab) in enumerate(tab.iter_query_tables()):
            qreg, qstart, qend = parse_region(qname)
            sreg, sstart, send = parse_region(qtab.sname[0])  # FIXME for multiple case
            qoffset = qstart or 0
            soffset = sstart or 0
            plt.figure(figsize=(8., 8.))
            ax = plt.gca()
            ax.set_xlim(sstart, send)  # FIXME for start, end were unavailable
            ax.set_ylim(qstart, qend)  # FIXME for start, end were unavailable
            for row in qtab.itertuples():
                linestyle = '-'
                lw = 2.
                color = ('r' if row.is_reverse else 'b')
                alpha = match_rate(1. - 1. * row.edit_distance / row.aln_len)
                alpha = min(max(0, alpha), 1)
                logging.info('Plotting %s (alpha: %s)', ((row.sstart, row.send), (row.qstart, row.qend)), alpha)
                ax.add_line(Line2D((soffset + row.sstart, soffset + row.send), (qoffset + row.qstart, qoffset + row.qend),
                    linestyle=linestyle, color=color, linewidth=lw, alpha=alpha))
            plt.savefig(output)
