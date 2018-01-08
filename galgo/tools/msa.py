from __future__ import absolute_import
from __future__ import print_function
from ..utils import Counter, cached_property
from ..io import open
from ..bioseq import color_term_dna, color_term_aa, mask_same_bases, Fasta, get_aln_pos_text, mask_unknown_seq, SeqLocator, contig2chain, AlignmentAnalyzer, pos2text, get_default_mask_char, STOP_AA
from builtins import zip
from argtools import command, argument
import itertools
import logging
import pandas as pd


def msa_consensus_longest_common(seqs):
    """
    Rules:
    - select most common base for each position other than deletion character ('-')
    - if multiple top bases exist, one of those is selected using Counter(bases).most_common(1)
    - if most common

    >>> seq1 = '----AAA--AAAATTT----'
    >>> seq2 = '-AA-AAA--AACACTT-GGG'
    >>> seq3 = 'G---AAATTAACATTTAGGT'

                GAA AAATTAACATTTAGGT  # the last T is undetermined by definition

    >>> ''.join(msa_consensus_longest_common([seq1, seq2, seq3]))[:-1]
    'GAAAAATTAACATTTAGG'

    >>> ''.join(msa_consensus_longest_common([seq1, seq2, seq3]))[-1] in 'GT'  # differed for python2 or python3
    True
    """
    for bases in zip(*seqs):
        bases = filter(lambda x: x !='-', bases)  # no bases other than padding '-'
        counts = Counter(bases)
        if not counts:  # no bases other than padding '-'
            yield ''
        else:
            cons_base = counts.most_common(1)[0][0]
            yield cons_base


def msa_consensus_filled(seqs, del_as_var=False):
    """
    >>> seq1 = '----AAA--AAAATTT--N-'
    >>> seq2 = '-AA-AAA--AACACTT-GGG'
    >>> seq3 = 'G---AAATTAACATTTAGGT'

                GAANAAATTAANANTTAGGN   # consensus is 

    >>> ''.join(msa_consensus_filled([seq1, seq2, seq3]))
    'GAANAAATTAANANTTAGGN'

    >>> ''.join(msa_consensus_filled([seq1, seq2, seq3], del_as_var=True))
    'NNNNAAANNAANANTTNNGN'
    """
    filter_chars = 'N'
    if not del_as_var:
        filter_chars += '-'

    for bases in zip(*seqs):
        bases = filter(lambda x: x not in filter_chars, bases)  # no bases other than padding '-'
        counts = Counter(bases)
        if not counts:
            yield 'N'
        elif len(counts) == 1:
            base = tuple(counts.keys())[0]
            if base == '-':
                yield 'N'
            else:
                yield base
        else:
            yield 'N'


@command.add_sub
@argument('fasta', help='all contigs should be the same lengths')
@argument('-n', '--name', default='new_contig', help='name of new contig')
@argument('-m', '--mode', default='longest_common', choices=['longest_common', 'filled'])
@argument('--mask-clips', action='store_true')
@argument('--del-as-var', action='store_true')
def msa_consensus(args):
    """ Make a consensus sequence from multiple sequence alignment

    mode:
        1) '-' is a variant
        2) '-' consumes a space
        3) 'N' is a variant
        4) Consensus base
                            1             2     3       4
        - longest_commmon : No          | No  | (Yes) | most common base
        - filled          : del_as_var  | Yes | No    | N if variant or del
    """
    with open(args.fasta) as fp:
        fasta = Fasta(fp)
        if args.mask_clips:
            seqs = list(contig.get_clip_mask_seq(mask_char='N') for contig in fasta)
        else:
            seqs = list(contig.seq for contig in fasta)

    if args.mode == 'longest_common':
        if args.del_as_var:
            logging.warning('del_as_var is ignored for `longest_common` mode')
        new_seq = ''.join(msa_consensus_longest_common(seqs))
    elif args.mode == 'filled':
        new_seq = ''.join(msa_consensus_filled(seqs, del_as_var=args.del_as_var))
    else:
        raise NotImplementedError

    print ('>' + args.name)
    print (new_seq)


@command.add_sub
@argument('fasta', help='all contigs should be the same lengths')
@argument('-f', '--from-name', default='1', help='name of new contig')
@argument('-t', '--to-name', default='2', help='name of new contig')
def msa2chain(args):
    """ Make chain file from a pair of aligned contigs (default pair is 1st and 2nd)
    """
    with open(args.fasta) as fp:
        fasta = Fasta(fp)
        if args.from_name.isdigit():
            from_contig = fasta.contigs[int(args.from_name) - 1]
        else:
            from_contig = fasta.get(args.from_name)
        if args.to_name.isdigit():
            to_contig = fasta.contigs[int(args.to_name) - 1]
        else:
            to_contig = fasta.get(args.to_name)

        chain = contig2chain(from_contig, to_contig)
        print ('chain', *chain['header'], sep=' ')
        for row in chain['data']:
            print (*row, sep='\t')


# TODO realign option
@command.add_sub
@argument('fasta', help='all contigs should be the same lengths')
@argument('-n', '--names', help='comma separated list of names to calculate')
@argument('-r', '--redundant', action='store_true', help='calculate for all the combinations')
@argument('--ctx-size', type=int, default=3, help='when set this, emit CTX base context sequences for edits')
def msa2edit(args):
    """ Emit pairwise sequence difference from multiple aligned contigs

    Emission records:
        name1 :
        name2 :
        pos  : 123  # pos in alignment
        pos1 : 121  # last end in contig 1 (last pos if seq is deletion)
        pos2 : 122  # last end in contig 2 (last pos if seq is deletion)
        seq1 : A
        seq2 : -
        ctx_size : 3
        ctx1 : AAAAAAA
        ctx2 : AAA-AAA

    # TODO collapsing neighbor variants
    """
    assert args.ctx_size > 0
    with open(args.fasta) as fp:
        fasta = Fasta(fp)
        if args.names:
            names = args.names.split(',')
        else:
            names = fasta.names

        if args.redundant:
            pairs = itertools.product(names, names)
        else:
            pairs = itertools.combinations(names, 2)
        ctx_size = args.ctx_size

        print ('name1', 'name2', 'msa_pos', 'pos1', 'pos2', 'base1', 'base2', 'ctx_size', 'ctx1', 'ctx2', sep='\t')
        for n1, n2 in pairs:
            c1 = fasta.get(n1)
            aln1 = c1.seq.upper()
            c2 = fasta.get(n2)
            aln2 = c2.seq.upper()
            aa = AlignmentAnalyzer(aln1, aln2)
            for pos, b1, b2 in aa.iter_edits():
                ctx = aa.get_context(pos, left=args.ctx_size, right=args.ctx_size+1)
                pos1 = ctx['last_end1']
                pos2 = ctx['last_end2']
                ctx1 = ctx['ctx1']
                ctx2 = ctx['ctx2']
                print (n1, n2, pos, pos1, pos2, b1, b2, ctx_size, ctx1, ctx2, sep='\t')


# TODO realign option
@command.add_sub
@argument('fasta', help='all contigs should be the same lengths')
@argument('-n', '--names', help='comma separated list of names to calculate')
@argument('-r', '--redundant', action='store_true', help='calculate for all the combinations')
def msa2dist(args):
    """ Calculate pairwise edit distance from multiple aligned contigs

        name1 :
        name2 :
        distance : distance
        ins : ins from name1 to name2
        del : del from name1 to name2
        mut : mutation
    """
    with open(args.fasta) as fp:
        fasta = Fasta(fp)
        if args.names:
            names = args.names.split(',')
        else:
            names = fasta.names

        if args.redundant:
            pairs = itertools.product(names, names)
        else:
            pairs = itertools.combinations(names, 2)

        print ('name1', 'name2', 'alen', 'len1', 'len2', 'distance', 'ins', 'del', 'mut', sep='\t')
        for n1, n2 in pairs:
            c1 = fasta.get(n1)
            aln1 = c1.seq.upper()
            c2 = fasta.get(n2)
            aln2 = c2.seq.upper()
            aa = AlignmentAnalyzer(aln1, aln2)

            ins = del_ = mut = 0
            for pos, b1, b2 in aa.iter_edits():
                if b1 == b2:
                    logging.warning('not an edit: (%s, %s, %s)', pos, b1, b2)
                elif b1 == '-':
                    ins += 1
                elif b2 == '-':
                    del_ += 1
                elif b1 != b2:
                    mut += 1
            dist = ins + del_ + mut
            alen = len(aln1)
            len1 = len(aln1.replace('-', ''))
            len2 = len(aln2.replace('-', ''))
            print (n1, n2, alen, len1, len2, dist, ins, del_, mut, sep='\t')


class MSAPileup:
    """
    >>> p1 = MSAPileup(['A', 'A', '*'])
    >>> p1.width, p1.is_variant
    (1, False)
    >>> p1.masks
    (False, False, True)

    >>> p2 = MSAPileup(['A', 'T', '*'])
    >>> p2.width, p2.is_variant
    (1, True)
    >>> p2.masks
    (False, False, True)

    >>> MSAPileup.merge(p1, p2).masks
    (False, False, True)

    >>> p3 = MSAPileup(['*', 'T', 'G'])
    >>> p3.width, p3.is_variant
    (1, True)
    >>> p3.masks
    (True, False, False)

    >>> p4 = MSAPileup.merge(p1, p2, p3)
    >>> p4.width
    3
    >>> p4.bases
    ['AA*', 'ATT', '**G']
    >>> p4.masks
    (False, False, False)
    """

    def __init__(self, bases, pos=None, names=None, mask_char='*'):
        self.bases = bases
        self.width = len(bases[0])
        self.base_counts = Counter(self.bases)
        self.mask_char = mask_char
        self.pos = pos
        self.names = names

    @cached_property
    def masks(self):
        mask_base = self.mask_char * self.width
        return tuple(b == mask_base for b in self.bases)

    @cached_property
    def is_variant(self):
        if self.width > 1:
            warnings.warn('Variant determination for pileup with width > 1 may be inaccurate')
        if self.mask_char * self.width in self.base_counts:
            return len(self.base_counts) > 2
        else:
            return len(self.base_counts) > 1

    @staticmethod
    def merge(*pileups):
        head = pileups[0]
        bases = [''.join(bs) for bs in zip(*(p.bases for p in pileups))]
        return MSAPileup(bases, pos=head.pos, names=head.names, mask_char=head.mask_char)


def iter_msa_pileup(fasta, names=None, variant_only=False, merge_variant=False, mode='dna', mask_char=None, pos_offset=0):
    mask_char = mask_char or get_default_mask_char(mode)
    contigs = fasta.contigs
    if names:
        contig_map = dict((c.name, c) for c in contigs)
        contigs = [contig_map[name] for name in names]
    else:
        names = fasta.names

    #seq_iters = [iter(c.get_clip_mask_seq(mask_char=mask_char), None) for c in contigs]
    seq_iters = [iter(mask_unknown_seq(c.seq.upper(), mask_char=mask_char, mode=mode)) for c in contigs]
    pos = pos_offset  # default is 0 start
    buffer = []  # TODO review is need for merge_variant logic

    while 1:
        bases = [next(it, None) for it in seq_iters]
        #logging.info(bases)
        if bases[0] is None:
            break

        pileup = MSAPileup(bases, pos, names, mask_char=mask_char)
        if merge_variant:
            if (pileup.is_variant and (not buffer or buffer[0].masks == pileup.masks)):  # mask positions are the same
                buffer.append(pileup)
            elif buffer:
                p = MSAPileup.merge(*buffer)
                yield p
                buffer = []
            else:
                buffer.append(pileup)
        elif variant_only:
            if pileup.is_variant:
                assert len(set(pileup.bases) - set([mask_char])) > 1
                yield pileup
        else:
            yield pileup
        pos += 1
    if buffer:
        p = MSAPileup.merge(*buffer)
        yield p


def msa2df(fasta):
    tab = pd.DataFrame()
    for contig in fasta.contigs:
        contig.name
        seq = contig.get_clip_mask_seq()
        tab[contig.name] = pd.Series(list(seq))
    return tab


@command.add_sub
@argument('fasta', help='all contigs should be the same lengths')
@argument('-c', '--color', action='store_true')
@argument('-m', '--mask', action='store_true')
@argument('-u', '--mask-unknown-bases', action='store_true')
@argument('-p', '--show-polymorphic', action='store_true', help='show only positions with difference')
@argument('--vertical', action='store_true')
@argument('-s', '--start', type=int)
@argument('-e', '--end', type=int)
@argument('--mode', choices=['dna', 'aa'], default='dna')
@argument('-S', '--input-stop-char', default=None, help='Stop codon charactor for input')
#@argument('-n', '--ref-name', help='default is consensus')
#@argument('-i', '--ref-index', type=int, help='if specified, ref_name is ignored')   # consensus
def msa_view(args):
    """ Make a consensus sequence from multiple aligned contigs
    """
    mode = args.mode
    input_stop_char = args.input_stop_char
    def decorate(seq, ref=None, columns=None):
        if input_stop_char is not None:
            seq = seq.replace(input_stop_char, STOP_AA)
            if ref is not None:
                ref = ref.replace(input_stop_char, STOP_AA)
        if args.mask_unknown_bases:
            seq = mask_unknown_seq(seq, mask_char='_', mode=mode)
        if columns is not None:
            seq = ''.join(seq[i] for i in columns)
            if ref:
                ref = ''.join(ref[i] for i in columns)
        if args.mask and ref is not None:
            seq = mask_same_bases(seq, ref)
        if args.color:
            if mode == 'dna':
                seq = color_term_dna(seq)
            if mode == 'aa':
                seq = color_term_aa(seq)
        return seq

    with open(args.fasta) as fp:
        columns = None
        fasta = Fasta(fp, load_now=True)
        it = iter(fasta.contigs)
        if args.show_polymorphic:
            args.vertical = True
            columns = [p.pos for p in iter_msa_pileup(fasta, variant_only=True, mode=mode)]
            it = iter(fasta.contigs)

        contig = next(it, None)
        if contig is None:
            logging.warning('No contigs were found.')
            return 1
        ref = contig.seq.upper()
        if (args.start or args.end) and columns is None:
            columns = list(range(len(ref)))
        if args.start:
            columns = list(filter(lambda i: i >= (args.start - 1), columns))
        if args.end:
            columns = list(filter(lambda i: i <= (args.end - 1), columns))

        if args.vertical:
            pos_list = [i for i, _ in enumerate(ref, 1)]
            if columns is not None:
                pos_list = [pos_list[i] for i in columns]
            for coord in pos2text(pos_list, vertical=True, skip_digit=True):
                print (coord)
        else:
            ref1 = ref[args.start:args.end]
            offset = 0 if args.start is None else args.start - 1
            pos_text = get_aln_pos_text(ref1, offset=offset, count_del=True)
            print (pos_text['number'], sep='\t')
            print (pos_text['indicator'], sep='\t')
            pos_text = get_aln_pos_text(ref1, offset=offset, count_del=False)
            print (pos_text['number'], sep='\t')
            print (pos_text['indicator'], sep='\t')

        print (decorate(ref, columns=columns), 0, contig.name_line, sep='\t')
        for idx, contig in enumerate(it, 1):
            seq = contig.seq.upper()
            #pos_text = get_aln_pos_text(seq)
            #print (pos_text['number'], sep='\t')
            #print (pos_text['indicator'], sep='\t')
            #seq = ' ' * contig.lclip + seq[contig.lclip : len(seq) - contig.rclip] + ' ' * contig.rclip
            print (decorate(seq, ref=ref, columns=columns), idx, contig.name_line.replace('\t', ' '), sep='\t')


@command.add_sub
@argument('fasta', help='all contigs should be the same lengths')
@argument('--layout', help='layout file')
@argument('--fig-x', type=float, default=8.3)
@argument('--fig-y', type=float, default=11.7)
@argument('--title', default='MSA layout plot')
@argument('--adjust-right', type=float, default=.8)
def msa_layout_plot(args):
    """
    layout is tsv file
        name:
        show_name:
        type: [gene,exon]
        parent: name
        show_offset:
        show_diff: [01]   # annotate difference
        # show_variant: [01]
        color: e.g. 'red'
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.patches as patches
    matplotlib.rcParams['figure.figsize'] = (args.fig_x, args.fig_y)  # A4
    matplotlib.rcParams['font.size'] = 10.
    matplotlib.rcParams['axes.titlesize'] = 'medium'
    matplotlib.rcParams['axes.labelsize'] = 'medium'
    matplotlib.rcParams['axes.labelpad'] = 8.
    matplotlib.rcParams['xtick.labelsize'] = 'small'
    matplotlib.rcParams['ytick.labelsize'] = 'small'
    matplotlib.rcParams['legend.fontsize'] = 'small'
    tab = pd.read_table(args.layout)
    out = '{0}.pdf'.format(args.layout)
    logging.info('Create %s', out)

    with open(args.fasta) as fp:
        fasta = Fasta(fp)
        name_contigs = dict(zip(fasta.names, fasta.contigs))
        fig = plt.figure()
        #ax = fig.add_subplot(111, aspect='equal')
        ax = plt.gca()
        y_margin = 100
        offset = 100
        height = 50
        y_ids = []
        id_names = {}

        for row in tab.itertuples():
            print (row)
            name = row.name
            if name not in name_contigs:
                logging.warning('%s is not included in MSA', name)
                continue
            parent = row.parent if isinstance(row.parent, int) else 0
            if parent and (parent not in id_names or id_names[parent] not in name_contigs):
                logging.warning('%s is not included in MSA', name)
                continue
            if not parent:
                y_index = len(y_ids)
                id_names[row.id] = name
                y_ids.append(row.id)
            else:
                y_index = y_ids.index(parent)

            color = row.color
            seq = mask_unknown_seq(name_contigs[name].seq)
            loc = SeqLocator(seq)
            end = len(seq)
            y = y_index * y_margin
            for s, e, t in loc.blocks:
                x = s
                width = e - s
                print (x, y, width, height)
                ax.add_patch(patches.Rectangle((x, y), width, height, facecolor=color))
            if not parent:
                ax.text(end + offset, y + y_margin*.5, row.show_name)
            # TODO
            #if args.show_variant:
            #    for var_info in msa.get_variants():
            #        s, e = var_info['start'], var_info['end']
            #        var_info['types'].index(var_info)
            #    for s, e, var_type in ().bases:
            #        x = s
            #        width = e - s
            #        print (x, y, width, height)
            #        ax.add_patch(patches.Rectangle((x, y), width, height, facecolor=color))
        ax.set_aspect('equal')
        ax.set_xlim(0, end + offset)
        ax.set_ylim(y + y_margin , - y_margin)
        ax.get_yaxis().set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.xaxis.set_ticks_position('bottom')
        ax.spines['left'].set_visible(False)
        ax.spines['right'].set_visible(False)
        #ax.set_frame_on(False)
        ax.set_title(args.title)
        plt.subplots_adjust(left=.05, right=args.adjust_right)
        plt.savefig(out, transparent=True)
