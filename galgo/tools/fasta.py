from __future__ import absolute_import
from __future__ import print_function
from collections import namedtuple, defaultdict
from argtools import command, argument
from ..utils import Counter, iter_tabs, safediv
from ..bioseq import Fasta
from ..samutil import parse_region
from .. io import open
import re
import os
import logging


class Contig(object):
    __slots__ = ('name_line', 'name', '_seqs')

    def __init__(self, name_line=''):
        assert not name_line.startswith('>'), 'please remove charcter ">"'
        self.name_line = name_line
        self.name = name_line.split(' ', 2)[0]
        self._seqs = []

    def append(self, seq):
        self._seqs.append(seq)

    def get_seq(self):
        return ''.join(self._seqs)

def iter_contigs(fp):
    contig = None
    try:
        while fp:
            line = next(fp).rstrip()
            if line.startswith('>'):
                if contig is not None:
                    yield contig
                contig = Contig(line[1:])
            else:
                contig.append(line)
    except StopIteration:
        if contig is not None:
            yield contig


@command.add_sub
@argument('fasta', nargs='?', default='/dev/stdin')
@argument('-w', '--width', type=int, default=50)
@argument('-u', '--uppercase', action='store_true')
def fa_fold(args):
    """ Set line breaks with fixed width
    """
    width = args.width
    with open(args.fasta) as fp:
        for contig in iter_contigs(fp):
            print ('>', contig.name_line, sep='')
            seq = contig.get_seq()
            if args.uppercase:
                seq = seq.upper()
            if not width:
                print (seq)
            else:
                reg = len(seq) // width
                last = len(seq) % width
                if reg:
                    for s in xrange(0, reg*width, width):
                        print (seq[s:s + width])
                if last:
                    print (seq[-last:])

@command.add_sub
@argument('fasta', nargs='?', default='/dev/stdin')
@argument('--use-nonn', action='store_true', help='Use Non N bases for gc/at_ratio calculation')
def fa_profile(args):
    """ Analyse contigs
    """
    print ('name', 'length', 'gc_ratio', 'at_ratio', 'nonn_ratio', 'gc_nonn_ratio', 'count_A', 'count_C', 'count_G', 'count_T', 'count_N', sep='\t')
    with open(args.fasta) as fp:
        for contig in iter_contigs(fp):
            seq = contig.get_seq()
            length = len(seq)
            count = Counter(seq.upper())
            cA = count['A']
            cC = count['C']
            cG = count['G']
            cT = count['T']
            cN = count['N']
            cNonN = length - cN
            gc_ratio = 1. * (cC + cG) / length
            at_ratio = 1. * (cA + cT) / length
            nonn_ratio = 1. * cNonN / length
            gc_nonn_ratio = safediv(1. * (cC + cG), cNonN, 0)
            print (contig.name, length,
                   '{0:.3f}'.format(gc_ratio),
                   '{0:.3f}'.format(at_ratio),
                   '{0:.3f}'.format(nonn_ratio),
                   '{0:.3f}'.format(gc_nonn_ratio),
                   cA, cC, cG, cT, cN,
                   sep='\t')


@command.add_sub
@argument('fasta', nargs='?', default='/dev/stdin')
@argument('-m', '--min-run', type=int, default=4)
def fa_stats_hompol(args):
    """ Analyse contigs
    """
    print ('name', 'length', 'min_run',
           'count', 'max_len', 'min_len',
           'len_tot', 'len_AT', 'len_CG', 'len_A', 'len_T', 'len_C', 'len_G',
           'ratio', 'ratio_AT', 'ratio_CG', 'ratio_A', 'ratio_T', 'ratio_C', 'ratio_G',
           sep='\t')

    patterns = dict((base, '%s{%s,}' % (base, args.min_run)) for base in 'ACGT')
    pattern_res = dict((base, re.compile(patterns[base])) for base in 'ACGT')

    def get_ratio(a, b):
        return round(1. * a / b, 3)

    for contig in iter_contigs(open(args.fasta)):
        seq = contig.get_seq().upper()
        length = len(seq)
        base_lens = {}
        for base, reg in pattern_res.items():
            base_lens[base] = [m.end() - m.start() for m in reg.finditer(seq)]

        tot_len = sum(sum(base_lens[base]) for base in 'ACGT')
        count = sum(len(base_lens[base]) for base in 'ACGT')
        max_len = max(max(base_lens[base] or [0]) for base in 'ACGT')
        min_len = min(min(base_lens[base] or [args.min_run]) for base in 'ACGT')
        AT_len = sum(base_lens['A']) + sum(base_lens['T'])
        CG_len = sum(base_lens['C']) + sum(base_lens['G'])

        print (contig.name,
               length,
               args.min_run,
               count,
               max_len,
               min_len,
               tot_len,
               AT_len,
               CG_len,
               sum(base_lens['A']),
               sum(base_lens['T']),
               sum(base_lens['C']),
               sum(base_lens['G']),
               get_ratio(tot_len, length),
               get_ratio(AT_len, length),
               get_ratio(CG_len, length),
               get_ratio(sum(base_lens['A']), length),
               get_ratio(sum(base_lens['T']), length),
               get_ratio(sum(base_lens['C']), length),
               get_ratio(sum(base_lens['G']), length),
               sep='\t')

# deplicated
@command.add_sub
@argument('fasta', nargs='?', default='/dev/stdin')
@argument('--names', nargs='*', help='list of fasta names in order of fasta')
def fa_rename_all(args):
    """ Replace sequence names

    The number of names should be the same as that of sequences.
    """
    names = iter(args.names)

    for line in open(args.fasta):
        if line.startswith('>'):
            tokens = line.rstrip('\r\n').split(' ', 1)
            name = next(names)
            print ('>', name, ' ', ''.join(tokens[1:]), sep='')
        else:
            print (line, end='')  # avoiding doubled return


# TODO make fasta a Fasta instance
def _select_or_omit(fasta, names=None, regexp=None, omit=False, omits=None):
    name_set = set(names or [])
    if regexp:
        pat = re.compile(regexp)
    else:
        pat = None

    def is_target(name):
        matched = bool(name in name_set or (pat and pat.match(name)))
        if matched:
            logging.info('%s is matched. to be omit? %s', name, omit)
        if omit:
            return not matched
        else:
            return matched

    for line in open(fasta):
        if line.startswith('>'):
            logging.debug('%s', line.rstrip())
            tokens = line.rstrip('\r\n').split(' ', 1)
            name = tokens[0][1:]
            rest = ''.join(tokens[1:])
            if is_target(name):
                print ('>', name, ' ', rest, sep='')
                emit = True
            else:
                emit = False
        elif emit:
            print (line, end='')


@command.add_sub
@argument('fasta')
@argument.exclusive(
    argument('-n', '--names', nargs='*', help='list of fasta names to select'),
    argument('-N', '--name-file', help='list of fasta names to select'),
    argument('-r', '--regexp', help='match patterns for name'),
)
def fa_select(args):
    """ Select sequences by names

    O(n) but not required faidx
    """
    if args.name_file:
        names = [line.rstrip() for line in open(args.name_file)]
    else:
        names = args.names
    _select_or_omit(args.fasta, names=names, regexp=args.regexp)


@command.add_sub
@argument('fasta')
@argument.exclusive(
    argument('-n', '--names', nargs='*', help='list of fasta names to omit'),
    argument('-N', '--name-file', help='list of fasta names to omit'),
    argument('-r', '--regexp', help='match patterns for name'),
)
def fa_omit(args):
    """ Omit sequenses by names

    O(n) but not required faidx
    """
    if args.name_file:
        names = [line.rstrip() for line in open(args.name_file)]
    else:
        names = args.names
    _select_or_omit(args.fasta, names=names, regexp=args.regexp, omit=True)


@command.add_sub
@argument('fasta', nargs='?', default='/dev/stdin')
def fa_uniq(args):
    """ Discard duplicated names and show warnings
    """
    have_seen = set()
    emit = True
    with open(args.fasta) as fp:
        for line in fp:
            if line.startswith('>'):
                tokens = line.rstrip('\r\n').split(' ', 1)
                name = tokens[0][1:]
                if name in have_seen:
                    logging.warning('Discarded duplicated entry %s (%s)', name, line.rstrip('\r\n'))
                    emit = False
                else:
                    emit = True
                    print (line, end='')
                have_seen.add(name)
            elif emit:
                print (line, end='')


@command.add_sub
@argument('fasta', nargs='?', default='/dev/stdin')
@argument.exclusive(
    argument('-f', '--from-name'),  # TODO grouping with to_name
    argument('-n', '--names', nargs='*'),
    argument('-N', '--name-file'),
    argument('-M', '--name-map-file', help='tsv file of from-name, to-name maps'),
)
@argument('-t', '--to-name')
def fa_rename(args):
    """ Rename a sequence

    O(n) but not required faidx
    """
    if args.from_name:
        from_name = args.from_name
        to_name = args.to_name

        for line in open(args.fasta):
            if line.startswith('>'):
                name = line.rstrip('\r\n').split(' ', 2)[0][1:]
                if name == from_name:
                    print ('>', to_name, sep='')
                else:
                    print ('>', name, sep='')
            else:
                print (line, end='')
        return 0

    if args.names:
        names = iter(args.names)
        def get_name(current_name):
            return next(names)
    elif args.name_file:
        with open(args.name_file) as fp:
            names = [row[0] for row in iter_tabs(fp)]
            names = iter(names)
        def get_name(current_name):
            return next(names)
    elif args.name_map_file:
        with open(args.name_map_file) as fp:
            name_map = dict((row[0], row[1]) for row in iter_tabs(fp))
        def get_name(current_name):
            return name_map[current_name]
    else:
        raise NotImplementedError

    for line in open(args.fasta):
        if line.startswith('>'):
            tokens = line.rstrip('\r\n').split(' ', 1)
            name = get_name(tokens[0][1:])
            print ('>', name, ' ', ''.join(tokens[1:]), sep='')
        else:
            print (line, end='')  # avoiding doubled return


@command.add_sub
@argument('fasta', nargs='?', default='/dev/stdin')
@argument('--prefix', default='')
@argument('--suffix', default='')
def fa_subname(args):
    """ Substitute fasta names
    """
    prefix = args.prefix
    suffix = args.suffix

    for line in open(args.fasta):
        if line.startswith('>'):
            name = line.rstrip('\r\n').split(' ', 2)[0][1:]
            name = prefix + name + suffix
            print ('>', name, sep='')
        else:
            print (line, end='')


@command.add_sub
@argument('fastq', nargs='?', default='/dev/stdin')
@argument('--prefix', default='')
@argument('--suffix', default='')
def fq_subname(args):
    """ Substitute fastq names
    """
    prefix = args.prefix
    suffix = args.suffix

    for i, line in enumerate(open(args.fastq)):
        if i % 4 == 0:
            name = line.rstrip('\r\n').split(' ', 2)[0][1:]
            name = prefix + name + suffix
            print ('@', name, sep='')
        else:
            print (line, end='')


@command.add_sub
@argument('fasta')
@argument('outdir')
def pb_fa_split(args):
    """ Split fasta to create pacbio input.fofn
    """
    logging.info('Making directory: %s', args.outdir)
    os.makedirs(args.outdir)  # this should be error on found
    run_lines = defaultdict(list)
    with open(args.fasta) as fp:
        for line in fp:
            line = line.rstrip('\r\n')
            if line.startswith('>'):
                run_name = line[1:].split('/')[0]
            run_lines[run_name].append(line)

    # TODO create pool of writers, or buffering
    for run_name in run_lines:
        path = os.path.join(args.outdir, '{0}.fa'.format(run_name))
        with open(path, '+a') as fp:
            for line in run_lines[run_name]:
                print (line, file=fp)


@command.add_sub
@argument('fastq', nargs='?', default='/dev/stdin')
@argument('--prefix', default='')
@argument('--suffix', default='')
@argument('--zmw-start', default=0)
def pb_fq_subname(args):
    """ Substitute fastq names
    """
    prefix = args.prefix
    suffix = args.suffix
    zmw = args.zmw_start
    with open(args.fastq) as fp:
        try:
            while 1:
                line1 = next(fp)
                name = line1.rstrip().split(' ', 2)[0][1:]
                seq = next(fp).rstrip()
                plus = next(fp).rstrip()
                qseq = next(fp).rstrip()
                length = len(seq)
                line1 = '@{}{}{}/{}/{}_{}'.format(prefix, name, suffix, zmw, 0, length)
                print (line1)
                print (seq)
                print (plus)
                print (qseq)
                zmw += 1
        except StopIteration:
            pass

@command.add_sub
@argument('fasta', nargs='?', default='/dev/stdin')
@argument('-u', '--uppercase', action='store_true')
@argument('-H', '--no-header', dest='header', action='store_false')
def fa2tsv(args):
    """ Substitute fasta names
    """
    sep = '\t'
    if args.header:
        print ('name', 'seq', sep=sep)
    with open(args.fasta) as fp:
        for contig in iter_contigs(fp):
            name = contig.name
            seq = contig.get_seq()
            if args.uppercase:
                seq = seq.upper()
            print (name, seq, sep=sep)


@command.add_sub
@argument('fasta', nargs='?', default='/dev/stdin')
@argument('-q', '--quality', default='I')
def fa2fq(args):
    """ Convert fasta to fastq
    """
    bq = args.quality
    class Fastq(object):
        __slots__ = ('name', '_seqs')
        def __init__(self, name=''):
            self.name = name
            self._seqs = []

        def append(self, seq):
            self._seqs.append(seq)

        def emit(self):
            seq = ''.join(self._seqs)
            print ('@' + self.name)
            print (seq)
            print ('+')
            print (bq * len(seq))

    with open(args.fasta) as fp:
        fastq = None
        try:
            while fp:
                line = next(fp).rstrip()
                if line.startswith('>'):
                    if fastq is not None:
                        fastq.emit()
                    name = line[1:]
                    fastq = Fastq(name)
                else:
                    fastq.append(line)
        except StopIteration:
            if fastq is not None:
                fastq.emit()


@command.add_sub
@argument('fastq', nargs='?', default='/dev/stdin')
def fq2fa(args):
    """ Convert fastq to fasta
    """
    with open(args.fastq) as fp:
        try:
            while fp:
                line = next(fp).rstrip()
                rname = line[1:]
                print ('>' + rname)
                seq = next(fp).rstrip() # seq
                print (seq)
                next(fp) # +
                next(fp) # quality
        except StopIteration:
            pass

@command.add_sub
@argument('fasta', default='/dev/stdin')
@argument('-r', '--regions', nargs='+', required=True)
@argument('--mask-char', default='N')
def fa_mask(args):
    """ Mask regions in fasta sequences
    """
    mask_regions = defaultdict(list) # chrom: [(start, end)]
    for region in args.regions:
        c, s, e = parse_region(region)
        mask_regions[c].append((s, e))
        logging.info('Mask region: %s [%s, %s]', c, s, e)
    assert len(args.mask_char) == 1

    with open(args.fasta) as fasta:
        fa = Fasta(fasta)
        for contig in fa:
            seq = contig.seq
            name = contig.name
            print ('>' + contig.name_line)
            for s, e in mask_regions[name]:  # replace all seqs
                logging.info('Masking %s [%s, %s]', name, s, e)
                mask = args.mask_char * (e - s)
                seq = ''.join((seq[:s], mask, seq[e:]))
            print (seq)


@command.add_sub
@argument('fasta', nargs='?', default='/dev/stdin')
@argument('-n', '--names', nargs='+', required=True)
@argument('-e', '--edits', nargs='*', help='sequence of pos:alts')
@argument('--is-compact', action='store_true', default=False, help='input sequence is assumed compact ("-" is omitted)')
def fa_edit(args):
    """ Edit fasta sequence

    pos list (1-based) and alt bases
    ex: 123:AAA-T        (ref)
    ex: 123:AA-TT:AAA-T  ()
    """
    names = args.names
    name_set = set(names)
    contigs = list(iter_contigs(open(args.fasta)))
    logging.info('names: %s, edits: %s', names, args.edits)

    Edit = namedtuple('Edit', 'pos ref alt')
    def parse_edit(edit):
        els = edit.split(':')
        pos = int(els[0])
        if len(els) == 2:
            ref = None
            alt = els[1]
        elif len(els) == 3:
            ref = els[1]
            alt = els[2]
        else:
            raise Exception('Invalid edit {0}'.format(edit))
        if args.is_compact:
            alt = alt.replace('-', '')
            ref = ref and ref.replace('-', '')
        return Edit(pos, ref, alt)

    edits = list(sorted(map(parse_edit, args.edits), key=lambda e: e.pos))

    for contig in contigs:
        print ('>', contig.name_line, sep='')
        seq = contig.get_seq()
        if contig.name in name_set:
            bufs = []
            seq_pos = 1
            for edit in edits:
                bufs.append(seq[seq_pos - 1 : edit.pos - 1])  # prev seq
                bufs.append(edit.alt)
                if edit.ref is None:
                    ref_len = alt_len
                else:
                    ref_len = len(edit.ref)  # assume same length

                seq_pos = edit.pos + ref_len

            if seq_pos - 1 <= len(seq):
                bufs.append(seq[seq_pos - 1 :])

            len0 = len(seq)
            seq = ''.join(bufs)
            logging.info('seq length (%s -> %s)', len0, len(seq))
        print (seq)


_rev_maps = dict(zip('ATCGatcgNn', 'TAGCtagcNn'))

def _revcomp(seq):
    return ''.join(_rev_maps[c] for c in reversed(seq))


@command.add_sub
@argument('fasta', nargs='?', default='/dev/stdin', help='all contigs should be the same lengths')
@argument('-c', '--contig-names', nargs='*', default=None, help='list of reverse strand')
@argument('-a', '--add-suffix', action='store_true', help='add suffix to reversed')
def fa_reverse(args):
    """ Get reverse complements of sequences
    """
    names = set(args.contig_names)

    for contig in iter_contigs(args.fasta):
        reverse = contig.name in names

        name = contig.name
        if args.add_suffix:
            name = contig.name + ('.R' if reverse else '.F')

        print ('>' + name)
        seq = contig.seq
        print (_revcomp(seq) if reverse else seq)


@command.add_sub
@argument('regex')
@argument('fasta', nargs='?', default='/dev/stdin', help='all contigs should be the same lengths')
def fa_grep(args):
    """ Emit matched coordinates and sequence as BED format

    # TODO include context sequence

    Making N bed:
    $ galgo fa_grep 'N+' test.fa | cut -f-3

    Making nonN bed:
    $ galgo fa_grep '[^N]+' test.fa | cut -f-3
    """
    logging.info('pattern: %s', args.regex)
    pattern = re.compile(args.regex)
    for contig in iter_contigs(open(args.fasta)):
        seq = contig.get_seq().upper()
        for m in pattern.finditer(seq):
            print (contig.name, m.start(), m.end(), m.group(), sep='\t')
