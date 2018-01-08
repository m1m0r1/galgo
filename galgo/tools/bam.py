from __future__ import print_function
from argtools import command, argument
import pysam
from pysam import Samfile
import re
import logging
from collections import namedtuple
from builtins import filter, zip, map
from itertools import groupby, chain
from ..utils import fill_text, lfill_text, Counter, cached_property
from ..bioseq import color_term_dna, get_aln_pos_text, dna_revcomp, Fasta, aln2pos, pos2text
from ..samutil import LocalInfo, sam_intervals, Read
from ..samutil import aln2cigar, Cigar, SamInfo
from operator import attrgetter
from tqdm import tqdm


@command.add_sub
@argument('bamfile')
@argument('-l', '--is-list', action='store_true')
@argument('-t', '--tags', help='comma separated tags (e.g. RG)')
def bam_header2tab(args):
    """
    USAGE: bam_haeder2tab {bam_list} -l -t RG
    """
    tag1_set = set(args.tags.split(','))
    #tags = args.tags.split(',')
    fnames = []
    fn_recs = []
    if args.is_list:
        bamfiles = [line.rstrip() for line in open(args.bamfile)]
    else:
        bamfiles = [args.bamfile]
    for bam in bamfiles:
        sam = pysam.Samfile(bam)
        fnames.append(sam.filename)
        fn_recs.append((sam.filename, sam.header))
        sam.close()

    # header
    # {tag1: [{tag2: cont}]}
    def gather_tags(fn_recs):
        """
        Returns: set(['PG', 'RG.ID', 'RG.SM', ...]), {'PG': ..., 'RG.ID': ..., ...})
        """
        tag_set = set()
        out_recs = []
        for fn, recs in fn_recs:
            for tag1, rec1s in recs.items():
                if tag1_set and tag1 not in tag1_set:
                    continue
                for rec1 in rec1s:
                    out_rec = {'fname': fn}
                    if isinstance(rec1, dict):
                        for tag2, rec2 in rec1.items():
                            tag = tag1 + '.' + tag2
                            tag_set.add(tag)
                            out_rec[tag] = rec2
                    else:
                        tag = tag1
                        tag_set.add(tag)
                        out_rec[tag] = rec2
                    out_recs.append(out_rec)
        return tag_set, out_recs

    tag_set, out_recs = gather_tags(fn_recs)
    keys = list(sorted(tag_set)) + ['fname']

    print (*keys, sep='\t')
    for rec in out_recs:
        out = [rec.get(key, '.') for key in keys]
        print (*out, sep='\t')


# This can be used for mapped mate only (so maybe useless).
@command.add_sub
@argument('bamfile')
@argument('regions', nargs='*')
@argument('-o', '--output', default='/dev/stdout')
@argument('--ext-only', action='store_true', help='only emit external mate not found in each region')
#@argument('-u', '--use-unmapped-mate', action='store_true', help='emit unmapped mate')
def bam_with_mate(args):
    """ Emit record with unpaired mate (unmapped mate is ignored)

    - Ignore paring of optical duplicate
    - Ignore paring of secondary and supplementary alignments
    # - Ignore paring if mate has no reference id
    - Ignore pairng if mate has unmapped
    """
    if args.output.endswith('.bam'):
        mode = 'wb'
    else:
        mode = 'wh'
    sam = pysam.Samfile(args.bamfile)
    out = pysam.Samfile(args.output, mode=mode, template=sam)

    use_um = args.use_unmapped_mate

    def is_mating_target(rec):
        return \
           not rec.is_duplicate and \
           not rec.is_secondary and \
           not rec.is_supplementary and \
           not rec.is_unmapped and \
           not rec.mate_is_unmapped

    with sam:
        if args.regions:
            its = (sam.fetch(region=region) for region in args.regions)
        else:
            its = [sam]

        for i, it in enumerate(its):
            region_name = args.regions and args.regions[i]
            wait_recs = {}
            for rec in it:
                if not args.ext_only:
                    out.write(rec)
                if is_mating_target(rec):
                    if rec.qname in wait_recs:
                        del wait_recs[rec.qname]     # mating for this read name is resolved
                    else:
                        wait_recs[rec.qname] = rec   # add this read name to mate waiting list

            logging.info('%s mates will be additionally emited for region %s', len(wait_recs), region_name)
            for rec in sorted(wait_recs.values(), key=lambda r: r.pos):
                try:
                    mrec = sam.mate(rec)
                    if mrec:
                        out.write(mrec)
                except ValueError as e:
                    logging.warning('%s', rec)
                    logging.warning('Mate not found for %s (%s)', rec.qname, e)


@command.add_sub
@argument('bamfile')
@argument('regions', nargs='+')
@argument('-o', '--output', default='/dev/stdout')
@argument('-O', '--only-emit-outside', action='store_true', help='Only emit reads outside of specified regions')
@argument('--search-strategy', choices=['mate_api', 'mate_scan'], default='mate_scan', help='mate_scan is faster')
def bam_extract_region_cover(args):
    """ Emit read pairs that at least one of those are mapped in the specified regions

    This method is to complement regional reads.

    1. search pairs in region (some unmapped reads were contained)
    2. search mapped mates
    3. search unmapped mates

    Following reads are ignored.
    - optical duplicate
    - secondary alignment
    - supplementary alignment

    Usage1:
        $ galgo bam_extract_region_cover all.bam chr1 chr2 -o ext.bam

    Usage2:
        $ samtools view all.bam chr1 chr2 -o ext1.bam
        $ galgo bam_extract_region_cover all.bam chr1 chr2 -o ext2.bam --only-emit-outside
        $ samtools cat -o ext.bam ext1.bam ext2.bam

    TODO need write thread
    """

    if args.output.endswith('.bam'):
        mode = 'wb'
    else:
        mode = 'wh'
    emit_inside = not args.only_emit_outside
    sam = pysam.Samfile(args.bamfile)
    out = pysam.Samfile(args.output, mode=mode, template=sam)

    def is_mating_target(rec):
        return \
           not rec.is_duplicate and \
           not rec.is_secondary and \
           not rec.is_supplementary

    def gen_unmapped(sam):
        last_ref = sam.references[-1]
        last_pos  = sam.lengths[-1] - 1
        next(sam.fetch(last_ref, last_pos), None)
        return sam.fetch(until_eof=True)

    with sam:
        its = (sam.fetch(region=region) for region in args.regions)

        wait_recs = {}
        # 1. search pairs in region (some unmapped reads were contained)
        for i, it in enumerate(its):
            region_name = args.regions and args.regions[i]
            for rec in it:
                if is_mating_target(rec):
                    if rec.qname in wait_recs:   # mating for this read name is resolved
                        mrec = wait_recs.pop(rec.qname)
                        if emit_inside:
                            out.write(rec)   # TODO this part seems a bit slow
                            out.write(mrec)
                    else:
                        wait_recs[rec.qname] = rec   # add this read name to mate waiting list

        # 2. search mapped mates
        cands = (r for r in wait_recs.values() if not r.mate_is_unmapped)
        cands = sorted(cands, key=lambda r: (r.next_reference_id, r.next_reference_start))
        cands = list(cands)

        logging.info('%s mates will be searched from mapped reads', len(cands))
        def iter_mapped_mates1(cands, ref_name):  # using mate api
            for rec in cands:
                mrec = sam.mate(rec)
                yield rec, mrec

        def iter_mapped_mates2(cands, ref_name):  # using full scan
            cand_d = {c.qname: c for c in cands}
            poss = [c.next_reference_start for c in cands]
            start = min(poss)
            end = max(poss) + 1
            logging.info('Searching reads from %s:%s-%s', ref_name, start, end)
            for mrec in sam.fetch(reference=ref_name, start=start, end=end):
                if mrec.qname in cand_d:
                    rec = cand_d.pop(mrec.qname)
                    yield rec, mrec

        # TODO cleverer search strategy is demanded
        if args.search_strategy == 'mate_api':
            iter_mapped_mates = iter_mapped_mates1
        elif args.search_strategy == 'mate_scan':
            iter_mapped_mates = iter_mapped_mates2

        for next_ref_id, cands in groupby(cands, lambda r: r.next_reference_id):
            next_ref_name = sam.getrname(next_ref_id)
            cands = list(cands)
            logging.info('Search mate from %s ... (%s remained)', next_ref_name, len(cands))
            for rec, mrec in tqdm(iter_mapped_mates(cands, next_ref_name)):
                if emit_inside:
                    out.write(rec)
                out.write(mrec)
                del wait_recs[rec.qname]

        # 3. search unmapped mates
        logging.info('%s mates will be searched from unmapped reads', len(wait_recs))
        it = gen_unmapped(sam)
        for mrec in tqdm(it):
            rec = wait_recs.pop(mrec.qname, None)
            if rec:
                if emit_inside:
                    out.write(rec)
                out.write(mrec)

        if wait_recs:
            logging.warning('Mating of %s reads were not resolved.', len(wait_recs))
        else:
            logging.info('Mating of all the reads has been resolved.')


@command.add_sub
@argument('bamfile')
@argument('bamout')
@argument('-u', '--unpaired')
@argument('-f', '--uncount-flag', type=lambda x: int(x, 16), default=0x900, help='uncount-flag for existing reads; default is secondary or supplementary reads')
def bam_discard_unpaired(args):
    sam = pysam.Samfile(args.bamfile)
    read1s = set()
    read2s = set()
    uncount_flag = args.uncount_flag
    for rec in sam:
        if rec.flag & uncount_flag:
            continue
        if rec.is_read1:
            read1s.add(rec.qname)
        else:
            read2s.add(rec.qname)

    paired = read1s & read2s
    unpaired = (read1s - paired) | (read2s - paired)
    logging.info('Paired reads: %s', len(paired))
    logging.info('Unpaired reads: %s', len(unpaired))

    sam.reset()
    out = pysam.Samfile(args.bamout, mode='wb', template=sam)
    out2 = pysam.Samfile(args.unpaired, mode='wb', template=sam) if args.unpaired else None
    for rec in sam:
        if rec.qname in paired:
            # workaround for avoiding duplicated emission of the same record
            if rec.is_read1:
                if rec.qname in read1s:
                    out.write(rec)
                    read1s.remove(rec.qname)
            else:
                if rec.qname in read2s:
                    out.write(rec)
                    read2s.remove(rec.qname)
        elif out2:
            out2.write(rec)


@command.add_sub
@argument('bamfile')
@argument('-l', '--list')
@argument('--max-edit', type=int, help='use NM tag')
@argument('--max-overhang', type=int, help='set limit to overhang length between read and reference')
@argument('-o', '--output', default='/dev/stdout')
def bam_filter(args):
    sam = pysam.Samfile(args.bamfile)
    if args.output.endswith('.bam'):
        mode = 'wb'
    else:
        mode = 'wh'
    out = pysam.Samfile(args.output, mode=mode, template=sam)
    max_edit = args.max_edit
    max_overhang = args.max_overhang

    saminfo = SamInfo(sam)
    it = map(saminfo.get_read_info, sam)
    if args.list:
        white_list = set([name.rstrip() for name in open(args.list)])
        cond = lambda read: read.rec.qname in white_list
        it = filter(cond, it)

    if args.max_edit is not None:
        cond = lambda read: read.edit <= max_edit
        it = filter(cond, it)

    if args.max_overhang is not None:
        cond = lambda read: read.overhang <= max_overhang
        it = filter(cond, it)

    for read in it:  # fetch only primary and non-supplementary reads
        out.write(read.rec)


@command.add_sub
@argument('bamfile')
@argument('-r', '--region')
@argument('-m', '--with-mate')
def bam_profile(args):
    """
    """
    sam = pysam.Samfile(args.bamfile)
    attrs = ['name', 'read1', 'unmapped', 'suppl', 'num_frags', 'qlen', 'rname', 'start', 'end', 'alen', 'mapq', 'reverse', 'lclip', 'rclip', 'edit', 'nins', 'ndel']
    if args.with_mate:
        attrs.extend('mate_miss mate_rname mate_pos mate_invert mate_back'.split(' '))
    getter = attrgetter(*attrs)
    with sam:
        print (*(attrs + ['tags_md']), sep='\t')
        for rec in sam.fetch(args.region):
            r = Read(rec)
            print (*(getter(r) + (dict((k, v) for k, v in rec.tags if k != 'MD'),)), sep='\t')


# TODO sort by read_name, read_length, ref_identity, ref_diversity <= weighted identity (using base probability as weight is ok)
# show read properties (mismatch on reads
@command.add_sub
@argument('bam')
@argument('region', nargs='+')  # TODO currently, need this option
@argument('-r', '--reference')
@argument('-c', '--color-base', action='store_true')
@argument('-m', '--mask-ref-base', action='store_true')
@argument('--order', help='sort orders: select from qname, edit, nins, ndel, alen, <ref_pos>'.format())
@argument('-R', '--reverse', action='store_true')
@argument('--vertical', action='store_true')
@argument('-p', '--show-polymorphic', action='store_true')
def bam_aln_view(args):
    """
    """
    sam = Samfile(args.bam)
    skip_flag = 0x004
    sep = ' ' * 10
    sep = '\t'
    def decorate(aln, ref_aln=None, indicator=None, columns=None):
        if ref_aln and args.mask_ref_base:
            #aln = ''.join(('.' if a == r and a != ('-' or ' ') else a) for a, r in zip(aln, ref_aln))
            aln = ''.join(('.' if a == r else a) for a, r in zip(aln, ref_aln))  # '-' is also shown as '.'
        if columns:
            aln = ''.join(aln[i] for i in columns)
            if indicator:
                indicator = ''.join(indicator[i] for i in columns)
        if indicator:
            # if aln is blank, fill with indicator instead
            aln = ''.join(i if a == ' ' else a for a, i in zip(aln, indicator))
        if args.color_base:
            aln = color_term_dna(aln)
        return aln
    if args.show_polymorphic:
        args.vertical = True

    def get_sorted(read_alns):
        if args.order is None:
            return read_alns
        if args.order.isdigit():
            pos = int(args.order) - 1
            lpos = pos - loc.start
            return sorted(read_alns, key=lambda x: (x[0].end >= pos, x[0].start > pos, x[1][lpos], x[0].start, x[0].end))  # reverse
        order = args.order
        assert order in 'qname edit edit_ratio nins ndel alen'.split(' ')
        return sorted(read_alns, key=lambda x: getattr(x[0], order), reverse=args.reverse)

    logging.info('Target region: %s', args.region)
    with sam:
        for iv in sam_intervals(sam, regions=(args.region or None)):
            loc = LocalInfo(sam, iv, fasta=args.reference, skip_flags=skip_flag)
            columns = None
            if args.show_polymorphic:
                # show polymorphic columns
                columns = []
                for i, info in enumerate(loc.iter_column_info()):
                    #if len(set(info['bases']) - set(['N', ' '])) == 1:
                    #    columns.append(i)
                    if len(set(info['bases']) - set(['N', ' '])) > 1:
                        columns.append(i)
            print ('start', loc.start)
            print ('end', loc.end)
            print ('align length', loc.align_length)
            ref_aln = ''.join(loc.get_ref_align())
            ref_dec = decorate(ref_aln, columns=columns)
            indicator = None
            # show vertical pos text
            if args.vertical:
                pos_list = aln2pos(ref_aln, offset=loc.start+1)
                if columns:
                    pos_list = [pos_list[i] for i in columns]
                for coord in pos2text(pos_list, vertical=True, skip_digit=True):
                    #if columns:
                    #    coord = ''.join(coord[i] for i in columns)
                    print (coord)
            else:
                pos_txt = get_aln_pos_text(ref_aln, offset=loc.start)
                indicator = pos_txt['indicator']
                print (pos_txt['number'])
                print (indicator)

                indicator = indicator.replace('|', ':')  # for visibility

                # if data is sorted by bases, show the sorted position
                if args.order and args.order.isdigit():
                    order_pos = int(args.order) - 1
                    lpos1 = loc.get_left_offset(order_pos)
                    #lpos1 = order_pos - loc.start #+ 1
                    if indicator[lpos1] in ' :':
                        indicator = indicator[:lpos1] + '|' + indicator[lpos1+1:]

            print (ref_dec, iv.contig, 'mapq', 'AS', 'clip', 'edit', 'nins', 'ndel', 'alen', 'edit_ratio', sep=sep)
            it = loc.iter_read_aligns()
            it = get_sorted(it)
            for read, aln in it:
                read_aln = decorate(''.join(aln), ref_aln, indicator=indicator, columns=columns)
                clip_status = ('1' if read.has_lclip else '0') + ('1' if read.has_rclip else '0')
                print (read_aln, '',
                       read.mapq,
                       read.tags.get('AS', '.'),
                       clip_status, read.edit, read.nins, read.ndel, read.alen,
                       '{0:.2f}'.format(read.edit_ratio),
                       read.rec.qname,
                       sep=sep)



class SamConcat(object):
    # TODO concat other headers

    def __init__(self, bams):
        self._bams = bams
        self._refs = []
        self._lengths = []
        self._refset = set()

        for fn in self._bams:
            with Samfile(fn) as sam:
                for ref, length in zip(sam.references, sam.lengths):
                    if ref not in self._refset:
                        self._refs.append(ref)
                        self._lengths.append(length)
                        self._refset.add(ref)

        logging.info('Ref lengths: %s', len(self._refs))

    def write(self, output):
        if output.endswith('.bam'):
            mode = 'wb'
        else:
            mode = 'wh'
        out = pysam.Samfile(output, mode=mode, reference_names=self._refs, reference_lengths=self._lengths)

        for fn in self._bams:
            with Samfile(fn) as sam:
                tid_map = {sam.gettid(ref): out.gettid(ref) for ref in sam.references}  # {org_tid: out_tid}
                for rec in sam:
                    rec.reference_id = tid_map.get(rec.reference_id, -1)
                    rec.next_reference_id = tid_map.get(rec.reference_id, -1)
                    out.write(rec)



# TODO rescue other header contents (e.g. @RG)
# TODO fai file for set reference order
@command.add_sub
@argument('bams', nargs='+')
@argument('-o', '--output', default='/dev/stdout')
def bam_cat(args):
    """ Concatenate bamfiles in order with different references
    """
    concat = SamConcat(args.bams)
    concat.write(args.output)


@command.add_sub
@argument('bam')
@argument('-r', '--region', help='region of target bam file')
@argument('-s', '--source-bam')  # TODO pair of fastq
@argument('-o', '--output', default='/dev/stdout')
def bam_fill_seq(args):
    """ Fill empty sequence with known seqs
    """
    if not args.source_bam:
        source_bam = args.bam
    else:
        source_bam = args.source_bam
    logging.info('Loading samfile: %s', source_bam)
    src_seqs = {1: {}, 2: {}}

    src = pysam.Samfile(source_bam)
    with src:
        for rec in src:
            if rec.is_supplementary:  # skip supplementary alignment
                continue
            if rec.is_secondary:  # skip supplementary alignment
                continue
            if rec.query_sequence is None:  # empty
                continue
            if rec.is_read2:
                src_seqs[2][rec.qname] = (rec.query_sequence, rec.query_qualities, rec.is_reverse)
            else:
                src_seqs[1][rec.qname] = (rec.query_sequence, rec.query_qualities, rec.is_reverse)

    logging.info('Loaded read1 : %s', len(src_seqs[1]))
    logging.info('Loaded read2 : %s', len(src_seqs[2]))

    sam = Samfile(args.bam)
    if args.output.endswith('.bam'):
        mode = 'wb'
    else:
        mode = 'wh'
    out = pysam.Samfile(args.output, mode=mode, template=sam)

    if args.region:
        it = sam.fetch(region=args.region)
    else:
        it = sam

    for rec in it:
        qname = rec.qname
        if rec.query_sequence is None:  # only fill when empty
            ret = src_seqs[2 if rec.is_read2 else 1].get(rec.qname)
            if ret is not None:
                seq, qual, is_rev = ret
                if is_rev != rec.is_reverse:
                    seq = dna_revcomp(seq)
                    if qual is not None:
                        qual = list(reversed(qual))
                cigar = Cigar(rec.cigartuples)
                seq = cigar.hard_clip_seq(seq)
                if qual is not None:
                    qual = cigar.hard_clip_seq(qual)
                rec.query_sequence = seq  # refill
                rec.query_qualities = qual

        out.write(rec)


@command.add_sub
@argument('bam')
@argument('-r', '--region', help='region of target bam file')
@argument('--max-nm', default=4, type=int)
@argument('--max-depth', default=8000, type=int)
@argument.exclusive(
    argument('--summary', action='store_true'),
    argument('--read-count', action='store_true'),
)
def bam_depth_with_nm(args):
    """
    * unmapped is discarded
    * both clipped is discarded
    * end clipped is included
    * multimap is included
    * stratified with NM

    default mode:
        pos is 1-based

    summary mode:
        covered
    """
    sam = Samfile(args.bam)
    if args.region:
        c, s, e = parse_region(args.region)
        it = sam.pileup(reference=r, start=s, end=e, max_depth=args.max_depth)
    else:
        it = sam.pileup(max_depth=args.max_depth)
    sam_info = SamInfo(sam)

    def cond(prec):
        rec = prec.alignment
        if rec.is_unmapped:
            return False
        read = sam_info.get_read_info(rec)
        if read.overhang > 0:
            return False
        return True

    max_key = 'NM_more'
    nm_keys = ['NM' + str(nm) for nm in range(args.max_nm+1)] + [max_key]
    def get_key(prec):
        rec = prec.alignment
        nm = rec.get_tag('NM')
        if nm < args.max_nm:
            return 'NM' + str(nm)
        else:
            return max_key

    header = ['contig', 'pos'] + nm_keys
    def iter_table(it):
        Record = namedtuple('Record', header)
        for pcol in it:
            ps = filter(cond, pcol.pileups)
            counts = Counter(map(get_key, ps))
            yield Record(pcol.reference_name, pcol.pos+1, *(counts[k] for k in nm_keys))

    summary_header = ['contig', 'length', 'covered'] + nm_keys
    def iter_summary(it):
        """ NMx is the number of covered position with at least a read whose edit distance to the refernece is under x.
        """
        Record = namedtuple('Record', summary_header)
        def get_min_nm(row):
            for k in nm_keys:
                if getattr(row, k) > 0:
                    return k

        it1 = iter_table(it)
        for contig, rows in groupby(it1, lambda row: row.contig):
            length = sam_info.get_length(contig)
            counts = Counter([get_min_nm(row) for row in rows])
            nm_counts = [counts[k] for k in nm_keys]
            covered = sum(nm_counts)
            yield Record(contig, length, covered, *nm_counts)

    read_count_header = ['contig', 'length', 'total'] + nm_keys
    def iter_read_counts(it):
        """ NMx is the number of reads whose edit distance to the refernece is under x.
        """
        Record = namedtuple('Record', read_count_header)

        it1 = iter_table(it)
        for contig, rows in groupby(it1, lambda row: row.contig):
            length = sam_info.get_length(contig)
            rows = list(rows)
            counts = {}
            for k in nm_keys:
                counts[k] = sum(getattr(row, k) for row in rows)

            nm_counts = [counts[k] for k in nm_keys]
            total = sum(nm_counts)
            yield Record(contig, length, total, *nm_counts)

    if args.summary:
        logging.info('Emit coverage summary')
        print (*summary_header, sep='\t')
        for row in iter_summary(it):
            print (*row, sep='\t')
    elif args.read_count:
        logging.info('Emit read counts')
        print (*read_count_header, sep='\t')
        for row in iter_read_counts(it):
            print (*row, sep='\t')
    else:
        print (*header, sep='\t')  # header
        for row in iter_table(it):
            print (*row, sep='\t')


_MappingBlock = namedtuple('MappingBlock', 'start end is_del local_start local_end')
class MappingConverter(object):
    '''
    Illustration with examples

    pos: 0-based coordinate

                         1         2         3         4         5         6         7
    pos        0123456789012345678901234567890123456789012345678901234567890123456789012345
                                                          1
    local_pos            0123 456                     7 8901         234
                                                           1
    last_end             01234456                     677890111111111123

    >>> aln = '----------AGGG-GTT---------------------C-CCTA---------TTA-------------------'

                         AGGG GTT  # (0, '7M') => (10, '4M1D3M')
                         AGGG GTT  # (0, '4H7M') => (10, '4H4M1D3M')
                     AAAAAGGG GTT  # (0, '4S7M') => (10, '4S4M1D3M')
                     AAAAAGG- GTT  # (0, '4S3M1D3M') => (10, '4S3M2D3M')
                         AGGG GTT  # (0, '3M4I4M') => (10, '3M4I1M1D3M')
                           +  # insertion at next
                         AGGG GTT  # (0, '4M4I3M') => (10, '4M4I1D3M')
                            +
                          GGG GTT  # (1, '3M4I1M2I2M') => (11, '3M4I1D1M2I2M')
                            + +
                          GGG GTT                     C -CTA         TC    # (1, 5M4I2M1D5M) => (11, '3M1D2M4I1M21D1M2D3M9D2M4H')
                               +
                          GGG GTT                     C -CTA         TC    # (1, 30H5M4I2M3I1D5M4H) => (11, '30H3M1D2M4I1M21D1M3I2D3M9D2M4H')
                               +                      +
                          GGG GTT                     C -CTA         TCATTTT    # (1, 30H5M4I2M3I1D6M4S) => (11, '30H3M1D2M4I1M21D1M3I2D3M9D3M4S')
                               +                      +

    >>> mc = MappingConverter(aln)

    >>> lposs = [0, 3, 4, 7, 11, 15]
    >>> [mc.get_pos(lpos) for lpos in lposs]
    [10, 13, 15, 39, 44, -1]

    >>> mc.get_ref_cigar(3, 6).to_str()
    '3D'
    >>> mc.get_ref_cigar(3, 17).to_str()
    '7D4M1D2M'

    >>> x, y = (mc.convert(0, Cigar.parse('7M'))); x, y.to_str()
    (10, '4M1D3M')
    >>> x, y = (mc.convert(0, Cigar.parse('4H7M'))); x, y.to_str()
    (10, '4H4M1D3M')
    >>> x, y = (mc.convert(0, Cigar.parse('4S7M'))); x, y.to_str()
    (10, '4S4M1D3M')
    >>> x, y = (mc.convert(0, Cigar.parse('4S3M1D3M'))); x, y.to_str()
    (10, '4S3M2D3M')
    >>> x, y = (mc.convert(0, Cigar.parse('3M4I4M'))); x, y.to_str()
    (10, '3M4I1M1D3M')
    >>> x, y = (mc.convert(0, Cigar.parse('4M4I3M'))); x, y.to_str()
    (10, '4M4I1D3M')
    >>> x, y = (mc.convert(1, Cigar.parse('3M4I1M2I2M'))); x, y.to_str()
    (11, '3M4I1D1M2I2M')
    >>> x, y = (mc.convert(1, Cigar.parse('5M4I2M1D5M4H'))); x, y.to_str()
    (11, '3M1D2M4I1M21D1M2D3M9D2M4H')
    >>> x, y = (mc.convert(1, Cigar.parse('30H5M4I2M3I1D5M4H'))); x, y.to_str()
    (11, '30H3M1D2M4I1M21D1M3I2D3M9D2M4H')
    >>> x, y = (mc.convert(1, Cigar.parse('30H5M4I2M3I1D6M4S'))); x, y.to_str()
    (11, '30H3M1D2M4I1M21D1M3I2D3M9D3M4S')

                         1         2         3         4         5         6         7         8         9        10        11        12        13        14        15
               0....5....0....5....0....5.7..0....5...90....5....0....5.7..0...45....01...5....0....5....0....5....0....5....0....5....0....5....0....5....0....5....0
    >>> aln = 'TTTATTTATTTATTTATTTAT-----TTTTTTAAGATGGA-----------------GTCTCGCTTTGTTGC---------------CCAGGCTGGAGTGCAG--------------------------------------TGGCGTGATC'

               ....................-     -.............                 ...............               ................                                      ..    # 6H20M2D21M1I25M
                                                                               +1
                        20M         7D        13M           17D          8M    1I  7M       15D             16M                        38D                   2M

    >>> mc = MappingConverter(aln)
    >>> x, y = mc.convert(0, Cigar.parse('6H20M2D21M1I25M')); x, ' '.join(y.to_str_tuple())
    (0, '6H 20M 7D 13M 17D 8M 1I 7M 15D 16M 38D 2M')

    '''
    def __init__(self, ref_aln):
        """ alignment of reference sequence to objective coordinate
        """
        self._blocks = []
        s = 0
        local_s = 0
        for op, l in aln2cigar(ref_aln):
            local_l = l if op == Cigar.M else 0
            b = _MappingBlock(s, s + l, op == Cigar.D, local_start=local_s, local_end=local_s + local_l)
            self._blocks.append(b)
            s += l
            local_s += local_l

    def get_pos(self, lpos):  # TODO binary search?
        for b in self._blocks:
            if b.is_del:
                continue
            if b.local_end <= lpos:
                continue
            offset = (lpos - b.local_start)
            return b.start + offset
        return -1

    def get_ref_cigar(self, start, end):
        """
        start: aln coordinate
        end: aln coodrinate
        """
        cigar = Cigar()
        for b in self._blocks:
            if b.end < start:
                continue
            s = max(start, b.start)
            e = min(end, b.end)
            op = Cigar.D if b.is_del else Cigar.M
            cigar.add((op, e - s))
            if end < b.end:
                break
        return cigar

    def convert(self, lpos, cigar):
        poss = []
        out_cigar = Cigar()
        cur_lpos = lpos   # reference local position (not a position on sequence!)
        ib = iter(self._blocks)
        block = next(ib)
        while block.local_end <= cur_lpos:  # skip blocks
            block = next(ib, None)

        for op, l in cigar:
            if op != Cigar.I:  # emit past blocks if operation is not an I (which should be emited just after previous cigar operation)
                while block and block.local_end < cur_lpos:
                    if block.is_del:
                        out_cigar.add((Cigar.D, block.end - block.start))
                    else:
                        out_cigar.add((Cigar.M, block.end - block.start))
                    block = next(ib, None)
            if op in (Cigar.S, Cigar.H):
                out_cigar.add((op, l))
                continue
            if op in (Cigar.M, Cigar.D):
                if block:
                    poss.append((block.start + (cur_lpos - block.local_start)))  # offset at cur_lpos
                    cur_lend = cur_lpos + l
                while block and block.local_start < cur_lend:
                    if block.is_del:
                        out_cigar.add((Cigar.D, block.end - block.start))
                    else:
                        length = min(block.local_end, cur_lend) - max(block.local_start, cur_lpos)
                        if length:
                            #out_cigar.add((Cigar.M, length))
                            out_cigar.add((op, length))
                    if block.local_end <= cur_lend:  # when this block is end, read next
                        block = next(ib, None)
                    else:
                        break
                cur_lpos += l
                continue
            if op == Cigar.I:
                out_cigar.add((Cigar.I, l))
                continue
            raise NotImplementedError

        return (poss[0], out_cigar)


# Warning when incompatible
class _CigarChecker:
    def __init__(self):
        self._dat = {}

    def check(self, rec, pos, cigar, org_cigar, mc, q_aln):
        end = mc.get_pos(rec.pos + cigar.ref_length)
        end = None if end == -1 else end
        if (rec.qname, pos) in self._dat and self._dat[rec.qname, pos]['cigar'] != cigar:
            old = self._dat[rec.qname, pos]
            logging.warning('Cigar differed for %s, %s', rec.qname, pos)
            logging.warning('org %s %s %s', rec.pos, org_cigar, org_cigar.query_length)
            logging.warning('new %s %s %s', pos, cigar, cigar.query_length)
            logging.warning('old %s %s %s', pos, old['cigar'], old['cigar'].query_length)
            logging.warning('read %s', rec.seq)
            logging.warning('ref     %s-%s %s', pos, end, q_aln.seq[pos:end])
            logging.warning('ref_old %s-%s %s', pos, old['end'], old['seq'])
        self._dat[rec.qname, pos] = {'end': end, 'seq': q_aln.seq[pos:end], 'cigar': cigar}

#TODO
# - take multiple fasta and names
# - set mate information for given msas
@command.add_sub
@argument('bam')
@argument('-f', '--msa-fastas', required=True, nargs='+')
@argument('-n', '--refnames', nargs='+')
@argument('-o', '--output', default='/dev/stdout')
@argument('--check', action='store_true')
@argument('--skip-flag', type=lambda x: int(x, 16), default=0x000)
@argument('--keep-rest', action='store_true', help='keep rest of reference specified in msa_fastas')
def bam_surject_msa(args):
    """
    Caveats:
    - flags are remained as original statuses
    - remaining original values for MD, NM, and AS tags
    - mate are given as unmapped
    - same records are emited
    """
    skip_flag = args.skip_flag
    sam = Samfile(args.bam)
    fasta = Fasta(open(args.msa_fasta))
    mapped_ref_set = set(sam.references)

    # setup output
    if args.refnames is None:
        refnames = ['consensus{0}'.format(i) for i in xrange(len(args.msa_fastas))]
    else:
        refnames = args.refnames
    assert len(refnames) == len(args.msa_fastas), 'The number of refnames should be the same as that of msa_fastas.'

    logging.info('Loading MSA fastas')
    logging.info('Skip flag: %s', args.skip_flag)
    fastas = []
    ref_lens = []
    target_ref_set = set()
    for fn in args.msa_fastas:
        with open(fn) as fp:
            fasta = Fasta(fp)
            fastas.append(fasta)
            if len(fasta.contigs) == 0:
                logging.error('Fasta file %s has no contigs', fn)
                raise Exception('No contigs')
            ref_lens.append(len(fasta.contigs[0]))
            target_ref_set.update(fasta.names)

    rest_refs = [r for r in sam.references if r not in target_ref_set]
    logging.info('%s are included in surjection targets.', len(target_ref_set))
    logging.info('%s are not included in surjection targets.', len(rest_refs))
    if args.keep_rest:
        logging.info('Rest of reference will be kept in surjected BAM file')
        org_ref_len_map = dict(zip(sam.references, sam.lengths))
        refnames.extend([r for r in rest_refs])
        ref_lens.extend([org_ref_len_map[r] for r in rest_refs])
        fastas.extend([None for r in rest_refs])

    logging.info('Setting up output BAMs')
    if args.output.endswith('.bam'):
        mode = 'wb'
    else:
        mode = 'wh'
    out = pysam.Samfile(args.output, mode=mode, reference_names=[refname], reference_lengths=[ref_length])

    # iteration
    for refname, fasta in zip(refnames, fastas):
        out_tid = out.gettid(refname)
        if fasta is None:
            logging.info('Transfering %s', refname)
            src_tid = sam.gettid(refname)
            for rec in sam.fetch(reference=refname):
                if rec.flag & skip_flag:
                    continue
                a = rec.__copy__()
                a.reference_id = out_tid
                if a.next_reference_id != src_tid:  # pair on the same refs
                    a.next_reference_id = out_tid
                else:
                    a.next_reference_id = -1     # unpair
                    a.next_reference_start = -1
                out.write(a)
            continue
        logging.info('Surjecing to %s', refname)
        query_refs = fasta.names
        cc = _CigarChecker() if args.check else None
        for qref in query_refs:
            if qref not in mapped_ref_set:
                logging.warning('%s is not found in original BAM file', qref)
                continue
            #a = pysam.AlignedSegment()
            a = rec.__copy__()
            #print (rec)
            if not rec.is_unmapped:
                org_cigar = Cigar(rec.cigartuples)
                pos, cigar = mc.convert(rec.pos, org_cigar)
                if org_cigar.query_length != cigar.query_length:
                    logging.error('Invalid cigar conversion for %s', rec.qname)
                    logging.error('org %s %s %s', rec.pos, org_cigar, org_cigar.query_length)
                    logging.error('new %s %s %s', pos, cigar, cigar.query_length)
                    s1 = pos
                    e1 = mc.get_pos(rec.pos + cigar.ref_length)
                    logging.error('ref %s-%s %s', s1, e1, mc.get_ref_cigar(s1, e1))
                    logging.error('read %s', rec.seq)
                    logging.error('qref %s', q_aln.seq[s1:e1])
                    raise Exception('Incompatible Cigar')
                cc and cc.check(rec, pos, cigar, org_cigar, mc, q_aln)
                a.cigar = cigar.values
                a.reference_start = pos
            a.reference_id = out_tid
            a.next_reference_id = -1     # this is required
            a.next_reference_start = -1  # this is required
            #a.flag = rec.flag
            #orec.seq = '*'
            #print (orec)
            out.write(a)


@command.add_sub
@argument('bam')
@argument('-o', '--output', default='/dev/stdout')
def bam_uniq(args):
    """
    * BAM file should be sorted in (tid, pos)
    * (qname, pos, is_unmapped, is_read_2, cigar) is checked
    * if multiple records exist, primary alignment is selected
    * scores are not changed
    """
    sam = Samfile(args.bam)

    # setup output
    if args.output.endswith('.bam'):
        mode = 'wb'
    else:
        mode = 'wh'
    out = pysam.Samfile(args.output, mode=mode, template=sam)
    it = sam  # TODO region

    def get_key(rec):
        return (rec.qname, rec.pos, rec.is_unmapped, rec.is_read2, rec.cigar)

    def get_best_rec(recs):
        for rec in recs:
            if not rec.is_secondary and not rec.is_supplementary:
                return rec
        return rec  # No primary alignments were found

    for (tid, pos), recs in groupby(it, lambda rec: (rec.tid, rec.pos)):  # assume position sorted
        recs1 = sorted(recs, key=get_key)   # manual sort by key is needed
        for key, recs2 in groupby(recs1, get_key):
            rec = get_best_rec(recs2)
            out.write(rec)


@command.add_sub
@argument('bam')
@argument('-o', '--output', default='/dev/stdout')
def bam_add_fraginfo(args):
    """ Add information of mapped fragments

    * BAM file should be sorted by qname
    * This method does not assume an unmapped record which has mapping candidates of another references

    -----------------------------------------------------------------------------------------
    Tag  Type  Description
    -----------------------------------------------------------------------------------------
     # XI   Z     candidate id for multimap alternatives
     XC   Z     List of alternative candidates (RNAME,POS,STRAND,CIGAR,IS_SUPPLEMENTARY;...)
     YC   i     Number of alternative candidates
     XM   Z     List of mate candidates (RNAME,POS,STRAND,CIGAR,IS_SUPPLEMENTARY;...)
     YM   i     Number of mate candidates
    -----------------------------------------------------------------------------------------
    """

    if args.output.endswith('.bam'):
        mode = 'wb'
    else:
        mode = 'wh'

    def get_map_info(rec):
        rname = sam.getrname(rec.tid)
        strand = '+-'[rec.is_reverse]
        is_sup = '01'[rec.is_supplementary]
        return '{},{},{},{},{}'.format(rname, rec.pos, strand, rec.cigarstring, is_sup)

    def annotate(recs):
        recs = list(recs)
        alt_cands = []
        mate_cands = []
        same_ref_segs = []
        part_recs = {1: [], 2: []}

        if len(recs) < 2:  # no annotation required
            return

        for rec in recs:
            part = 2 if rec.is_read2 else 1
            info = get_map_info(rec)
            part_recs[part].append((rec, info))

        if len(part_recs[1]) < 2 and len(part_recs[2]) < 2:  # no annotation required
            return

        # annotate here
        for i, rec in enumerate(recs, 1):
            part = 2 if rec.is_read2 else 1
            mate_part = 2 if part == 1 else 1
            cand_infos = [info for r, info in part_recs[part] if r != rec]  # avoid to add self
            mate_infos = [info for r, info in part_recs[mate_part]]
            tags = rec.get_tags() + [
                # ('XI', i, 'Z'),
                ('XC', ';'.join(cand_infos), 'Z'),
                ('YC', len(cand_infos), 'i'),
                ('XM', ';'.join(mate_infos), 'Z'),
                ('YM', len(mate_infos), 'i'),
            ]
            rec.set_tags(tags)

    with Samfile(args.bam) as sam:
        out = pysam.Samfile(args.output, mode=mode, template=sam)
        it = sam
        for qname, recs in groupby(it, lambda x: x.qname):
            try:
                recs = list(recs)
                annotate(recs)
                for rec in recs:
                    out.write(rec)
            except Exception as e:
                logging.error('Error for %s', qname)
                for rec in recs:
                    logging.error('%s', str(rec))
                raise
