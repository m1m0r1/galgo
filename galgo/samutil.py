# -*- coding: utf-8 -*-
import logging
from operator import attrgetter
from .interval import Interval, IntervalMixin
import pysam
from .utils import fill_text, cached_property, attrdict
from . import csamutil
from builtins import zip
import re

# deprecated
def _parse_region_old(region):
    """
    Args:
        region: 1-based region
    Returns:
        (contig, start, end) 0-based coordinate
    """
    sp = region.split(':')
    contig = sp[0]
    if len(sp) == 1:
        return (contig, None, None)
    sp = sp[1].split('-')
    start = int(sp[0].replace(',', ''))
    if len(sp) == 1:
        return (contig, start - 1, None)
    return (contig, start - 1, int(sp[1].replace(',', '')))


_region_pat2 = re.compile('^([\w*:]+):(\d+)-(\d+)$')
_region_pat1 = re.compile('^([\w*:]+):(\d+)$')
_region_pat0 = re.compile('^([\w*:]+)$')

def parse_region(region):
    """
    Args:
        region: 1-based region
    Returns:
        (contig, start, end) 0-based coordinate

    >>> parse_region('abc*011')
    ('abc*011', None, None)
    >>> parse_region('abc:abc00001')
    ('abc:abc00001', None, None)
    >>> parse_region('abc:abc00001:1001')
    ('abc:abc00001', 1000, None)
    >>> parse_region('abc00001:1001')
    ('abc00001', 1000, None)
    >>> parse_region('abc:abc00001:1001-2000')
    ('abc:abc00001', 1000, 2000)
    >>> parse_region('1:1001-2000')
    ('1', 1000, 2000)
    >>> parse_region('1:1001')
    ('1', 1000, None)
    >>> parse_region('1')
    ('1', None, None)
    """
    m = _region_pat2.match(region)
    if m:
        c, s, e = m.groups()
        return (c, int(s) - 1, int(e))
    m = _region_pat1.match(region)
    if m:
        c, s = m.groups()
        return (c, int(s) - 1, None)
    m = _region_pat0.match(region)
    if m:
        c, = m.groups()
        return (c, None, None)
    raise ValueError(region)


def sam_intervals(sam, regions=None):
    """
    sam: pysam.Samfile
    """
    ref_lens = dict(zip(sam.references, sam.lengths))
    if regions is None:
        # ':' can be accidentally contained in refenrece names (e.g. HLA:HLA00001)
        for contig in sam.references:
            start = 0
            end = ref_lens[contig]
            yield Interval(contig, start, end, None)
    else:
        for r in regions:
            contig, start, end = parse_region(r)
            if start is None:
                start = 0
            if end is None:
                end = ref_lens[contig]
            yield Interval(contig, start, end, None)


import re
_blank_pat1 = re.compile('^\s+$')
def _is_blank(st):
    return bool(_blank_pat1.match(st))


# This does not work for multimapped records
class LocalInfo:
    def __init__(self, sam, iv, fasta=None, skip_flags=0x004):
        self.contig = iv.contig
        self.start = iv.start
        self.end = iv.end
        if fasta:
            fasta = pysam.Fastafile(fasta)
            ref_bases = fasta.fetch(reference=iv.contig, start=iv.start, end=iv.end).upper()
        else:
            ref_bases = '?' * (iv.end - iv.start)
        self._ref_bases = ref_bases
        self._reads = [Read(rec) for rec in sam.fetch(reference=iv.contig, start=iv.start, end=iv.end) if not rec.flag & skip_flags]

        self._bases_liss = bases_liss = []   # e.g. [['A', 'GC', '-', ...]]
        for read in self._reads:
            bases_lis = read.get_bases_list(iv.start, iv.end)
            if bases_lis is None:  # is it ok?
                continue
            bases_liss.append(bases_lis)

        self._pos_infos = pos_infos = []
        #logging.info(bases_liss[56])
        for bases_lis in zip(*bases_liss):
            length = max(1, max(len(b) for b in bases_lis))
            pos_infos.append({'length': length, 'bases': bases_lis})
        self._max_lens = lens = [r['length'] for r in pos_infos]
        self.align_length = sum(lens)

    def iter_column_info(self):
        """
        Yields: {
            'bases': ['A', 'T', 'A', '-', '-', '-', ...],  # all the bases should have 1 bp width
        }
        """
        #logging.info(self._max_lens)
        #logging.info(self._bases_liss)
        assert len(self._max_lens) == len(self._pos_infos)
        for l, pos_info in zip(self._max_lens, self._pos_infos):
            base_bundles = [fill_text(bs, l, char='-') for bs in pos_info['bases']]
            # TODO distinguish empty and padding
            #logging.info((l, bases_lis))
            for base_column in zip(*base_bundles):
                yield {'bases': list(base_column)}

    def get_left_offset(self, pos):  # pos0
        assert self.start <= pos < self.end
        offset = pos - self.start
        return sum(self._max_lens[:offset], 0)

    def get_ref_seq(self):
        return ''.join(self._ref_bases)

    def get_ref_align(self):
        filled = [fill_text(bs, l, char='-') for l, bs in zip(self._max_lens, self._ref_bases)]   # fill missing bases
        return filled

    def iter_read_aligns(self):  # TODO sorted by position bases
        for read, bases_lis in zip(self._reads, self._bases_liss):
            if bases_lis is None:
                continue
            filled = [fill_text(bs, l, char=(' ' if _is_blank(bs) else '-')) for l, bs in zip(self._max_lens, bases_lis)]   # fill missing bases
            yield (read, filled)

    def iter_align_pairs(self):
        """
        Yields: (Read, (aln, ref_aln))
        """
        ref_bases = self._ref_bases
        for read, bases_lis in zip(self._reads, self._bases_liss):
            if bases_lis is None:
                continue
            b1s = []
            b2s = []
            for b1, b2 in zip(bases_lis, ref_bases):
                fill_len = max(len(b1), len(b2))
                #b1 = fill_text(b1, fill_len, char=(' ' if _is_blank(b1) else '-'))
                #b2 = fill_text(b2, fill_len, char=(' ' if _is_blank(b2) else '-'))
                b1 = fill_text(b1, fill_len, char=(' ' if _is_blank(b1) else '-'))
                b2 = fill_text(b2, fill_len, char=(' ' if _is_blank(b2) else '-'))
                b1s.append(b1)
                b2s.append(b2)
            seq1 = ''.join(b1s)
            seq2 = ''.join(b2s)
            yield read, (seq1, seq2)


class Cigar(csamutil.Cigar):
    """
    >>> Cigar.parse('2I3M3D3M').values == [(Cigar.I, 2), (Cigar.M, 3), (Cigar.D, 3), (Cigar.M, 3)]
    True
    >>> Cigar.parse('3M2I3M2S').values == [(Cigar.M, 3), (Cigar.I, 2), (Cigar.M, 3), (Cigar.S, 2)]
    True

    >>> Cigar([(Cigar.I, 2), (Cigar.M, 3), (Cigar.D, 3), (Cigar.M, 3)]).to_str()
    '2I3M3D3M'
    >>> Cigar([(Cigar.M, 3), (Cigar.I, 2), (Cigar.M, 3), (Cigar.S, 2)]).to_str()
    '3M2I3M2S'

    >>> ' '.join(Cigar([(Cigar.I, 2), (Cigar.M, 3), (Cigar.D, 3), (Cigar.M, 3)]).to_str_tuple())
    '2I 3M 3D 3M'
    >>> ' '.join(Cigar([(Cigar.M, 3), (Cigar.I, 2), (Cigar.M, 3), (Cigar.S, 2)]).to_str_tuple())
    '3M 2I 3M 2S'

    >>> Cigar.parse('1I1I3M2D1D3M').contract().to_str()
    '2I3M3D3M'
    >>> Cigar.parse('3M2I3M2S').contract().to_str()
    '3M2I3M2S'

    Without hard clip length
    >>> Cigar().read_length
    0
    >>> Cigar.parse('4S3M').read_length
    7
    >>> Cigar.parse('10H3D3M').read_length
    13
    >>> Cigar.parse('3M3D2M').read_length
    5
    >>> Cigar.parse('15M2D3I4M').read_length
    22

    With hard clip length
    >>> Cigar().query_length
    0
    >>> Cigar.parse('4S3M').query_length
    7
    >>> Cigar.parse('10H3D3M').query_length
    3
    >>> Cigar.parse('3M3D2M').query_length
    5
    >>> Cigar.parse('15M2D3I4M').query_length
    22

    Reference length
    >>> Cigar().ref_length
    0
    >>> Cigar.parse('4S3M').ref_length
    3
    >>> Cigar.parse('10H3D3M').ref_length
    6
    >>> Cigar.parse('3M3D2M').ref_length
    8
    >>> Cigar.parse('15M2D3I4M').ref_length
    21

    >>> c = Cigar()
    >>> for op, l in reversed(list(Cigar.parse('1I1I3M2D1D3M'))):
    ...     c.prepend((op, l))
    >>> str(c)
    '1I1I3M2D1D3M'

    >>> c = Cigar()
    >>> for op, l in Cigar.parse('1I1I3M2D1D3M'):
    ...     c.append((op, l))
    >>> str(c)
    '1I1I3M2D1D3M'

    >>> c = Cigar()
    >>> for op, l in Cigar.parse('1I1I3M2D1D3M'):
    ...     c.add((op, l))
    >>> str(c)
    '2I3M3D3M'

    >>> Cigar.parse('101M').clips
    (0, 0)
    >>> Cigar.parse('10M1I90M').clips
    (0, 0)
    >>> Cigar.parse('10H100M10S').clips
    (10, 10)
    >>> Cigar.parse('62H39M').clips
    (62, 0)
    >>> Cigar.parse('39M62H').clips
    (0, 62)
    >>> Cigar.parse('20S61M20H').clips
    (20, 20)

    >>> Cigar.parse('10M1I90M').has_clip()
    False
    >>> Cigar.parse('10H100M10S').has_clip()
    True
    >>> Cigar.parse('62H39M').has_clip()
    True

    >>> Cigar.parse('2M4I2M').hard_clip_seq('ATCGATCG')
    'ATCGATCG'
    >>> Cigar.parse('2S4I2S').hard_clip_seq('ATCGATCG')
    'ATCGATCG'
    >>> Cigar.parse('2H4I2H').hard_clip_seq('ATCGATCG')  # only clip when hard clip
    'CGAT'
    >>> Cigar.parse('2H4I2H').hard_clip_seq([0, 1, 2, 3, 4, 5, 6, 7])  # quality scores
    [2, 3, 4, 5]

    # >>> Cigar.parse('2H' '4I' '2H').is_consistent_with(Cigar.parse('2S' '4I' '2S'))
    """


_re_del = re.compile('-+')

def aln2cigar(aln, ref_aln=None):
    """
    >>> assert aln2cigar('NNN') == Cigar.parse('3M')
    >>> assert aln2cigar('NNN') != Cigar.parse('4M')
    >>> assert aln2cigar('---') == Cigar.parse('3D')
    >>> assert aln2cigar('---NNN') == Cigar.parse('3D3M')
    >>> assert aln2cigar('NNN---') == Cigar.parse('3M3D')
    >>> assert aln2cigar('---NNN--') == Cigar.parse('3D3M2D')
    >>> assert aln2cigar('NNN---NN') == Cigar.parse('3M3D2M')

    >>> aln1 = '--------ATATGGGCCATCT'
    >>> assert aln2cigar(aln1) == Cigar.parse('8D13M')
    >>> aln2 = 'ATATATATATACGGG--ATAT'
    >>> assert aln2cigar(aln2) == Cigar.parse('15M2D4M')

    >>> aln2cigar(aln1, aln2).to_str()
    '8D7M2I4M'
    >>> aln2cigar(aln2, aln1).to_str()
    '8I7M2D4M'
    """
    if ref_aln is None:
        return _aln2cigar1(aln)
    return _aln2cigar2(aln, ref_aln)

def _aln2cigar1(aln):
    cigar = Cigar()
    s = 0
    e = len(aln)
    for m in _re_del.finditer(aln):
        s1 = m.start()
        e1 = m.end()
        #print (m, s1, e1)
        if s < s1:
            cigar.append((Cigar.M, s1 - s))
        cigar.append((Cigar.D, e1 - s1))
        s = e1
    if s < e:
        cigar.append((Cigar.M, e - s))
    return cigar

def _aln2cigar2(aln, ref_aln):
    cigar = Cigar()
    op = None
    l = 0
    for a1, a2 in zip(aln, ref_aln):
        if a1 == '-':
            if op == Cigar.D:
                l += 1
                continue
            if op is not None:
                cigar.append((op, l))
            op = Cigar.D
            l = 1
        elif a2 == '-':
            if op == Cigar.I:
                l += 1
                continue
            if op is not None:
                cigar.append((op, l))
            op = Cigar.I
            l = 1
        else:
            if op == Cigar.M:
                l += 1
                continue
            if op is not None:
                cigar.append((op, l))
            op = Cigar.M
            l = 1
    cigar.append((op, l))
    return cigar


class SamInfo(csamutil.SamInfo):
    """

    >>> sam_info = SamInfo(attrdict(references=['chr1'], lengths=[1000], get_tid=lambda x: 1))
    >>> rec_tmpl = dict(tid=1, qname='1', pos=10, aend=110, is_unmapped=0)

    >>> sam_info.get_length('chr1')
    1000
    >>> sam_info.get_length_tid(1)
    1000

    >>> r1 = sam_info.get_read_info(attrdict(rec_tmpl, cigartuples=Cigar.parse('100M').values))
    >>> r1.left_overhang, r1.right_overhang, r1.overhang
    (0, 0, 0)
    >>> r2 = sam_info.get_read_info(attrdict(rec_tmpl, cigartuples=Cigar.parse('5S100M').values))
    >>> r2.left_overhang, r2.right_overhang, r2.overhang
    (5, 0, 5)
    >>> r3 = sam_info.get_read_info(attrdict(rec_tmpl, cigartuples=Cigar.parse('11S100M').values))
    >>> r3.left_overhang, r3.right_overhang, r3.overhang
    (10, 0, 10)
    >>> r4 = sam_info.get_read_info(attrdict(rec_tmpl, cigartuples=Cigar.parse('11S100M12H').values))
    >>> r4.left_overhang, r4.right_overhang, r4.overhang
    (10, 12, 12)
    """


class Read(IntervalMixin):
    has_lclip = None
    has_rclip = None
    lclip = None
    rclip = None
    qlen = None
    qalen = None

    def __init__(self, rec):
        self._rec = rec
        self.name = rec.qname
        self.start = rec.pos
        self.end = rec.pos if rec.aend is None else rec.aend
        self.rec = rec   # pysam.AlignedRead
        self._seq = rec.seq

        self.unmapped = int(rec.is_unmapped)
        self.nins = 0
        self.ndel = 0
        if not self.unmapped:
            seq_len = rec.query_length         # previously rec.rlen
            qstart = rec.query_alignment_start  # previously rec.qstart
            qend = rec.query_alignment_end      # previously rec.qend
            self.qlen = rec.query_alignment_length   # aligned length of query sequence
            self.cigar = Cigar(rec.cigar)
            self.lclip, self.rclip = self.cigar.clips  # soft clip or hard clip
            self.has_lclip = self.lclip > 0
            self.has_rclip = self.rclip > 0
            #self.has_lclip = (qstart > 0)
            #self.has_rclip = (qend < seq_len)

            for op, length in self.cigar:
                if op == Cigar.D:
                    self.ndel += length
                elif op == Cigar.I:
                    self.nins += length

        self.mapq = rec.mapq
        self.alen = 0 if self.end is None else self.end - self.start # aligned length of reference reference sequence
        self.tags = dict(rec.get_tags())
        self.edit = self.tags.get('NM', None)

        #self.mismatches = 0
        self.reverse = int(rec.is_reverse)
        self.suppl = int(rec.is_supplementary)
        self.read1 = int(rec.is_read1)
        self.contig = self.rname = self._rec.reference_name

        # set mate pattern
        # ---------->  ....   # right unmapped
        # ...   <----------  # left unmapped
        self.is_left = int(not rec.is_reverse)  # mate is expected to be in opposite side
        self.mate_miss = int(rec.mate_is_unmapped)
        self.mate_tid = rec.next_reference_id
        self.mate_pos = rec.next_reference_start
        self.tlen = tlen = rec.tlen
        # whether mate is far is determined from tlen and deviation of tlen distribution
        #self.mate_end = rec.pos + tlen if tlen > 0 else (rec.pos is rec.end is None else rec.aend) - tlen # is ok?
        self.mate_invert = 0
        self.mate_back = 0
        self.mate_jump = 0

        if self.mate_miss:
            return
        self.mate_jump = int(rec.tid != self.mate_tid)
        self.mate_invert = int(rec.mate_is_reverse == rec.is_reverse)
        if not rec.is_reverse:
            self.mate_back = int(rec.pnext < rec.pos)
        else:
            self.mate_back = int(rec.pos < rec.pnext)

    @property
    def mate_rname(self):
        return rec.next_reference_name

    @property
    def both_clipped(self):
        return self.has_lclip and self.has_rclip

    @cached_property
    def _bases(self):
        return _get_aln_columns(self._seq, self.rec.cigar)

    @cached_property
    def _quals(self):
        return _get_aln_columns(tuple(self.rec.query_qualities), self.rec.cigar, is_qual=True)

    @cached_property
    def fragmented(self):
        return self.suppl or 'SA' in self.tags

    @cached_property
    def num_frags(self):
        sa = self.tags.get('SA')
        if sa is None:
            return 1
        else:
            return 1 + sa.count(';')

    @cached_property
    def edit_ratio(self):
        return 1. * self.edit / self.alen

    def _get_sublist(self, lis, start, end, missing=' '):
        start_offset = start - self.start
        end_offset = end - self.end
        if not lis:
            return None
        if start_offset > 0:
            lis = lis[start_offset:]
        else:
            lis = [missing] * (-start_offset) + lis
        if end_offset < 0:
            lis = lis[:end_offset]
        else:
            lis = lis + [missing] * end_offset
        return lis

    def get_bases_list(self, start, end, missing=' '):
        return self._get_sublist(self._bases, start, end, missing=missing)

    def get_quals_list(self, start, end, missing=None):
        return self._get_sublist(self._quals, start, end, missing=missing)


def _get_aln_columns(seq, cigar, is_qual=False, del_char=None, n_char=None):
    '''
    Returns:
        [(A|T|C|G|-|*)+] at each columns

    >>> _get_aln_columns('AATCAGTA', Cigar.parse('2I3M3D3M').values)
    ['T', 'C', 'A', '-', '-', '-', 'G', 'T', 'A']
    >>> _get_aln_columns('AATCAGTACC', Cigar.parse('3M2I3M2S').values)
    ['A', 'A', 'TCA', 'G', 'T', 'A']
    >>> _get_aln_columns([1, 2, 3, 4, 5, 6, 7, 8, 9, 10], Cigar.parse('3M2I3M2S').values, is_qual=True)
    [(1,), (2,), (3, 4, 5), (6,), (7,), (8,)]
    '''
    if is_qual:
        del_char = del_char or (20,)
        n_char = n_char or (20,)
    else:
        del_char = del_char or '-'
        n_char = n_char or '*'
    if seq is None:
        return None
    # mask low qualities
    #seq = ''.join(s if q >= 10 else '.' for (s, q) in zip(seq, rec.qual))
    offset = 0
    bases = []
    for (op, length) in cigar:   # TODO (show clip?)
        if op in (0, 7, 8):  # M, X, =
            if is_qual:
                bases.extend([(x,) for x in seq[offset:offset + length]])
            else:
                bases.extend(list(seq[offset:offset + length]))
            offset += length
        elif op == 1:   # I
            if bases:   # inserted seq. at start position will be ignored
                if is_qual:
                    bases[-1] = bases[-1] + tuple(seq[offset:offset + length])  # add to prev sequence
                else:
                    bases[-1] = bases[-1] + seq[offset:offset + length]  # add to prev sequence
            offset += length
        elif op == 2:   # D
            bases.extend([del_char] * length)
        elif op == 3:   # N
            bases.extend([n_char] * length)  # TODO this is rare case, but need revision
        elif op == 4:   # S
            offset += length
        elif op == 5:   # H
            pass
        elif op == 6:   # P
            raise NotImplementedError
    return bases


class ReadCounterBase(object):
    __slots__ = ('start', 'end', 'window',)

    def __init__(self, start, window=1):
        self.start = start
        self.end = start + window

    def add(self, rec):
        """
        Add SAM Record
        """
        raise NotImplementedError


class ReadCounter(ReadCounterBase):
    def __init__(self, start, window):
        super(ReadCounter, self).__init__(start, window=window)
        self.count = 0       # qstart is contained
        self.clip = 0
        self.lclip = 0
        self.rclip = 0
        self.bclip = 0
        self.reverse = 0
        self.mapq1 = 0
        self.mate = MateInfoCounter()
        self.covlen = 0   # covered length within the bin
        self.covlen_mapq1 = 0

    def add(self, rec):
        self.count += 1
        is_lclip = (rec.qstart > 0)
        is_rclip = (rec.qend < rec.rlen)
        is_mapq1 = (rec.mapq <= 1)
        covlen = min(self.end, rec.aend) - max(self.start, rec.pos)
        self.clip += int(is_lclip or is_rclip)
        self.lclip += int(is_lclip)
        self.rclip += int(is_rclip)
        self.bclip += int(is_lclip and is_rclip)
        self.reverse += rec.is_reverse
        self.mapq1 += int(is_mapq1)
        self.mate.add(rec)
        self.covlen += covlen
        if is_mapq1:
            self.covlen_mapq1 += covlen

    def __str__(self):
        return '\t'.join((str(self.start), str(self.end), str(self.__dict__)))


class MateInfoCounter(ReadCounterBase):
    """ Count mate info
    """
    __slots__ = ('_tlen_far', 'unmapped', 'jumped', 'far', 'overlap', 'penetrate', 'rr', 'ff')
    attrs = tuple(filter(lambda x: not x.startswith('_'), __slots__))
    _getter = attrgetter(*attrs)

    def values(self):
        return self._getter(self)

    def items(self):
        return dict(zip(self.attrs, self.values()))

    def __init__(self, tlen_far=700):
        self._tlen_far = tlen_far
        self.unmapped = self.jumped = self.far = self.overlap = self.penetrate = self.rr = self.ff = 0

    def add(self, rec):
        if rec.mate_is_unmapped:
            self.unmapped += 1
            return
        if rec.rnext != rec.tid:
            self.jumped += 1
            return

        # check orientation
        if not rec.is_reverse:
            if not rec.mate_is_reverse:
                self.ff += 1
            elif rec.tlen > self._tlen_far:
                self.far += 1
            elif rec.pnext < rec.aend:
                self.overlap += 1
            elif rec.tlen < 0:
                self.penetrate += 1
        else:
            if rec.mate_is_reverse:
                self.rr += 1
            elif - rec.tlen > self._tlen_far:
                self.far += 1
            elif rec.aend - rec.pnext < - 2 * rec.tlen:  # adhoc
                self.overlap += 1
            elif rec.aend < rec.pnext:
                self.penetrate += 1

    def __str__(self):
        return str(self.items())


class BreakPointCounter(ReadCounterBase):
    def __init__(self, start, window):
        self.start = start
        self.end = start + window
        self.lclips = []

    def add(self, rec):
        self.count += 1
        is_lclip = (rec.qstart > 0)
        is_rclip = (rec.qend < rec.rlen)
        is_clip = is_lclip or is_rclip
        is_mapq1 = (rec.mapq <= 1)
        if is_lclip and self.start <= rec.pos < self.end:
            lclips.append()
        self.lclip += int(is_lclip)
        self.rclip += int(is_rclip)
        self.bclip += int(is_lclip and is_rclip)
        self.reverse += rec.is_reverse
        self.mapq1 += int(is_mapq1)
        self.mate.add(rec)
        self.covlen += covlen
        if is_mapq1:
            self.covlen_mapq1 += covlen

    def __str__(self):
        return '\t'.join((str(self.start), str(self.end), str(self.__dict__)))


class ReadCountGenerator(object):
    def __init__(self, sam, rname, start=0, end=None, window=50, mass='middle', skip_flag=0x904, counter_cls=ReadCounter):
        #self._samit = samit
        if end is None:
            rlens = dict(zip(sam.references, sam.lengths))
            end = rlens[rname]
        self._samit = sam.fetch(reference=rname, start=start, end=end)
        self.counters = []
        self.wstart = start if start % window == 0 else start // window * window # 110 -> 100, 100 -> 100, 149 -> 100  # window: 50
        self.wend = end if end % window == 0 else (end // window + 1) * window   # 110 -> 150, 100 -> 100, 149 -> 150  # window: 50
        self.cstart = self.wstart
        self.cend = self.wstart
        self.window = window
        self.skip_count = 0
        self.skip_flag = skip_flag  # unmapped, secondary or supplementary
        get_read_masses = {
            'start': lambda rec: (rec.pos,),
            'end':   lambda rec: (rec.aend - 1,),   # end position of the alignment (note that aend points one past the last aligned residue)
            'middle': lambda rec: ((rec.pos + rec.aend - 1) / 2.,),   # middle point of the alignment
            'overlap': lambda rec: range(rec.pos, rec.aend, window),   # one point per window overlaped for each alignment
        }
        self.get_read_masses = get_read_masses[mass]
        self._bulk_size = 200
        self._counter_cls = counter_cls

    def _flush(self):
        while self.counters:
            yield self._dequeue()

    def _enqueue(self):
        if self.wend is None:
            enque_size = self._bulk_size
        else:
            enque_size = min(self._bulk_size, (self.wend - self.cend) // self.window)
        self.counters.extend([self._counter_cls(self.cend + self.window * i, self.window) for i in xrange(enque_size)])
        self.cend += self.window * enque_size

    def _dequeue(self):
        self.cstart += self.window
        return self.counters.pop(0)

    def _should_skip(self, rec):
        if rec.flag & self.skip_flag:
            self.skip_count += 1
            return True

    def __iter__(self):
        try:
            while 1:
                rec = next(self._samit)
                if self._should_skip(rec):
                    continue

                start = rec.pos  # 0-based
                end = rec.aend
                while self.cend < start:
                    for counter in self._flush():
                        yield counter

                    self._enqueue()
                    #yield self._dequeue()

                while self.cend < end and self.cend < self.wend:
                    self._enqueue()

                while self.counters and self.counters[0].end < start:
                    yield self._dequeue()

                masses = self.get_read_masses(rec)
                for mass in masses:
                    rec_index = int(mass - self.cstart) // self.window
                    if 0 <= rec_index < len(self.counters):
                        self.counters[rec_index].add(rec)

        except StopIteration:
            pass
        except AssertionError as e:
            logging.error('Invalid record was found: (pos: %s, aend: %s)', rec.pos, rec.aend)
            logging.info(rec)
            raise

        for counter in self._flush():
            yield counter
