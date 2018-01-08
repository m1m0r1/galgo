from __future__ import print_function
import logging
from collections import namedtuple, defaultdict
from itertools import groupby, chain
from operator import attrgetter, itemgetter
from .utils import deprecated


class IntervalMixin(object):
    """
    >>> iv0 = Interval.with_range(100, 200)
    >>> iv1 = Interval('a', 100, 200, None)
    >>> iv2 = Interval('a', 150, 250, None)
    >>> iv3 = Interval('b', 150, 250, None)
    >>> iv4 = Interval('a', 200, 250, None)
    >>> iv5 = Interval('a', 80, 100, None)
    >>> iv6 = Interval('a', 80, 101, None)

    >>> iv0.contains(iv1)
    True
    >>> iv0.intersects(iv2)
    True
    >>> iv0.intersects(iv3)
    True
    >>> iv0.contains(iv3)
    False
    >>> iv0.intersects(iv4)
    False
    >>> iv0.intersects(iv4, adjacent=True)
    True
    >>> iv2.intersects(iv3)
    False
    >>> iv5.intersects(iv0)
    False
    >>> iv5.leftside_of(iv0)
    True
    >>> iv5.rightside_of(iv0)
    False
    >>> iv5.leftside_of(iv0)
    True
    >>> iv6.leftside_of(iv0)
    False
    >>> iv5.rightside_of(iv0)
    False

    >>> iv0.length
    100

    >>> iv0.intersection(iv2)
    Interval(contig=None, start=150, end=200, data=None)
    >>> iv0.intersection(iv4)  # returns nothing

    >>> iv0.intersection(iv4, adjacent=True)
    Interval(contig=None, start=200, end=200, data=None)

    >>> iv0.intersection(iv6)
    Interval(contig=None, start=100, end=101, data=None)
    >>> iv1.move(start=-150, end=150)
    Interval(contig='a', start=-50, end=350, data=None)
    >>> iv1.move(start=-150, end=150).trim()
    Interval(contig='a', start=0, end=350, data=None)
    >>> iv1.move(start=-150, end=150).trim(start=50, end=200)
    Interval(contig='a', start=50, end=200, data=None)
    >>> iv1.iv
    Interval(contig='a', start=100, end=200, data=None)
    """

    @property
    def iv(self):
        return Interval(self.contig, self.start, self.end, None)   # get interval

    @property
    def length(self):
        return self.end - self.start

    def comparable(self, iv):
        return self.contig is None \
                or iv.contig is None \
                or self.contig == iv.contig

    def contains(self, iv):
        return self.comparable(iv) and (self.start <= iv.start and iv.end <= self.end)

    def intersects(self, iv, adjacent=False):
        return self.comparable(iv) and (
                (self.start <= iv.end and iv.start <= self.end)
                if adjacent
                else (self.start < iv.end and iv.start < self.end)
               )

    def intersection(self, iv, adjacent=False):
        if self.intersects(iv, adjacent=adjacent):
            return Interval(None, max(self.start, iv.start), min(self.end, iv.end), None)

    def leftside_of(self, iv):
        return self.comparable(iv) and self.end <= iv.start

    def rightside_of(self, iv):
        return iv.leftside_of(self)

_missing = object()
class Interval(namedtuple('Interval', 'contig,start,end,data'), IntervalMixin):
    __slots__ = ()

    @classmethod
    def with_range(cls, start, end, contig=None, data=None):
        return cls(contig, start, end, data)

    def move(self, start=0, end=0, data=_missing):
        return Interval(self.contig, self.start + start, self.end + end,
                self.data if data is _missing else data)

    def trim(self, start=0, end=float('inf'), data=_missing):
        return Interval(self.contig, max(self.start, start), min(self.end, end),
                self.data if data is _missing else data)

    def replace(self, contig=_missing, start=_missing, end=_missing, data=_missing):
        return Interval(
                self.contig if contig is _missing else contig,
                self.start if start is _missing else start,
                self.end if end is _missing else end,
                self.data if data is _missing else data)


def min_vacant_integer(lis):
    """
    >>> min_vacant_integer([])
    0
    >>> min_vacant_integer([1, 2, 3])
    0
    >>> min_vacant_integer([0, 1, 3])
    2
    >>> min_vacant_integer([0, 1, 2, 3])
    4
    """
    iset = set(lis)
    i = 0
    while 1:
        if i in iset:
            i += 1
            continue
        else:
            return i

# this may be generalized from interval_cover
@deprecated
def get_containers(ivs):
    """ O(n)

    >>> ivs = [
    ...     Interval('a', 100, 250, None),
    ...     Interval('a', 150, 250, None),
    ...     Interval('b', 150, 250, None),
    ...     Interval('a', 200, 250, None),
    ...     Interval('a', 80, 100, None),
    ...     Interval('a', 80, 101, None),
    ... ]
    >>> ivs = get_containers(ivs)
    >>> ivs['a'] == Interval('a', 80, 250, None)
    True
    >>> ivs['b'] == Interval('b', 150, 250, None)
    True
    """
    starts = {}
    ends = {}
    for iv in ivs:
        start = starts.get(iv.contig)
        if start is None or iv.start < start:
            starts[iv.contig] = iv.start
        end = ends.get(iv.contig)
        if end is None or iv.end > end:
            ends[iv.contig] = iv.end

    return dict((c, Interval(c, starts[c], ends[c], None)) for c in starts)


# TODO rename to interval_levels
def interval_pileup(ivs):
    """
    >>> ivs = [
    ...     Interval('a', 100, 250, None),
    ...     Interval('a', 150, 250, None),
    ...     Interval('b', 150, 250, None),
    ...     Interval('a', 200, 250, None),
    ...     Interval('a', 80, 100, None),
    ...     Interval('a', 80, 101, None),
    ... ]
    >>> interval_pileup(ivs)
    [0, 1, 0, 2, 0, 1]
    """
    # maybe O(n^2) algorithm

    iv_levels = []
    def add_level(iv):
        isecs = [lv1 for iv1, lv1 in iv_levels if iv1.intersects(iv)]
        lv = min_vacant_integer(isecs)
        iv_levels.append((iv, lv))
        return lv

    return [add_level(iv) for iv in ivs]


def interval_cluster(ivs, adjacent=False):
    """
    >>> ivs = [
    ...     Interval('a', 80, 100, None),
    ...     Interval('a', 80, 101, None),
    ...     Interval('a', 100, 250, None),
    ...     Interval('a', 150, 250, None),
    ...     Interval('a', 250, 251, None),
    ...     Interval('b', 150, 270, None),
    ...     Interval('b', 160, 200, None),
    ...     Interval('b', 260, 300, None),
    ... ]
    >>> list(cl for cl, iv in interval_cluster(ivs))
    [1, 1, 1, 1, 2, 3, 3, 3]
    >>> list(cl for cl, iv in interval_cluster(ivs, adjacent=True))
    [1, 1, 1, 1, 1, 2, 2, 2]
    """
    cluster_id = 0
    ivs = check_sorted(ivs)
    for contig, ivs in groupby(ivs, attrgetter('contig')):
        cover_iv = None
        for iv in ivs:
            if cover_iv is None or not iv.intersects(cover_iv, adjacent=adjacent):
                cluster_id += 1
                cover_iv = iv
            else:
                cover_iv = interval_cover((cover_iv, iv))
            yield (cluster_id, iv)


def interval_cover(ivs, validate=False):
    """
    >>> ivs = [
    ...     Interval('a', 80, 100, None),
    ...     Interval('a', 80, 101, None),
    ...     Interval('a', 100, 250, None),
    ...     Interval('a', 150, 250, None),
    ...     Interval('a', 250, 251, None),
    ... ]
    >>> interval_cover(ivs)
    Interval(contig='a', start=80, end=251, data=None)
    """
    ivs = tuple(ivs)
    if validate:
        assert len(set(iv.contig for iv in ivs)) == 1
    contig = ivs[0].contig
    start = min(iv.start for iv in ivs)
    end = max(iv.end for iv in ivs)
    return Interval(contig, start, end, None)


# TODO interval tree
class IntervalStore(object):
    """
    >>> store = IntervalStore([
    ...     Interval('a', 80, 100, None),
    ...     Interval('a', 80, 101, None),
    ...     Interval('a', 150, 250, None),
    ...     Interval('b', 250, 251, None),
    ... ])
    >>> coverings = store.get_coverings(Interval('a', 90, 200, None))
    >>> coverings[0].replace(data=None)
    Interval(contig='a', start=90, end=101, data=None)
    >>> coverings[0].data
    [Interval(contig='a', start=80, end=100, data=None), Interval(contig='a', start=80, end=101, data=None)]
    >>> coverings[1].replace(data=None)
    Interval(contig='a', start=150, end=200, data=None)
    >>> coverings[1].data
    [Interval(contig='a', start=150, end=250, data=None)]
    """

    def __init__(self, ivs):
        # just a list structure
        # maybe interval_tree in the future
        self._ivs = ivs if isinstance(ivs, list) else list(ivs)

    def iter_intersections(self, iv, adjacent=False):
        for iv1 in self._ivs:  # TODO this will be O(n)
            if iv.intersects(iv1, adjacent=adjacent):
                yield iv1

    def get_coverings(self, target_iv):
        coverings = []
        isecs = list(self.iter_intersections(target_iv))
        if not isecs:
            return coverings
        isecs = sorted(isecs, key=lambda x: (x.contig, x.start))
        iv_cls = interval_cluster(isecs)

        for cluster_id, cluster_ivs in groupby(iv_cls, lambda x: x[0]):
            ivs = list(iv for _, iv in cluster_ivs)
            iv_cover = interval_cover(ivs)
            trimmed = iv_cover.trim(start=target_iv.start, end=target_iv.end, data=ivs)
            coverings.append(trimmed)
        return coverings


# TODO set contig order?
def check_sorted(ivs, check_end=False):
    """ Confirm chrom ordered, sorted in start position, and equally binned.
    """
    contig_done = set()

    for contig, ivs in groupby(ivs, attrgetter('contig')):
        assert contig not in contig_done, 'contig: {0} is appeared again!'.format(contig)
        contig_done.add(contig)
        prev_start = -float('inf')
        prev_end = -float('inf')
        for iv in ivs:
            assert prev_start <= iv.start, 'Not sorted in start position!'
            if check_end:
                assert prev_end <= iv.end, 'Not sorted in end position!'
            prev_start = iv.start
            prev_end = iv.end
            yield iv


class ChromCompare:
    def __init__(self, chrom_list):
        self._order = {}
        for i, name in enumerate(chrom_list):
            self._order[name] = i

    def get_order(self, chrom):
        return self._order[chrom]


def parse_bed(bedfile, offset=0):
    """
    Yields: Interval
    """
    with open(bedfile) as fp:
        for line in fp:
            if line.startswith('#'):
                continue
            row = line.rstrip('\r\n').split('\t')
            yield Interval(row[0], int(row[1]) - offset, int(row[2]), row)


def iter_bin_interval_flanks(intervals, bin_size=50, lflank=0, rflank=None, validate=True):
    """
    >>> ivs = [
    ...     Interval('chr1', 1000, 1050, 1),
    ...     Interval('chr1', 1050, 1100, 2),
    ...     Interval('chr1', 1150, 1200, 4),
    ...     Interval('chr1', 1200, 1250, 5),
    ...     Interval('chr1', 1250, 1300, 6),
    ...     Interval('chr1', 1300, 1350, 7),
    ...     Interval('chr1', 1400, 1450, 9),
    ...     Interval('chr1', 1450, 1500, 10),
    ...     Interval('chr2', 1000, 1050, 1),
    ... ]

    >>> for rec in iter_bin_interval_flanks(ivs, bin_size=50):
    ...     print (rec.iv.contig, rec.iv.start, rec.iv.data, tuple(iv.data for iv in rec.flanks), tuple(iv.data for iv in rec.adds), tuple(iv.data for iv in rec.removes))
    chr1 1000 1 (1,) (1,) ()
    chr1 1050 2 (2,) (2,) (1,)
    chr1 1150 4 (4,) (4,) (2,)
    chr1 1200 5 (5,) (5,) (4,)
    chr1 1250 6 (6,) (6,) (5,)
    chr1 1300 7 (7,) (7,) (6,)
    chr1 1400 9 (9,) (9,) (7,)
    chr1 1450 10 (10,) (10,) (9,)
    chr2 1000 1 (1,) (1,) ()

    >>> for rec in iter_bin_interval_flanks(ivs, bin_size=50, lflank=200, rflank=200):
    ...     print (rec.iv.contig, rec.iv.start, rec.iv.data, tuple(iv.data for iv in rec.flanks), tuple(iv.data for iv in rec.adds), tuple(iv.data for iv in rec.removes))
    chr1 1000 1 (1, 2, 4) (1, 2, 4) ()
    chr1 1050 2 (1, 2, 4, 5) (5,) ()
    chr1 1150 4 (1, 2, 4, 5, 6, 7) (6, 7) ()
    chr1 1200 5 (1, 2, 4, 5, 6, 7) () ()
    chr1 1250 6 (2, 4, 5, 6, 7, 9) (9,) (1,)
    chr1 1300 7 (4, 5, 6, 7, 9, 10) (10,) (2,)
    chr1 1400 9 (5, 6, 7, 9, 10) () (4,)
    chr1 1450 10 (6, 7, 9, 10) () (5,)
    chr2 1000 1 (1,) (1,) ()
    """

    if rflank is None:
        rflank = bin_size
    assert lflank >= 0 and rflank >= 0
    assert lflank % bin_size == 0, 'The size of the left flank should be a multiple of the bin size'
    assert rflank % bin_size == 0, 'The size of the right flank should be a multiple of the bin size'

    Record = namedtuple('Record', 'iv flanks adds removes')

    for (contig, ivs) in groupby(intervals, itemgetter(0)):  # group by contig
        emit_index = 0
        flanks = []
        last_iv = Interval(contig=contig, start=float('inf'), end=float('inf'), data=None)  # dummy record for avoid emission after ivs consumed
        adds = []
        removes = []

        for iv in chain(ivs, [last_iv]):
            while emit_index < len(flanks) and flanks[emit_index].start + rflank < iv.end:
                while flanks and flanks[0].start < flanks[emit_index].start - lflank:
                    removes.append(flanks.pop(0))
                    emit_index -= 1

                logging.debug('%s', (flanks[emit_index].start, emit_index, iv.start, flanks, adds, removes))
                if emit_index >= len(flanks):
                    break
                yield Record(iv=flanks[emit_index], flanks=tuple(flanks), adds=tuple(adds), removes=tuple(removes))
                emit_index += 1
                adds = []
                removes = []

            if iv is not last_iv:
                adds.append(iv)
                flanks.append(iv)


def align_intervals(ivs1, ivs2, chrom_cmp):
    """
    >>> ivs1 = [
    ...     Interval('chr1', 1000, 1050, '1_1'),
    ...     Interval('chr1', 1000, 1050, '1_1b'),
    ...     Interval('chr1', 1050, 1100, '1_2'),
    ...     Interval('chr1', 1150, 1200, '1_4'),
    ...     Interval('chr1', 1150, 1200, '1_4b'),
    ...     Interval('chr1', 1200, 1250, '1_5'),
    ...     Interval('chr1', 1250, 1300, '1_6'),
    ...     Interval('chr1', 1300, 1350, '1_7'),
    ...     Interval('chr1', 1400, 1450, '1_9'),
    ...     Interval('chr1', 1450, 1500, '1_10'),
    ...     Interval('chr2', 1000, 1050, '2_1'),
    ...     Interval('chr2', 1050, 1100, '2_2'),
    ...     Interval('chr2', 1100, 1150, '2_3'),
    ...     Interval('chr3', 1000, 1050, '3_1'),
    ... ]
    >>> ivs2 = [
    ...     Interval('chr1', 1050, 1100, '1_2'),
    ...     Interval('chr1', 1150, 1200, '1_4'),
    ...     Interval('chr1', 1200, 1250, '1_5'),
    ...     Interval('chr1', 1250, 1300, '1_6'),
    ...     Interval('chr1', 1300, 1350, '1_7'),
    ...     Interval('chr1', 1350, 1400, '1_8'),
    ...     Interval('chr1', 1400, 1450, '1_9'),
    ...     Interval('chr1', 1450, 1500, '1_10'),
    ...     Interval('chr1', 1500, 1550, '1_11'),
    ...     Interval('chr1', 1550, 1600, '1_12'),
    ...     Interval('chr3', 1000, 1050, '3_1'),
    ...     Interval('chr3', 1050, 1100, '3_2'),
    ... ]
    >>> l = list(align_intervals(ivs1, ivs2, ChromCompare(['chr1', 'chr2', 'chr3'])))
    >>> for iv in l:
    ...     print (iv.contig, iv.start, iv.end, iv.data[0], iv.data[1])
    chr1 1000 1050 1_1 None
    chr1 1000 1050 1_1b None
    chr1 1050 1100 1_2 1_2
    chr1 1150 1200 1_4 1_4
    chr1 1150 1200 1_4b None
    chr1 1200 1250 1_5 1_5
    chr1 1250 1300 1_6 1_6
    chr1 1300 1350 1_7 1_7
    chr1 1350 1400 None 1_8
    chr1 1400 1450 1_9 1_9
    chr1 1450 1500 1_10 1_10
    chr1 1500 1550 None 1_11
    chr1 1550 1600 None 1_12
    chr2 1000 1050 2_1 None
    chr2 1050 1100 2_2 None
    chr2 1100 1150 2_3 None
    chr3 1000 1050 3_1 3_1
    chr3 1050 1100 None 3_2
    """
    ivs1 = groupby(iter(ivs1), lambda x: x[0])
    ivs2 = groupby(iter(ivs2), lambda x: x[0])
    chrom1, iivs1 = next(ivs1, (None, None))
    chrom2, iivs2 = next(ivs2, (None, None))

    while 1:
        i1 = chrom_cmp.get_order(chrom1) if chrom1 is not None else -1
        i2 = chrom_cmp.get_order(chrom2) if chrom2 is not None else -1

        if i1 < i2:
            for iv1 in iivs1:
                yield Interval(chrom1, iv1.start, iv1.end, data=(iv1.data, None))
            chrom1, iivs1 = next(ivs1, (None, None))
            continue
        elif i2 < i1:
            for iv2 in iivs2:
                yield Interval(chrom2, iv2.start, iv2.end, data=(None, iv2.data))
            chrom2, iivs2 = next(ivs2, (None, None))
            continue
        elif i1 == i2 == -1:
            return

        iv1 = next(iivs1, None)
        iv2 = next(iivs2, None)
        while 1:
            if iv1 is None and iv2 is None:
                break
            elif iv1 is None:
                for iv in chain((iv2, ), iivs2):
                    yield Interval(iv.contig, iv.start, iv.end, data=(None, iv.data))
                break
            elif iv2 is None:
                for iv in chain((iv1, ), iivs1):
                    yield Interval(iv.contig, iv.start, iv.end, data=(iv.data, None))
                break

            o1 = (iv1.start, iv1.end)
            o2 = (iv2.start, iv2.end)
            if o1 == o2:
                yield Interval(chrom1, iv1.start, iv1.end, data=(iv1.data, iv2.data))
                iv1 = next(iivs1, None)
                iv2 = next(iivs2, None)
                continue
            elif o1 < o2:
                yield Interval(chrom1, iv1.start, iv1.end, data=(iv1.data, None))
                iv1 = next(iivs1, None)
                continue
            else:
                yield Interval(chrom2, iv2.start, iv2.end, data=(None, iv2.data))
                iv2 = next(iivs2, None)
                continue
        chrom1, iivs1 = next(ivs1, (None, None))
        chrom2, iivs2 = next(ivs2, (None, None))
