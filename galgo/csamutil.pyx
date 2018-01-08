#cython: profile=True
import re
from itertools import groupby
from .utils import cached_property
from pysam import AlignedSegment

_re_cigar = re.compile('(\d+)([MIDNSHPX=])')
_cigar_ops = 'MIDNSHPX='
_cigar_op_codes = dict((_op, _code) for _code, _op in enumerate(_cigar_ops))
cdef class Cigar(object):
    M = 0
    I = 1
    D = 2
    N = 3
    S = 4
    H = 5
    P = 6
    EQ = 7
    X = 8  # '='

    _read_len_tags  = set((M,    I, N, S, H, EQ, X))  # count for read length
    _query_len_tags = set((M,    I, N, S,    EQ, X))  # count for query length
    _ref_len_tags   = set((M, D,    N,       EQ, X))  # count for read length
    cdef:
        public list values

    def __init__(self, list values=None):
        self.values = [] if values is None else values[:]

    @classmethod
    def parse(cls, cigar_str):
        ret = []
        for m in _re_cigar.finditer(cigar_str):
            l, op = m.groups()
            ret.append((_cigar_op_codes[op], int(l)))
        return Cigar(ret)

    def to_str(self):
        return ''.join(str(l) + _cigar_ops[code] for code, l in self.values)

    def to_str_tuple(self):
        return tuple(str(l) + _cigar_ops[code] for code, l in self.values)

    def contract(self):
        return Cigar([(op, sum(x[1] for x in lis)) for op, lis in groupby(self.values, lambda x: x[0])])

    def __iter__(self):
        return iter(self.values)

    def __richcmp__(self, other, op):
        if op == 2:  # ==
            if isinstance(other, Cigar):
                return self.values == other.values
            else:
                return self.values == other
        elif op == 3:
            if isinstance(other, Cigar):
                return self.values != other.values
            else:
                return self.values != other
        else:
            raise TypeError

    @property
    def read_length(self):
        return sum((l for op, l in self.values if op in self._read_len_tags), 0)

    @property
    def query_length(self):
        return sum((l for op, l in self.values if op in self._query_len_tags), 0)

    @property
    def ref_length(self):
        return sum((l for op, l in self.values if op in self._ref_len_tags), 0)

    def prepend(self, tuple pair):
        self.values.insert(0, pair)

    def append(self, tuple pair):
        self.values.append(pair)

    def add(self, tuple pair):
        if self.values and self.values[-1][0] == pair[0]:
            self.values[-1] = (pair[0], self.values[-1][1] + pair[1])
        else:
            self.values.append(pair)

    def __str__(self):
        return self.to_str()

    @property
    def clips(self):
        return self.lclip, self.rclip

    def has_clip(self):
        return self.lclip > 0 or self.rclip > 0

    @property
    def lclip(self):
        if self.values[0][0] in (self.S, self.H):
            return self.values[0][1]
        return 0

    @property
    def rclip(self):
        if self.values[-1][0] in (self.S, self.H):
            return self.values[-1][1]
        return 0

    def hard_clip_seq(self, seq, hard_only=True):
        if not self.values:
            return seq
        l = self.values[0][1] if self.values[0][0] == self.H else 0
        r = self.values[-1][1] if self.values[-1][0] == self.H else 0
        return seq[l:len(seq) - r]


NINT = -1
cdef class SamReadInfo(object):
    cdef:
        SamInfo _info
        public rec
        Cigar _cigar
        int _lclip, _rclip, _left_overhang, _right_overhang

    def __cinit__(self):
        self._lclip = NINT
        self._rclip = NINT
        self._left_overhang = NINT
        self._right_overhang = NINT

    def __init__(self, rec, _saminfo):
        self._info = _saminfo
        self.rec = rec

    @property
    def edit(self):
        tag = dict(self.rec.get_tags())
        edit = tag.get('NM', 0)
        return edit

    @property
    def cigar(self):
        if self._cigar is None:
            self._cigar = Cigar(self.rec.cigartuples)
        return self._cigar

    @property
    def lclip(self):
        if self._lclip == NINT:
            self._lclip = self.cigar.lclip
        return self._lclip

    @property
    def rclip(self):
        if self._rclip == NINT:
            self._rclip = self.cigar.rclip
        return self._rclip

    @property
    def overhang(self):   # TODO how to consider region filled with N bases ?
        cdef int lo = self.left_overhang
        cdef int ro = self.right_overhang
        return max(lo, ro)

    @property
    def left_overhang(self):
        cdef int lclip
        if self._left_overhang == NINT:
            lclip = self.cigar.lclip
            if self.rec.is_unmapped:
                self._left_overhang = 0
            elif lclip > 0:
                self._left_overhang = min(self.rec.pos, lclip)
            else:
                self._left_overhang = 0
        return self._left_overhang

    @property
    def right_overhang(self):
        cdef int rclip, rdist
        if self._right_overhang == NINT:
            if self.rec.is_unmapped:
                self._right_overhang = 0
            else:
                rclip = self.cigar.rclip
                if rclip:
                    rdist = self._info.get_length_tid(self.rec.tid) - self.rec.aend
                    self._right_overhang = min(rclip, rdist)
                else:
                    self._right_overhang = 0
        return self._right_overhang


cdef class SamInfo(object):
    cdef:
        _sam
        list _filters
        dict _lengths

    def __init__(self, sam):
        """
        sam: pysam.Samfile object
        """
        self._sam = sam
        self._filters = []
        self._lengths = {sam.get_tid(name): length for name, length in zip(sam.references, sam.lengths)}

    cpdef int get_length_tid(self, int tid):
        return self._lengths[tid]

    def get_length(self, refname):
        return self._lengths[self._sam.get_tid(refname)]

    def get_read_info(self, rec):
        return SamReadInfo(rec, _saminfo=self)

