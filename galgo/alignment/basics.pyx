from __future__ import print_function
import logging
import numpy as np
from ..matrix import get_banded_matrix
cimport numpy as np
from ..sselogsumexp import logsumexp

cdef class LinearGapScore(object):
    cdef public:
        int match, mismatch, gap1
        int gap2

    def __init__(self, match=1, mismatch=-2, gap1=-2, gap2=-2):
        self.match = match
        self.mismatch = mismatch
        self.gap1 = gap1
        self.gap2 = gap2

cdef class AffineGapScore(object):
    cdef public:
        int match, mismatch, gap1_open, gap1_ext, gap2_open, gap2_ext

    def __init__(self, match=1, mismatch=-2, gap1_open=-2, gap1_ext=-1, gap2_open=-2, gap2_ext=-1):
        self.match = match
        self.mismatch = mismatch
        self.gap1_open = gap1_open
        self.gap1_ext = gap1_ext
        self.gap2_open = gap2_open
        self.gap2_ext = gap2_ext


class NWAlignmentTrace:
    """
    >>> score = LinearGapScore(1, -2, -2, -2)
    >>> seq1 = 'ATATGCTATAT'
    >>> seq2 = 'ATATCTGATCT'
    >>> a = NWAlignmentTrace(seq1, seq2, move_score=score)
    >>> r = a.get_best()
    >>> r['aln1']
    'ATATGCT-ATAT'
    >>> r['aln2']
    'ATAT-CTGATCT'
    >>> r['score'] == 4 - 2 + 2 - 2 + 2 - 2 + 1
    True
    """

    move_score = LinearGapScore()

    def __init__(self, seq1, seq2, move_score=None, dtype=int):
        self.seq1 = seq1
        self.seq2 = seq2
        self.len1 = len(seq1)
        self.len2 = len(seq2)
        self.dtype = dtype
        if move_score:
            self.move_score = move_score
        self._fill_matrix()

    class moves:
        MATCH = 0
        GAP1 = 1
        GAP2 = 2

    def _fill_matrix(self):
        cdef int len1, len2, i1, i2
        cdef np.ndarray M
        len1 = self.len1
        len2 = self.len2
        self._M = M = np.zeros((len1 + 1, len2 + 1), dtype=self.dtype)
        self._trace = trace = np.zeros((len1 + 1, len2 + 1), dtype=int)
        sc = self.move_score
        moves = self.moves
        ops = (moves.MATCH, moves.GAP1, moves.GAP2)

        for i1 in xrange(1, len1 + 1):   # fill left
            M[i1, 0] = M[i1 - 1, 0] + sc.gap2
            trace[i1, 0] = moves.GAP2
        for i2 in xrange(1, len2 + 1):   # fill top
            M[0, i2] = M[0, i2 - 1] + sc.gap1
            trace[0, i2] = moves.GAP1

        for i1 in xrange(1, len1 + 1):
            r1 = self.seq1[i1-1]
            for i2 in xrange(1, len2 + 1):
                r2 = self.seq2[i2-1]
                m = sc.match if r1 == r2 else sc.mismatch
                scs = (
                    M[i1-1, i2-1] + m,
                    M[i1, i2-1] + sc.gap1,
                    M[i1-1, i2] + sc.gap2,
                )
                M[i1, i2] = max(scs)
                trace[i1, i2] = ops[np.argmax(scs)]  # TODO inline

    def get_best(self):
        """ Returns aln1, aln2, score
        """
        score = self._M[self.len1, self.len2]
        moves = self.moves
        raln1 = []
        raln2 = []
        i1 = self.len1
        i2 = self.len2
        seq1, seq2 = self.seq1, self.seq2
        while i1 > 0 or i2 > 0:
            op = self._trace[i1, i2]
            if op == moves.MATCH:
                a1 = seq1[i1-1]
                a2 = seq2[i2-1]
                i1 -= 1
                i2 -= 1
            elif op == moves.GAP1:
                a1 = '-'
                a2 = seq2[i2-1]
                i2 -= 1
            elif op == moves.GAP2:
                a1 = seq1[i1-1]
                a2 = '-'
                i1 -= 1
            else:
                raise Exception('undefined move')
            raln1.append(a1)
            raln2.append(a2)

        aln1 = ''.join(reversed(raln1))
        aln2 = ''.join(reversed(raln2))
        return {'score': score, 'aln1': aln1, 'aln2': aln2}


class NWAlignment(object):
    """
    >>> score = LinearGapScore(1, -2, -2, -2)
    >>> seq1 = 'ATATGCTATAT'
    >>> seq2 = 'ATATCTGATCT'
    >>> a = NWAlignment(seq1, seq2, move_score=score)
    >>> r = a.get_best()
    >>> r['aln1']
    'ATATGCT-ATAT'
    >>> r['aln2']
    'ATAT-CTGATCT'
    >>> r['score'] == 4 - 2 + 2 - 2 + 2 - 2 + 1
    True
    """

    move_score = LinearGapScore()

    def __init__(self, seq1, seq2, move_score=None, dtype=int):
        self.seq1 = seq1
        self.seq2 = seq2
        self.len1 = len(seq1)
        self.len2 = len(seq2)
        self.dtype = dtype
        if move_score:
            self.move_score = move_score
        self._fill_matrix()

    def _fill_matrix(self):
        len1 = self.len1
        len2 = self.len2
        self._M = M = np.zeros((len1 + 1, len2 + 1), dtype=self.dtype)
        sc = self.move_score

        for i1 in xrange(1, len1 + 1):   # fill left
            M[i1, 0] = M[i1 - 1, 0] + sc.gap2
        for i2 in xrange(1, len2 + 1):   # fill top
            M[0, i2] = M[0, i2 - 1] + sc.gap1

        for i1 in xrange(1, len1 + 1):
            r1 = self.seq1[i1-1]
            for i2 in xrange(1, len2 + 1):
                r2 = self.seq2[i2-1]
                m = sc.match if r1 == r2 else sc.mismatch
                M[i1, i2] = max(
                    M[i1-1, i2-1] + m,
                    M[i1, i2-1] + sc.gap1,
                    M[i1-1, i2] + sc.gap2,
                )

    def get_best(self):
        """ Returns aln1, aln2, score
        """
        score = self._M[self.len1, self.len2]
        sc = self.move_score
        seq1, seq2 = self.seq1, self.seq2
        M = self._M
        raln1 = []
        raln2 = []
        i1 = self.len1
        i2 = self.len2
        while i1 > 0 or i2 > 0:
            m = M[i1, i2]
            matchable = i1 > 0 and i2 > 0
            if matchable:
                matched = self.seq1[i1-1] == self.seq2[i2-1]
                m_score = sc.match if matched else sc.mismatch
            if matchable and m == M[i1-1, i2-1] + m_score:  # match
                a1 = seq1[i1-1]
                a2 = seq2[i2-1]
                i1 -= 1
                i2 -= 1
            elif i2 > 0 and m == M[i1, i2-1] + sc.gap1:  # gap1
                a1 = '-'
                a2 = seq2[i2-1]
                i2 -= 1
            elif i1 > 0 and m == M[i1-1, i2] + sc.gap2:  # gap2
                a1 = seq1[i1-1]
                a2 = '-'
                i1 -= 1
            else:
                raise Exception('undefined move')
            raln1.append(a1)
            raln2.append(a2)

        aln1 = ''.join(reversed(raln1))
        aln2 = ''.join(reversed(raln2))
        return {'score': score, 'aln1': aln1, 'aln2': aln2}


class SWAlignment:
    """
    >>> score = LinearGapScore(1, -2, -2, -2)
    >>> aln1 = 'ATATGCTGATAT'
    >>> aln2 = 'ATAT-CTGATCT'
    >>> seq1 = aln1.replace('-', '')
    >>> seq2 = aln2.replace('-', '')
    >>> a = SWAlignment(seq1, seq2, move_score=score)
    >>> r = a.get_best()
    >>> r['start1'], r['aln1']
    (1, 'ATATGCTGAT')
    >>> r['start2'], r['aln2']
    (1, 'ATAT-CTGAT')
    >>> r['score'] == 4 - 2 + 5 == 7
    True

    >>> aln1 = 'AAGATATGCTAATAT'
    >>> aln2 = '----TAT-CTAA---'
    >>> seq1 = aln1.replace('-', '')
    >>> seq2 = aln2.replace('-', '')
    >>> a = SWAlignment(seq1, seq2, move_score=score)
    >>> r = a.get_best()
    >>> r['start1'], r['aln1']
    (5, 'TATGCTAA')
    >>> r['start2'], r['aln2']
    (1, 'TAT-CTAA')
    >>> r['score'] == 3 - 2 + 4 == 5
    True
    """

    move_score = LinearGapScore()

    def __init__(self, seq1, seq2, move_score=None, dtype=int):
        self.seq1 = seq1
        self.seq2 = seq2
        self.len1 = len(seq1)
        self.len2 = len(seq2)
        self.dtype = dtype
        if move_score:
            self.move_score = move_score
        self._fill_matrix()

    def _fill_matrix(self):
        len1 = self.len1
        len2 = self.len2
        self._M = M = np.zeros((len1 + 1, len2 + 1), dtype=self.dtype)
        sc = self.move_score

        for i1 in xrange(1, len1 + 1):   # fill left
            M[i1, 0] = max(M[i1 - 1, 0] + sc.gap2, 0)
        for i2 in xrange(1, len2 + 1):   # fill top
            M[0, i2] = max(M[0, i2 - 1] + sc.gap1, 0)

        for i1 in xrange(1, len1 + 1):
            r1 = self.seq1[i1-1]
            for i2 in xrange(1, len2 + 1):
                r2 = self.seq2[i2-1]
                m = sc.match if (r1 == r2) else sc.mismatch
                M[i1, i2] = max(
                    M[i1-1, i2-1] + m,
                    M[i1, i2-1] + sc.gap1,
                    M[i1-1, i2] + sc.gap2,
                    0,
                )
        logging.debug(M)

    def get_best(self):
        """ Returns aln1, aln2, score, start1, start2
        """
        sc = self.move_score
        M = self._M
        raln1 = []
        raln2 = []
        i1, i2 = np.unravel_index(np.argmax(M), M.shape)  # get best end position
        last_score = score = M[i1, i2]
        while last_score > 0:
            matchable = i1 > 0 and i2 > 0
            if matchable:
                matched = self.seq1[i1-1] == self.seq2[i2-1]
                m_score = sc.match if matched else sc.mismatch
            if matchable and last_score == M[i1-1, i2-1] + m_score:  # match
                a1 = self.seq1[i1-1]
                a2 = self.seq2[i2-1]
                i1 -= 1
                i2 -= 1
            elif i2 > 0 and last_score == M[i1, i2-1] + sc.gap1:  # gap1
                a1 = '-'
                a2 = self.seq2[i2-1]
                i2 -= 1
            elif i1 > 0 and last_score == M[i1-1, i2] + sc.gap2:  # gap2
                a1 = self.seq1[i1-1]
                a2 = '-'
                i1 -= 1
            else:
                raise Exception('undefined move')
            raln1.append(a1)
            raln2.append(a2)
            last_score = M[i1, i2]

        start1 = i1 + 1
        start2 = i2 + 1
        aln1 = ''.join(reversed(raln1))
        aln2 = ''.join(reversed(raln2))
        return {'score': score, 'aln1': aln1, 'aln2': aln2, 'start1': start1, 'start2': start2}


class NWAlignmentBandedNaive:
    """
    >>> score = LinearGapScore(1, -2, -2, -2)
    >>> seq1 = 'ATATGCTATAT'
    >>> seq2 = 'ATATCTGATCT'
    >>> a = NWAlignmentBandedNaive(seq1, seq2, band1=3, band2=3, move_score=score)
    >>> r = a.get_best()
    >>> r['aln1']
    'ATATGCT-ATAT'
    >>> r['aln2']
    'ATAT-CTGATCT'
    >>> r['score'] == 4 - 2 + 2 - 2 + 2 - 2 + 1
    True
    """

    move_score = LinearGapScore()

    def __init__(self, seq1, seq2, move_score=None, band1=None, band2=None, dtype=int):
        """
        Requirements:
            b1 = band1 in (0, l1)
            b2 = band2 in (0, l2)
            b1' = band1 - (len1 - len2) in (0, l1)
            b2' = band2 - (len2 - len1) in (0, l2)

             b1
             .--.----
          b2 |+  .
             . +  .
             |. +  .
             | . +  .

                     .  -  .|
                      .  -  .
                       .  - | b'2
                     -------.
                         b'1
        """
        self.seq1 = seq1
        self.seq2 = seq2
        self.len1 = l1 = len(seq1)
        self.len2 = l2 = len(seq2)
        w = int(round(np.sqrt(max(l1, l2))))
        self.band1 = min(max(l1 - l2, w), l1) if band1 is None else band1
        self.band2 = min(max(l2 - l1, w), l2) if band2 is None else band2
        assert max(l1 - l2, 1) <= self.band1 <= l1
        assert max(l2 - l1, 1) <= self.band2 <= l2
        self.dtype = dtype
        if move_score:
            self.move_score = move_score
        self._fill_matrix()

    def _fill_matrix(self):
        l1 = self.len1
        l2 = self.len2
        b1 = self.band1
        b2 = self.band2
        w2 = b1 + b2
        self._B = B = np.zeros((l1 + 1, w2 + 1), dtype=self.dtype)
        sc = self.move_score
        MIN = -float('inf')
        for i1 in xrange(1, b1 + 1):
            j1 = i1
            j2 = b1 - i1
            B[j1, j2] = B[j1-1, j2+1] + sc.gap2
        for i2 in xrange(1, b2 + 1):
            j2 = b1 + i2
            B[0, j2] = B[0, j2-1] + sc.gap1

        for j1 in xrange(1, l1 + 1):
            i1 = j1
            r1 = self.seq1[i1-1]
            j2_from = max(0, b1 - j1 + 1)
            j2_to = w2 + 1 + min(0, l2 - b2 - j1)
            for j2 in xrange(j2_from, j2_to):
                i2 = j2 + i1 - b1
                r2 = self.seq2[i2-1]
                m = sc.match if r1 == r2 else sc.mismatch
                B[j1, j2] = max(
                    B[j1-1, j2] + m,
                    (B[j1, j2-1] + sc.gap1), #if j2 > max(0, b1 - j1) else MIN,
                    (B[j1-1, j2+1] + sc.gap2) if j2 < w2 else MIN,
                )

    def get_best(self):
        """ Returns aln1, aln2, score
        """
        l1 = self.len1
        l2 = self.len2
        b1 = self.band1
        b2 = self.band2
        w2 = b1 + b2
        j1 = l1
        j2 = b1 - (l1 - l2)
        i1 = j1
        i2 = j1 + j2 - b1
        sc = self.move_score
        B = self._B
        score = B[j1, j2]
        raln1 = []
        raln2 = []
        seq1, seq2 = self.seq1, self.seq2
        while i1 > 0 or i2 > 0:
            #logging.debug((i1, i2, j1, j2))
            m = B[j1, j2]
            matchable = i1 > 0 and i2 > 0
            if matchable:
                matched = self.seq1[i1-1] == self.seq2[i2-1]
                m_score = sc.match if matched else sc.mismatch
            if matchable and m == B[j1-1, j2] + m_score:  # match
                a1 = seq1[i1-1]
                a2 = seq2[i2-1]
                i1 -= 1
                i2 -= 1
                j1 = i1
                j2 = i2 - i1 + b1
            elif j2 > max(0, b1 - j1) and m == B[j1, j2-1] + sc.gap1:  # gap1
                a1 = '-'
                a2 = seq2[i2-1]
                i2 -= 1
                j2 = i2 - i1 + b1
            elif m == B[j1-1, j2+1] + sc.gap2:  # gap2
                a1 = seq1[i1-1]
                a2 = '-'
                i1 -= 1
                j1 = i1
                j2 = i2 - i1 + b1
            else:
                raise Exception('undefined move')
            raln1.append(a1)
            raln2.append(a2)

        aln1 = ''.join(reversed(raln1))
        aln2 = ''.join(reversed(raln2))
        return {'score': score, 'aln1': aln1, 'aln2': aln2}


class NWAlignmentBandedSparse(NWAlignmentBandedNaive):
    """
    >>> score = LinearGapScore(1, -2, -2, -2)
    >>> seq1 = 'ATATGCTATAT'
    >>> seq2 = 'ATATCTGATCT'
    >>> a = NWAlignmentBandedSparse(seq1, seq2, band1=3, band2=3, move_score=score)
    >>> r = a.get_best()
    >>> r['aln1']
    'ATATGCT-ATAT'
    >>> r['aln2']
    'ATAT-CTGATCT'
    >>> r['score'] == 4 - 2 + 2 - 2 + 2 - 2 + 1
    True
    """

    def _fill_matrix(self):
        l1 = self.len1
        l2 = self.len2
        b1 = self.band1
        b2 = self.band2
        width = b1 + b2
        from scipy.sparse import spdiags
        #data = np.zeros((width + 1, l1 + 1), dtype=self.dtype)
        #diags = np.arange(-b1, b2 + 1)
        #self._M = M = spdiags(data, diags, l1+1, l2+1).tolil()
        #self._M = M = spdiags(data, diags, l1+1, l2+1).todense()
        #self._M = M = spdiags(data, diags, l1+1, l2+1).todok()
        self._M = M = get_banded_matrix(b1, b2, l1+1, l2+1, dtype=self.dtype)
        sc = self.move_score
        MIN = -float('inf')
        for i1 in xrange(1, b1 + 1):
            M[i1, 0] = M[i1-1, 0] + sc.gap2
        for i2 in xrange(1, b2 + 1):
            M[0, i2] = M[0, i2-1] + sc.gap1

        for i1 in xrange(1, l1 + 1):
            r1 = self.seq1[i1-1]
            i2_from = max(0, i1 - b1)
            i2_to = min(i1 + b2, l2)
            for i2 in xrange(i2_from, i2_to + 1):
                r2 = self.seq2[i2-1]
                m = sc.match if r1 == r2 else sc.mismatch
                M[i1, i2] = max(
                    M[i1-1, i2-1] + m if i1 > max(0, i2 - b2) and i2 > max(0, i1 - b1) else MIN,
                    (M[i1, i2-1] + sc.gap1) if i2 > max(0, i1 - b1) else MIN,
                    (M[i1-1, i2] + sc.gap2) if i1 > max(0, i2 - b2) else MIN,
                )

    def get_best(self):
        """ Returns aln1, aln2, score
        """
        l1 = self.len1
        l2 = self.len2
        b1 = self.band1
        b2 = self.band2
        i1 = l1
        i2 = l2
        sc = self.move_score
        M = self._M
        score = M[i1, i2]
        raln1 = []
        raln2 = []
        seq1, seq2 = self.seq1, self.seq2
        while i1 > 0 or i2 > 0:
            assert max(0, i1 - b1) <= i2 <= min(l2, i1 + b2) + 1, 'i1: {} i2: {}'.format(i1, i2)
            assert max(0, i2 - b2) <= i1 <= min(l1, i2 + b1) + 1, 'i1: {} i2: {}'.format(i1, i2)
            m = M[i1, i2]
            matchable = i1 > 0 and i2 > 0
            if matchable:
                matched = self.seq1[i1-1] == self.seq2[i2-1]
                m_score = sc.match if matched else sc.mismatch
            if matchable and m == M[i1-1, i2-1] + m_score:  # match
                a1 = seq1[i1-1]
                a2 = seq2[i2-1]
                i1 -= 1
                i2 -= 1
            elif i2 > max(0, i1 - b1) and m == M[i1, i2-1] + sc.gap1:  # gap1
                a1 = '-'
                a2 = seq2[i2-1]
                i2 -= 1
            elif i1 > max(0, i2 - b2) and m == M[i1-1, i2] + sc.gap2:  # gap2
                a1 = seq1[i1-1]
                a2 = '-'
                i1 -= 1
            else:
                raise Exception('undefined move')
            raln1.append(a1)
            raln2.append(a2)

        aln1 = ''.join(reversed(raln1))
        aln2 = ''.join(reversed(raln2))
        return {'score': score, 'aln1': aln1, 'aln2': aln2}

NWAlignmentBanded = NWAlignmentBandedSparse
#NWAlignmentBanded = NWAlignmentBandedNaive

class NWGAlignment:
    """
    >>> score = AffineGapScore(match=1, mismatch=-2, gap1_open=-2, gap1_ext=-1, gap2_open=-2, gap2_ext=-1)
    >>> aln1 = 'AT--G'
    >>> aln2 = 'ATAAG'
    >>> seq1 = aln1.replace('-', '')
    >>> seq2 = aln2.replace('-', '')
    >>> a = NWGAlignment(seq1, seq2, move_score=score)
    >>> r = a.get_best()
    >>> (r['aln1'], r['aln2']) == (aln1, aln2)
    True
    >>> r['score'] == 2 - 2 - 1 + 1 == 0
    True

    >>> aln1 = 'ATATATATATATGGG--ATAT'
    >>> aln2 = '--------ATATGGGCCATCT'
    >>> seq1 = aln1.replace('-', '')
    >>> seq2 = aln2.replace('-', '')
    >>> a = NWGAlignment(seq1, seq2, move_score=score)
    >>> r = a.get_best()
    >>> (r['aln1'], r['aln2']) == (aln1, aln2)
    True
    >>> r['score'] == -8 + 7 - 2 - 1 + 2 - 2 + 1 == -3
    True

    >>> aln1 = 'ATATGCTGATAT-----'
    >>> aln2 = '-----CTGATCTNNNNG'
    >>> seq1 = aln1.replace('-', '')
    >>> seq2 = aln2.replace('-', '')
    >>> a = NWGAlignment(seq1, seq2, move_score=score)
    >>> r = a.get_best()
    >>> (r['aln1'], r['aln2']) == (aln1, aln2)
    True
    >>> r['score'] == -5 + 5 - 2 + 1 - 2 - 4 == -7
    True
    """

    move_score = AffineGapScore()

    def __init__(self, seq1, seq2, move_score=None, dtype=int):
        self.seq1 = seq1
        self.seq2 = seq2
        self.len1 = len(seq1)
        self.len2 = len(seq2)
        self.dtype = dtype
        self._MIN = -60000 - self.move_score.mismatch * (self.len1 + self.len2)
        if move_score:
            self.move_score = move_score
        self._fill_matrix()

    def _fill_matrix(self):
        len1 = self.len1
        len2 = self.len2
        self._M = M = np.zeros((len1 + 1, len2 + 1), dtype=self.dtype)
        self._G1 = G1 = np.zeros((len1 + 1, len2 + 1), dtype=self.dtype)
        self._G2 = G2 = np.zeros((len1 + 1, len2 + 1), dtype=self.dtype)
        sc = self.move_score

        MIN = self._MIN
        M[0, 1:] = MIN
        M[1:, 0] = MIN
        G1[1:, 0] = MIN
        G2[0, 1:] = MIN

        for i1 in xrange(1, len1 + 1):
            G2[i1, 0] = G2[i1-1, 0] + sc.gap2_ext
        for i2 in xrange(1, len2 + 1):
            G1[0, i2] = G1[0, i2-1] + sc.gap1_ext

        for i1 in xrange(1, len1 + 1):
            r1 = self.seq1[i1-1]
            for i2 in xrange(1, len2 + 1):
                r2 = self.seq2[i2-1]
                m = sc.match if r1 == r2 else sc.mismatch
                M[i1, i2] = max(
                    M[i1-1, i2-1],
                    G1[i1-1, i2-1],
                    G2[i1-1, i2-1],
                ) + m
                G1[i1, i2] = max(
                    M[i1, i2-1] + sc.gap1_open,
                    G1[i1, i2-1] + sc.gap1_ext,
                )
                G2[i1, i2] = max(
                    M[i1, i2-1] + sc.gap2_open,
                    G2[i1, i2-1] + sc.gap2_ext,
                )

    def get_best(self):
        """ Returns aln1, aln2, score
        """
        sc = self.move_score
        M = self._M
        raln1 = []
        raln2 = []
        seq1 = self.seq1
        seq2 = self.seq2
        i1 = self.len1
        i2 = self.len2
        move_M, move_G1, move_G2 = 0, 1, 2
        scores = (self._M[i1, i2], self._G1[i1, i2], self._G2[i1, i2])
        last_score = total_score = max(scores)
        last_move = np.argmax(scores)
        logging.debug((last_score, last_move))
        while i1 > 0 or i2 > 0:
            #logging.debug((last_move, (self._M[i1, i2], self._G1[i1, i2], self._G2[i1, i2])))
            if last_move == move_M:
                a1 = seq1[i1-1]
                a2 = seq2[i2-1]
                raln1.append(a1)
                raln2.append(a2)
                m_score = sc.match if a1 == a2 else sc.mismatch
                i1 -= 1
                i2 -= 1
                matchable = i1 > 0 and i2 > 0
                last_score -= m_score
                if matchable and self._M[i1, i2] == last_score:
                    last_move = move_M
                elif self._G1[i1, i2] == last_score:
                    last_move = move_G1
                elif self._G2[i1, i2] == last_score:
                    last_move = move_G2
                else:
                    raise Exception('Invalid last move')
            elif last_move == move_G1:
                a2 = seq2[i2-1]
                raln1.append('-')
                raln2.append(a2)
                i2 -= 1
                matchable = i1 > 0 and i2 > 0
                if matchable and self._M[i1, i2] == last_score - sc.gap1_open:
                    last_score -= sc.gap1_open
                    last_move = move_M
                elif self._G1[i1, i2] == last_score - sc.gap1_ext:
                    last_score -= sc.gap1_ext
                    last_move = move_G1
                else:
                    raise Exception('Invalid last move')
            elif last_move == move_G2:
                a1 = seq1[i1-1]
                raln1.append(a1)
                raln2.append('-')
                i1 -= 1
                matchable = i1 > 0 and i2 > 0
                if matchable and self._M[i1, i2] == last_score - sc.gap2_open:
                    last_score -= sc.gap1_open
                    last_move = move_M
                elif self._G2[i1, i2] == last_score - sc.gap2_ext:
                    last_score -= sc.gap2_ext
                    last_move = move_G2
                else:
                    raise Exception('Invalid last move')
            else:
                raise Exception('undefined move')

        aln1 = ''.join(reversed(raln1))
        aln2 = ''.join(reversed(raln2))
        return {'score': total_score, 'aln1': aln1, 'aln2': aln2}


class DefaultReadMoveScore(object):
    def __init__(self, ins_prob=.1, del_prob=.1, mis_prob=.1):
        noins_prob = 1. - ins_prob
        assert noins_prob >= 0.
        self.ins_score = dict((x, np.log(ins_prob / 4.)) for x in 'ACGT')
        self.del_score = dict((x, np.log(noins_prob * del_prob)) for x in 'ACGT')
        self.trans_score = dict((x + x, np.log(noins_prob * (1. - mis_prob))) for x in 'ACGT')
        self.trans_score.update(dict((x + y, np.log(noins_prob * mis_prob / 3.)) for x in 'ACGT' for y in 'ACGT' if x != y))


class NWReadAlignment(NWAlignment):
    move_score = DefaultReadMoveScore()

    def __init__(self, tmpl_seq, read_seq, move_score=None):
        super(NWReadAlignment, self).__init__(tmpl_seq, read_seq, dtype=float, move_score=move_score)

    def _fill_matrix(self):
        cdef int len1, len2, i1, i2
        cdef float m, gap1, gap2
        cdef str r1, r2
        cdef np.ndarray M
        len1 = self.len1
        len2 = self.len2
        self._M = M = np.zeros((len1 + 1, len2 + 1), dtype=self.dtype)
        sc = self.move_score

        for i1 in xrange(1, len1 + 1):   # fill left
            r1 = self.seq1[i1 - 1]
            M[i1, 0] = M[i1 - 1, 0] + sc.del_score[r1]
        for i2 in xrange(1, len2 + 1):   # fill top
            r2 = self.seq2[i2 - 1]
            M[0, i2] = M[0, i2 - 1] + sc.ins_score[r2]

        for i1 in xrange(1, len1 + 1):
            r1 = self.seq1[i1-1]
            for i2 in xrange(1, len2 + 1):
                r2 = self.seq2[i2-1]
                m = sc.trans_score[r1 + r2]
                gap1 = sc.ins_score[r2]
                gap2 = sc.del_score[r1]
                M[i1, i2] = logsumexp(np.array((
                    M[i1-1, i2-1] + m,
                    M[i1, i2-1] + gap1,
                    M[i1-1, i2] + gap2,
                )))

    @property
    def score(self):
        return self._M[self.len1, self.len2]
