from __future__ import print_function
import logging
import numpy as np
from ..matrix import get_banded_matrix
from .basics import *


#TODO shadow class for doctest

class _NWAlignmentTrace(basics.NWAlignmentTrace):
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


class _NWAlignment(NWAlignment):
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


class _SWAlignment(SWAlignment):
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


class _NWAlignmentBandedNaive(NWAlignmentBandedNaive):
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


class _NWAlignmentBandedSparse(NWAlignmentBandedNaive):
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

#NWAlignmentBanded = NWAlignmentBandedSparse
#NWAlignmentBanded = NWAlignmentBandedNaive

class _NWGAlignment(NWGAlignment):
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
