from builtins import zip, range
from collections import namedtuple
from operator import attrgetter
import logging
from .utils import color_term, cached_property, fill_text, lfill_text, chunked, iter_counts
import re
import numpy as np

_dna_rev_maps = dict(zip('ATCGatcg', 'TAGCtagc'))

def dna_revcomp(seq):
    """
    >>> dna_revcomp('ATG-CgcatN')
    'NatgcG-CAT'
    """
    return ''.join(_dna_rev_maps.get(c, c) for c in reversed(seq))

_colored_dnas = {
        'A': color_term('A', 'green'),
        'C': color_term('C', 'cyan'),
        'G': color_term('G', 'yellow'),
        'T': color_term('T', 'red'),
        'a': color_term('a', 'green'),
        'c': color_term('c', 'cyan'),
        'g': color_term('g', 'yellow'),
        't': color_term('t', 'red'),
}
def color_term_dna(seq):
    return ''.join(_colored_dnas.get(a, a) for a in seq)

# See http://www.bioinformatics.nl/~berndb/aacolour.html (LESK)
_colored_aas = {
        x: y
        for x, y in (
            [(x, color_term(x, 'yellow')) for x in 'GAST'] +      # small nonpolar
            [(x, color_term(x, 'green')) for x in 'CVILPFYMW'] +  # hydrophobic
            [(x, color_term(x, 'purple')) for x in 'NQH'] +       # polar
            [(x, color_term(x, 'red')) for x in 'DE'] +           # neg. charged
            [(x, color_term(x, 'cyan')) for x in 'KR']            # pos. charged
        )
}
def color_term_aa(seq):
    return ''.join(_colored_aas.get(a, a) for a in seq)

# build codon table
STOP_AA = '$'
_aa_rnas = {
        'F': ('UUU', 'UUC',),
        'L': ('UUA', 'UUG', 'CUU', 'CUC', 'CUA', 'CUG',),
        'I': ('AUU', 'AUC', 'AUA',),
        'M': ('AUG',),                                     # Methionine start aa
        'V': ('GUU', 'GUC', 'GUA', 'GUG',),
        'S': ('UCU', 'UCC', 'UCA', 'UCG', 'AGU', 'AGC',),
        'P': ('CCU', 'CCC', 'CCA', 'CCG',),
        'T': ('ACU', 'ACC', 'ACA', 'ACG',),
        'A': ('GCU', 'GCC', 'GCA', 'GCG',),
        'Y': ('UAU', 'UAC',),
        STOP_AA: ('UAA', 'UAG', 'UGA',),
        'H': ('CAU', 'CAC',),
        'Q': ('CAA', 'CAG',),
        'N': ('AAU', 'AAC',),
        'K': ('AAA', 'AAG',),
        'D': ('GAU', 'GAC',),
        'E': ('GAA', 'GAG',),
        'C': ('UGU', 'UGC',),
        'W': ('UGG',),
        'R': ('CGU', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG',),
        'G': ('GGU', 'GGC', 'GGA', 'GGG',),
}
_rna_aas = {}
_dna_aas = {}
for _aa, _rnas in _aa_rnas.items():
    for _rna in _rnas:
        _rna_aas[_rna] = _aa
        _dna_aas[_rna.replace('U', 'T')] = _aa

# https://en.wikipedia.org/wiki/Amino_acid
_AminoAcid = namedtuple('AminoAcid', 'name,name3,name1,side_chains,polarity,charge,h_index,mass')
amino_acids = [
        _AminoAcid('Alanine',        'Ala', 'A', ('aliphatic',),           'nonpolar',    0., 1.8, 89.094),
        _AminoAcid('Arginine',       'Arg', 'R', ('basic',),               'basic polar', 1., -4.5, 174.203),
        _AminoAcid('Asparagine',     'Asn', 'N', ('acid', 'amide'),        'polar',       0., -3.5, 132.119),
        _AminoAcid('Asparatic acid', 'Asp', 'D', ('acid', 'amide'),        'acidic polar', -1., -3.5, 133.104),
        _AminoAcid('Cysteine',       'Cys', 'C', ('sulfur-containing',),   'polar',       0., 2.5, 121.154),
        _AminoAcid('Glutamic acid',  'Glu', 'E', ('acid', 'amide'),        'acidic polar', -1., -3.5, 147.131),
        _AminoAcid('Glutamine',      'Gln', 'Q', ('acid', 'amide'),        'polar',       0., -3.5, 146.146),
        _AminoAcid('Glycine',        'Gly', 'G', ('aliphatic',),           'nonpolar',    0., -0.4, 75.067),
        _AminoAcid('Histidine',      'His', 'H', ('basic',),               'basic polar', .9, -3.2, 155.156),
        _AminoAcid('Isoleucine',     'Ile', 'I', ('aliphatic',),           'nonpolar',    0., 4.5, 131.175),
        _AminoAcid('Leucine',        'Leu', 'L', ('aliphatic',),           'nonpolar',    0., 3.8, 131.175),
        _AminoAcid('Lysine',         'Lys', 'K', ('basic',),               'basic polar', 1., -3.9, 146.189),
        _AminoAcid('Methionine',     'Met', 'M', ('sulfur-containing',),   'nonpolar',    0., 1.9, 149.208),
        _AminoAcid('Phenylalanine',  'Phe', 'F', ('aromatic',),            'nonpolar',    0., 2.8, 165.192),
        _AminoAcid('Proline',        'Pro', 'P', ('cyclic',),              'nonpolar',    0., -1.6, 115.132),
        _AminoAcid('Serine',         'Ser', 'S', ('hydroxyl-containing',), 'polar',       0., -0.8, 105.093),
        _AminoAcid('Theonine',       'Thr', 'T', ('hydroxyl-containing',), 'polar',       0., -0.7, 119.12),
        _AminoAcid('Tryptophan',     'Trp', 'W', ('aromatic',),            'nonpolar',    0., -0.9, 204.228),
        _AminoAcid('Tyrosine',       'Tyr', 'Y', ('aromatic',),            'polar',       0., -1.3, 118.191),
        _AminoAcid('Valine',         'Val', 'V', ('aliphatic',),           'nonpolar',    0., 4.2, 117.148),
]

def dna_translate(seq):
    """
    >>> dna_translate('ACTGGTTGCTAT')
    'TGCY'
    >>> dna_translate(iter('ACTGGTTGCTAT'))
    'TGCY'
    """
    return ''.join(iter_dna_translate(seq))

def iter_dna_translate(seq):
    """ seq should be uppercased in advance
    """
    return (_dna_aas[t1 + t2 + t3] for t1, t2, t3 in chunked(seq, 3))


def mask_same_bases(seq, ref, mask_char='.'):
    """
    >>> mask_same_bases('ATTG', 'AT-C')
    '..TG'
    """
    return ''.join((mask_char if s == r else s) for s, r in zip(seq, ref))

_pat_lclip = re.compile('^(-*)[^-]')
_pat_rclip = re.compile('.*[^-](-*)$')

def get_lclip(seq):
    """
    >>> get_lclip('---AG---')
    3
    >>> get_lclip('A---')
    0
    >>> get_lclip('---')
    0
    """
    m = _pat_lclip.match(seq)
    return len(m.group(1)) if m else 0

def get_rclip(seq):
    """
    >>> get_rclip('---AG---')
    3
    >>> get_rclip('A---')
    3
    >>> get_rclip('---')
    0
    """
    m = _pat_rclip.match(seq)
    return len(m.group(1)) if m else 0

_MASK_CHAR_DNA = '*'
_MASK_CHAR_AA = '.'
def get_default_mask_char(mode):
    if mode == 'dna':
        return _MASK_CHAR_DNA
    if mode == 'aa':
        return _MASK_CHAR_AA

_unknown_pat_dna = re.compile('-*[Nn]+-*')
_unknown_pat_aa = re.compile('-*[X]+-*')
def mask_unknown_seq(seq, mask_char=None, mode='dna'):
    """
    >>> aln = '----------AGGGNNG-TT'
    >>> aln
    '----------AGGGNNG-TT'

    >>> mask_unknown_seq(aln)
    '**********AGGG**G-TT'

    >>> aln = '----------AGGG-GTT---NN--NNNNNNN-NNNN--C-CCTA---------TTANNNN--NNNNNNNN-----'
    >>> aln
    '----------AGGG-GTT---NN--NNNNNNN-NNNN--C-CCTA---------TTANNNN--NNNNNNNN-----'

    >>> mask_unknown_seq(aln)
    '**********AGGG-GTT*********************C-CCTA---------TTA*******************'

    >>> mask_unknown_seq(aln.replace('N', 'X'), mode='aa')
    '..........AGGG-GTT.....................C-CCTA---------TTA...................'
    """
    prev_end = 0
    tokens = []
    assert mode in ('dna', 'aa')
    mask_char = mask_char or get_default_mask_char(mode)
    if mode == 'dna':
        upat = _unknown_pat_dna
        seq1 = ''.join(('N', seq, 'N'))
    elif mode == 'aa':
        upat = _unknown_pat_aa
        seq1 = ''.join(('X', seq, 'X'))

    for m in upat.finditer(seq1):
        start = m.start()
        end = m.end()
        tokens.append(seq1[prev_end:start])
        tokens.append(mask_char * (end - start))
        prev_end = end

    tokens.append(seq1[prev_end:])
    return ''.join(tokens)[1:-1]  # remove ends


_SeqPoint = namedtuple('SeqPoint', 'is_block index inner_pos')
_SeqInterval = namedtuple('SeqInterval', 'is_block index start end')

class SeqLocator(object):
    """
    pos is 0-based coordinate

    last pos
    0 1 1
    : : :
    |A|-|
                         1         2         3         4         5         6         7
    pos        0123456789012345678901234567890123456789012345678901234567890123456789012345
    is_block   0         1       0                    1                 0
    index      0         0       1                    1                 2
    block_id             0                            1
    last_end             01234456                     011234555555555567
                                                           1
    cum_last_end         01234456                     677890111111111123
    >>> aln = '**********AGGG-GTT*********************C-CCTA---------TTA*******************'
    >>> loc = SeqLocater(aln)
    >>> loc.blocks
    [(10, 18, 'AGGG-GTT'), (39, 57, 'C-CCTA---------TTA')]

    >>> loc.get_context(0)
    (None, None, '*')
    >>> loc.get_context(9, left=3, right=3)
    (None, None, '******')
    >>> loc.get_context(10, left=3, right=3)
    (0, 0, '***AGG')
    >>> loc.get_context(13, left=3, right=3)
    (0, 3, 'AGGGGT')
    >>> loc.get_context(14, left=3, right=3)
    (0, 4, 'GGGGTT')
    >>> loc.get_context(45, left=4, right=4)
    (1, 5, 'CCTATTA*')
    >>> loc.get_context(56, left=3, right=3)
    (1, 7, 'ATTA**')

    >>> [tuple(iv) for iv in loc.intervals]   # blank or block intervals
    [(0, 0, 0, 10), (1, 0, 10, 18), (0, 1, 18, 39), (1, 1, 39, 57), (0, 2, 57, 76)]

    >>> x = loc.get_point(9)
    >>> (x.is_block, x.index, x.inner_pos)
    (0, 0, 9)

    >>> x = loc.get_point(13)
    >>> (x.is_block, x.index, x.inner_pos)
    (1, 0, 3)
    """
    _mask_char = '*'
    _mask_pat = re.compile('\*+')
    def __init__(self, aln):
        self.aln = aln
        prev_end = 0
        blocks = []   # [(s, e, aln)]
        for m in self._mask_pat.finditer(aln):
            start = m.start()
            end = m.end()
            if start > prev_end:
                blocks.append((prev_end, start, aln[prev_end:start]))
            prev_end = end
        if prev_end < len(self.aln):
            blocks.append((prev_end, len(self.aln), aln[prev_end:]))
        self.blocks = blocks

        # add intervals
        self.intervals = ivs = []   # [SeqInterval]
        prev_end = 0
        for idx, b in enumerate(blocks):
            ivs.append(_SeqInterval(0, idx, prev_end, b[0]))  # not a block
            ivs.append(_SeqInterval(1, idx, b[0], b[1]))      # a block
            prev_end = b[1]
        ivs.append(_SeqInterval(0, len(blocks), prev_end, len(aln)))   # last blank

        self._pos_maps = []
        for s, e, aln in blocks:
            is_regular = np.array(list(aln)) != '-'
            pos_map = np.cumsum(is_regular) - is_regular
            self._pos_maps.append(pos_map)
            logging.debug(pos_map)

    def get_point(self, pos):   # linear search
        """
        Returns: SeqPoint
        """
        assert 0 <= pos < len(self.aln)
        for iv in self.intervals:
            if pos < iv.end:
                offset = pos - iv.start
                return _SeqPoint(iv.is_block, iv.index, offset)
        raise Exception('This should not be called')

    def get_block_index(self, pos):  # linear search
        assert 0 <= pos < len(self.aln)
        for idx, block in enumerate(self.blocks):
            if pos < block[0]:
                break
            if pos < block[1]:
                return idx
        return None

    def get_context(self, pos, left=0, right=1):
        """ Returns: (block_index, block_last-end, context)
        """
        assert left >= 0 and right >= 0
        idx = self.get_block_index(pos)
        if idx is None:
            return (None, None, self._mask_char * left + self._mask_char * right)
        start, end, aln = self.blocks[idx]
        pos_map = self._pos_maps[idx]
        offset = pos - start
        block_last_end = pos_map[offset]
        lseq = aln[:offset].replace('-', '')
        rseq = aln[offset:].replace('-', '')
        lseq1 = (self._mask_char * left + lseq)[-left:] if left > 0 else ''
        rseq1 = (rseq + self._mask_char * right)[:right]
        return (idx, block_last_end, lseq1 + rseq1)

# deprecated
SeqLocater = SeqLocator

def aln2pos(aln, offset=0, count_del=False):
    """
    >>> offset = 3
    >>> aln2pos(('AAAAA' 'TTTTT' 'GGGGG' 'CCCCC'), offset=offset)
    [3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22]

    >>> offset = 3
    >>> aln2pos(('AA--AAA' 'TTTTT' 'G--GGGG' 'CCCCC'), offset=offset)
    [3, 4, 4, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 13, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22]
    """
    j = offset - 1
    poss = []
    for i in range(len(aln)):
        b = aln[i]
        if count_del or b != '-':
            j += 1
        poss.append(j)
    return poss


def _iter_abbrev_nums(nums):
    yield nums[0]
    for n, n0 in zip(nums[1:], nums):
        for i in range(len(n)):
            if n[i] == n0[i]:
                continue
            break
        yield ' ' * i + n[i:]
        #yield ''.join((' ' if d == d0 else d) for d, d0 in zip(n, n0))


def pos2text(poss, vertical=False, skip_digit=False):
    """
    >>> pos2text([1, 10, 11, 12, 12, 15])
    [' 1', '10', '11', '12', '12', '15']
    >>> pos2text([1, 10, 11, 12, 12, 15])
    [' 1', '10', '11', '12', '12', '15']
    >>> pos2text([1, 10, 11, 12, 12, 15, 25, 121], skip_digit=True)
    ['  1', ' 10', '  1', '  2', '  2', '  5', ' 25', '121']
    >>> pos2text([1, 10, 11, 12, 12, 15], vertical=True)
    [' 11111', '101225']
    >>> pos2text([1, 10, 11, 12, 12, 15], vertical=True, skip_digit=True)
    [' 1    ', '101225']
    """
    txts = [str(pos) for pos in poss]
    max_len = max(len(t) for t in txts) if txts else 0
    filled = [lfill_text(t, max_len) for t in txts]
    if skip_digit:
        filled = list(_iter_abbrev_nums(filled))
    if vertical:
        return [''.join(ts) for ts in zip(*filled)]
    else:
        return filled


def get_aln_pos_text(aln, offset=0, digit=5, step=10, count_del=False):
    """
    Args:
        start 0-based coordinate

    >>> offset = 10010054
    >>> t = get_aln_pos_text(('AAAAA' 'TTTTT' 'GGGGG' 'CCCCC'), offset=offset)
    >>> t['number']
    '      10061     10071'
    >>> t['indicator']
    '      |         |   '

    >>> t = get_aln_pos_text(('AA--AAA' 'TTTTT' 'G--GGGG' 'CCCCC'), offset=offset)
    >>> t['number']
    '        10061       10071'
    >>> t['indicator']
    '        |           |   '

    >>> t = get_aln_pos_text(('AA--AAA' 'TTTTT' 'G--GGGG' 'CCCCC'), offset=offset, count_del=True)
    >>> t['number']
    '      10061     10071   '
    >>> t['indicator']
    '      |         |       '

    >>> t = get_aln_pos_text('AAAAA', offset=offset)
    >>> t['number']
    '10055'
    >>> t['indicator']
    '|    '
    """
    poss = []  # poss to print
    offsets = []
    assert step > digit
    j = left = offset % step
    for i in range(len(aln)):
        b = aln[i]
        if count_del or b != '-':
            j += 1
            if j % step == 1:
                poss.append(offset - left + j)
                offsets.append(i)
    if not poss:
        poss.append(offset + 1)
        offsets.append(0)
    offsets.append(len(aln))   # last offset

    n = 10 ** digit
    o = offsets[0]
    buf_n = [' ' * o]
    buf_i = [' ' * o]
    for p, o1, o2 in zip(poss, offsets, offsets[1:]):
        num = str(p % n)
        o = o2 - o1
        buf_n.append(fill_text(num, o))
        buf_i.append(fill_text('|', o))
    return {'number': ''.join(buf_n), 'indicator': ''.join(buf_i)}


ChainHeader = namedtuple('ChainHeader', 'score tName tSize tStrand tStart tEnd qName qSize qStrand qStart qEnd id')

def contig2chain(from_contig, to_contig, id=1):
    """ Convert bed to chain file

    Args:
        bed: iterable of (chrom, start, end, length)
        rlenths : dictionary of {rname: length}
    Return:
        {'header': header tokens
         'data': [(size, dt, dq), ... , (size,)}

    See http://genome.ucsc.edu/goldenPath/help/chain.html for more details.

    header tokens:
        score -- chain score
        tName -- chromosome (reference sequence)
        tSize -- chromosome size (reference sequence)
        tStrand -- strand (reference sequence)
        tStart -- alignment start position (reference sequence)
        tEnd -- alignment end position (reference sequence)
        qName -- chromosome (query sequence)
        qSize -- chromosome size (query sequence)
        qStrand -- strand (query sequence)
        qStart -- alignment start position (query sequence)
        qEnd -- alignment end position (query sequence)
        id -- chain ID

    data tokens (the end of data token contains only first token):
        size -- the size of the ungapped alignment
        dt -- the difference between the end of this block and the beginning of the next block (reference sequence)
        dq -- the difference between the end of this block and the beginning of the next block (query sequence)

    >>> reader = ('>contig1', 'AAA-TTT--ATTTTGT', '>contig2', 'A---CGTATATTTGG-')
    >>> fasta = Fasta(reader)
    >>> result = contig2chain(fasta.get('contig1'), fasta.get('contig2'))
    >>> attrgetter('score', 'tName', 'tSize', 'tStart', 'tEnd', 'qName', 'qSize', 'qStart', 'qEnd', 'id')(result['header'])
    (0, 'contig1', 13, 0, 13, 'contig2', 12, 0, 12, 1)
    >>> result['data']
    [[1, 2, 0], [3, 0, 2], [6, 1, 0], [0]]
    """
    from_size = len(from_contig.seq.replace('-', ''))
    to_size = len(to_contig.seq.replace('-', ''))
    header = ChainHeader(
            score=0,
            tName=from_contig.name, tSize=from_size, tStrand='+', tStart=0, tEnd=from_size,
            qName=to_contig.name, qSize=to_size, qStrand='+', qStart=0, qEnd=to_size,
            id=id,
    )
    data = []
    states = []
    def iter_states(seq1, seq2):
        for s1, s2 in zip(seq1, seq2):
            if s1 == s2 == '-':
                yield None
            elif s1 == '-':
                yield 2
            elif s2 == '-':
                yield 1
            else:
                yield 0

    it = iter_states(from_contig.seq, to_contig.seq)
    it = (el for el in it if el is not None)
    it = ((st, count) for st, count in iter_counts(it))

    dat = []
    for st, count in it:
        if st == 0:
            dat.append([count])
        elif st == 1:
            if not dat:
                dat.append([0])
            assert len(dat[-1]), 'Unexpected last data: {0} for state: {1}'.format(dat, st)
            dat[-1].extend((count, 0))
        elif st == 2:
            if not dat:
                dat.append([0])
            assert len(dat[-1]), 'Unexpected last data: {0} for state: {1}'.format(dat, st)
            dat[-1].extend((0, count))
        else:
            raise Exception('Unexpected state: {0}'.format(st))

    if not dat or len(dat[-1]) == 3:
        dat.append([0])

    return {'header': header, 'data': dat}


#TODO use .fai file if possible
class Fasta(object):
    """
    >>> reader = ('>contig1', 'AAATTT', '>contig2', 'ACGT')
    >>> fasta = Fasta(reader)
    >>> fasta.get('contig2').seq
    'ACGT'
    """

    def __init__(self, reader, load_now=False):
        self._reader = iter(reader)
        self._iterator_called = False
        if load_now:
            self.contigs   # load now

    @cached_property
    def contigs(self):
        assert self._iterator_called == False, 'Do not iterate before call the method'
        return list(self)

    @cached_property
    def _name_contigs(self):
        return dict((c.name, c) for c in self.contigs)

    def __contains__(self, name):
        return name in self._name_contigs

    def get(self, name):
        return self._name_contigs[name]

    @property
    def names(self):
        return [contig.name for contig in self.contigs]

    @property
    def lengths(self):
        return [len(contig) for contig in self.contigs]

    def __iter__(self):
        self._iterator_called = True
        contig = None
        try:
            while self._reader:
                line = next(self._reader).rstrip()
                if line.startswith('>'):
                    if contig is not None:
                        yield contig
                    contig = FastaContig(line[1:])
                else:
                    contig.append(line)
        except StopIteration:
            if contig is not None:
                yield contig


class FastaContig(object):
    """
    >>> contig = FastaContig('contig description')
    >>> contig.name_line
    'contig description'
    >>> contig.name
    'contig'
    >>> contig.append('AAAA')
    >>> contig.append('TTTT')
    >>> contig.seq
    'AAAATTTT'
    >>> contig.append('GGGG')
    >>> contig.seq
    'AAAATTTTGGGG'
    >>> len(contig)
    12
    >>> contig.lclip
    0

    >>> contig.append('----')
    >>> contig.seq
    'AAAATTTTGGGG----'
    >>> contig.rclip
    4

    >>> contig.get_clip_mask_seq()
    'AAAATTTTGGGG****'
    >>> contig.append('C-CC----')
    >>> contig.get_clip_mask_seq()
    'AAAATTTTGGGG----C-CC****'

    >>> contig.prepend('----')
    >>> contig.get_clip_mask_seq()
    '****AAAATTTTGGGG----C-CC****'
    """

    def __init__(self, name_line=''):
        assert not name_line.startswith('>'), 'please remove charcter ">"'
        self.name_line = name_line
        self.name = name_line.split(' ', 2)[0]
        self._seqs = []

    def prepend(self, seq):
        self._seqs.insert(0, seq)
        self.__dict__.pop('seq', None)

    def append(self, seq):
        self._seqs.append(seq)
        self.__dict__.pop('seq', None)

    @cached_property
    def seq(self):
        return ''.join(self._seqs)

    def __len__(self):
        return len(self.seq)

    def get_clip_mask_seq(self, mask_char='*'):
        seq = mask_char * self.lclip + self.seq[self.lclip : len(self) - self.rclip] + mask_char * self.rclip
        return seq

    @property
    def lclip(self):
        return get_lclip(self.seq)

    @property
    def rclip(self):
        return get_rclip(self.seq)


class AlignmentAnalyzer(object):
    """
    0-based pos

    pos         01234567
    last_end1   01234456
    last_end2   01233345
    >>> aln1 = 'AAAA-AAC'
    >>> aln2 = 'AAA--ACT'
    >>> aa = AlignmentAnalyzer(aln1, aln2)
    >>> list(aa.iter_edits())
    [(3, 'A', '-'), (6, 'A', 'C'), (7, 'C', 'T')]

    >>> keys = 'last_end1 ctx1 last_end2 ctx2'.split(' ')
    >>> ctx = aa.get_context(3, left=4, right=4)
    >>> [ctx[k] for k in keys]
    [3, '*AAAAAAC', 3, '*AAAACT*']

    >>> ctx = aa.get_context(6, left=4, right=4)
    >>> [ctx[k] for k in keys]
    [5, 'AAAAAC**', 4, 'AAAACT**']
    """

    def __init__(self, aln1, aln2, validate=True):
        if validate:
            assert len(aln1) == len(aln2)
            assert aln1 == aln1.upper()
            assert aln2 == aln2.upper()
        self.aln1 = aln1
        self.aln2 = aln2
        _aln1 = np.array(list(aln1))
        _aln2 = np.array(list(aln2))
        self._is_edit = (_aln1 != _aln2)

    @cached_property
    def _loc1(self):
        return SeqLocator(self.aln1)

    @cached_property
    def _loc2(self):
        return SeqLocator(self.aln2)

    def iter_edits(self):
        """ (pos_msa, aln1, aln2)
        """
        for pos in np.argwhere(self._is_edit).flatten():
            yield (pos, self.aln1[pos], self.aln2[pos])

    def get_context(self, pos, left=0, right=1):
        ctxs = {}
        ctxs['block_idx1'], ctxs['last_end1'], ctxs['ctx1'] = self._loc1.get_context(pos, left=left, right=right)
        ctxs['block_idx2'], ctxs['last_end2'], ctxs['ctx2'] = self._loc2.get_context(pos, left=left, right=right)
        return ctxs
