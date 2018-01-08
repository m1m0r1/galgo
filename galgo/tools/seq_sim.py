from __future__ import print_function
from argtools import command, argument
from builtins import range, zip, filter
from collections import namedtuple
import logging
from .. import bioseq
from ..bioseq import dna_revcomp, Fasta
from ..samutil import parse_region
from tqdm import tqdm
import numpy as np

Fragment = namedtuple('Fragment', 'read1 read2')
class Read(namedtuple('Read', 'name seq qual')):
    def as_fastq(self):
        return '\n'.join([
            '@' + self.name,
            self.seq,
            '+',
            self.qual,
        ])

def get_qual_char(error_rate):
    return chr(33 + (- 10 * int(np.round(np.log10(error_rate)))))


DEFAULT_CHUNK_SIZE = 10000
def _gen_rand(fn):
    while 1:
        for v in fn():
            yield v

# TODO It is possible to use size option by using (chunk_size,) + size for size option of random function

def gen_integer(low, high=None, chunk_size=DEFAULT_CHUNK_SIZE):
    fn = lambda : np.random.random_integers(low, high=high, size=chunk_size)
    return _gen_rand(fn)

def gen_bernouille(p=0.5, chunk_size=DEFAULT_CHUNK_SIZE):
    q = 1. - p
    fn = lambda : np.random.choice([0, 1], size=chunk_size, p=[p, q])
    return _gen_rand(fn)

def gen_norm(mean=1, std=0, chunk_size=DEFAULT_CHUNK_SIZE):
    fn = lambda : np.random.normal(loc=mean, scale=std, size=chunk_size)
    return _gen_rand(fn)

def gen_poisson(lam=1.0, chunk_size=DEFAULT_CHUNK_SIZE):
    fn = lambda : np.random.poisson(lam=lam, size=chunk_size)
    return _gen_rand(fn)


class AlphaRot(object):
    """
    >>> arot = AlphaRot('ACGT')
    >>> arot.rot('A', 0)
    'A'
    >>> arot.rot('G', 3)
    'C'
    >>> arot.rot('T', 2)
    'C'
    >>> arot.seq_subst('ATCGATCGNN', [0, 8], [1, 2])
    'CTCGATCGNN'
    """
    def __init__(self, alphabet):
        self._index = {a: i for i, a in enumerate(alphabet)}
        self._as = alphabet
        self._len = len(alphabet)

    def rot(self, char, offset=0):
        i = self._index[char]
        return self._as[(i + offset) % self._len]

    def seq_subst(self, seq, poss, substs):
        for p, s in zip(poss, substs):
            char = seq[p]
            if char not in self._index:
                continue
            seq = seq[:p] + self.rot(char, s) + seq[p+1:]
        return seq


def _get_seq(seq, s, e, L):
    return ('N' * max(min(0, e) - s, 0)
            + seq[max(0, s) : max(min(e, L), 0)]
            + 'N' * max(e - max(s, L), 0))


class ReadPairSampler(object):
    """
    >>> sampler = ReadPairSampler('ACGTACGT', start=1, end=5)
    >>> sampler.get_seq(-5, 3)
    'NNNNNACG'
    >>> sampler.get_seq(-5, 12)
    'NNNNNACGTACGTNNNN'
    >>> sampler.get_seq(2, 7)
    'GTACG'
    >>> sampler.get_seq(-5, -1)
    'NNNN'
    >>> sampler.get_seq(9, 12)
    'NNN'

                           s      e
                           =======
    -------------->>>>>>>>>>>>>>>>>>>>>------------------
                  0                    L
                           =======
    --------------<<<<<<<<<<<<<<<<<<<<<------------------   (rev)
                 L        |      |    0
                          L-s    L-e
    """


    def __init__(self, seq, start=0, end=None, rlen=100, ilen=350, istd=70, error_rate=1e-3, prefix=''):
        self.start = start
        self.end = end = len(seq) if end is None else end
        self._seq_end = len(seq)
        self._len = end - start
        self.rlen = rlen
        self.ilen = ilen
        self.istd = istd
        self._seq = seq.upper()
        self._rev_seq = dna_revcomp(seq)
        self.error_rate = error_rate
        qual_char = get_qual_char(error_rate)
        self._qual = qual_char * rlen
        self._prefix = prefix

    def get_seq(self, s, e, reverse=False):
        L = self._seq_end
        if reverse:
            return _get_seq(self._rev_seq, L-e, L-s, L)
        else:
            return _get_seq(self._seq, s, e, L)

    def iter_samples(self, cov):
        rlen = self.rlen
        read_amount = np.random.poisson(int(1. * cov * self._len / (self.rlen * 2)))
        starts = gen_integer(self.start - self.ilen, self.end - 1)
        is_revs = gen_bernouille()
        ilens = gen_norm(self.ilen, self.istd)
        error_nums = gen_poisson(self.rlen * self.error_rate)
        error_poss = gen_integer(rlen-1)  # Note that this is only worked for fixed read length. Maybe it should sample wainting time for next variant.
        error_substs = gen_integer(1, 3)
        read_ids = iter(range(1, read_amount))

        # is_rev = 0
        # --1-->           <--2--
        # is_rev = 1
        # --2-->           <--1--
        arot = AlphaRot('ACGT')

        while 1:
            ilen = int(round(next(ilens)))
            if ilen < rlen:
                continue
            start = next(starts)
            is_rev = next(is_revs)
            end = start + ilen

            # seq1
            error_num = next(error_nums)
            error_poss1 = [next(error_poss) for i in range(error_num)]
            error_substs1 = [next(error_substs) for i in range(error_num)]

            # seq2
            error_num = next(error_nums)
            error_poss2 = [next(error_poss) for i in range(error_num)]
            error_substs2 = [next(error_substs) for i in range(error_num)]

            read_id = next(read_ids, None)

            if read_id:
                name = '{prefix}:{id}:{is_rev}'.format(prefix=self._prefix, id=read_id, is_rev=is_rev)
                name1 = name + '/1'
                name2 = name + '/2'
                logging.info((start, end, rlen))
                seq1 = self.get_seq(start, start + rlen)
                seq2 = self.get_seq(end - rlen, end, reverse=True)

                if error_poss1:
                    seq1 = arot.seq_subst(seq1, error_poss1, error_substs1)
                if error_poss2:
                    seq2 = arot.seq_subst(seq2, error_poss2, error_substs2)

                if is_rev:
                    seq1, seq2 = seq2, seq1

                read1 = Read(name1, seq1, self._qual)
                read2 = Read(name2, seq2, self._qual)
                f = Fragment(read1, read2)
                yield f
            else:
                break

@command.add_sub
@argument('fasta')
@argument('fastq1')
@argument('fastq2')
@argument('--regions', nargs='+')
@argument('-i', '--ilen', type=float, default=350)
@argument('-s', '--istd', type=float, default=70)
@argument('-r', '--rlen', type=int, default=100)
@argument('-x', '--coverage', type=float, default=20)
@argument('-e', '--error-rate', type=float, default=1e-3)
@argument('--seed', type=int, default=0)
@argument('--prefix', default='')
def gen_paired_reads(args):
    """
    """
    with open(args.fasta) as fa:
        fa = Fasta(fa)
        contigs = fa.contigs
    regions = args.regions or [c.name for c in contigs]
    np.random.seed(args.seed)
    logging.info('Ins len mean: %s', args.ilen)
    logging.info('Ins len std: %s', args.istd)
    logging.info('Read len: %s', args.rlen)
    logging.info('Coverage: %sx', args.coverage)
    logging.info('Error rate: %s', args.error_rate)
    logging.info('Random seed: %s', args.seed)

    with open(args.fastq1, 'w+') as read1, \
        open(args.fastq2, 'w+') as read2:

        for region in regions:
            rname, start, end = parse_region(region)
            contig = fa.get(rname)
            seq = contig.seq
            start = 0 if start is None else start
            end = len(seq) if end is None else end
            prefix = ':'.join((args.prefix, rname)) if args.prefix else rname
            sampler = ReadPairSampler(seq, start=start, end=end, rlen=args.rlen,
                    ilen=args.ilen, istd=args.istd, error_rate=args.error_rate, prefix=prefix)
            for frag in tqdm(sampler.iter_samples(args.coverage)):
                print (frag.read1.as_fastq(), file=read1)
                print (frag.read2.as_fastq(), file=read2)
