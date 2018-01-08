#!/usr/bin/env python

from __future__ import print_function
import re
import logging
from argtools import command, argument
from itertools import groupby, chain
from collections import namedtuple
from operator import attrgetter, itemgetter


@command
@argument('bed', nargs='?', default='/dev/stdin')
@argument('-s', '--bin-size', type=int, default=50)
@argument('-b', '--both-flank-bins', type=int, default=0)
@argument('-l', '--left-flank-bins', type=int, default=None)
@argument('-r', '--right-flank-bins', type=int, default=None)
@argument('--no-validate', dest='validate', action='store_false', default=True)
def binned_smooth(args):
    """ Smoothing each bedg value with an average value of flanking bins

    Currently, using gz file as input is slow. Following is recommended:
    zcat bedg.gz | <this program>
    """
    left = args.left_flank_bins
    right = args.right_flank_bins
    if left is None:
        left = args.both_flank_bins
    if right is None:
        right = args.both_flank_bins

    logging.info('Validation: %s', args.validate)
    logging.info('Bin size: %s bp', args.bin_size)
    logging.info('Left flank bins: %s (%s bp)', left, left * args.bin_size)
    logging.info('Right flank bins: %s (%s bp)', right, right * args.bin_size)

    smoother = BedgBinnedSmooth(bin_size=args.bin_size, left=left, right=right, validate=args.validate)
    it = read_bedg(args.bed)

    try:
        for (c, s, e, sig) in smoother(it):
            print (c, s, e, '{0:0.3f}'.format(sig), sep='\t')
    except AssertionError as e:
        logging.error('%s', e)
        return 1


BedgRec = namedtuple('BedgRecord', ('chrom', 'start', 'end', 'value'))
def read_bedg(file):
    re_skip = re.compile('#|browser|track')

    if file.endswith('.gz'):
        import gzip
        open = gzip.open
    else:
        import __builtin__
        open = __builtin__.open

    fp = open(file)
    try:
        line = next(fp)
        while re_skip.match(line):
            line = next(fp)

        while 1:
            row = line.rstrip().split()
            yield BedgRec(row[0], int(row[1]), int(row[2]), float(row[3]))
            line = next(fp)
    except StopIteration:
        pass
    finally:
        fp.close()


def with_validation(intervals, bin_size):
    """ Confirm chrom ordered, sorted in start position, and equally binned.
    """
    chrom_set_done = set()
    chrom = None
    start = -bin_size

    for (chrom, ivs) in groupby(intervals, itemgetter(0)):  # group by chrom
        assert chrom not in chrom_set_done, '{0} was appeared again!'.format(chrom)
        chrom_set_done.add(chrom)
        start = -bin_size
        for iv in ivs:
            assert start < iv.start, 'Not sorted in start position!'
            assert iv.end - iv.start == bin_size, 'Interval {0} != bin size: {1}'.format(iv.end - iv.start, bin_size)
            assert iv.start % bin_size == 0, 'Start position {0} should be divided by bin size: {1}'.format(iv.start, bin_size)
            start = iv.start
            yield iv


class BedgBinnedSmooth(object):
    """
    >>> data = [
    ...     BedgRec('chr1', 1000, 1050, 1),  # 20
    ...     BedgRec('chr1', 1050, 1100, 1),  # 21
    ...     BedgRec('chr1', 1100, 1150, 1),  # 22
    ...     BedgRec('chr1', 1200, 1250, 1),  # 24
    ...     BedgRec('chr1', 1950, 2000, 1),  # 39
    ...     BedgRec('chr1', 2000, 2050, 1),  # 40
    ...     BedgRec('chr3', 1000, 1050, 1),  # 20
    ...     BedgRec('chr3', 1050, 1100, 1),  # 21
    ...     BedgRec('chr2', 1000, 1050, 1),  # 20
    ...     BedgRec('chr2', 1050, 1100, 1),  # 21
    ... ]

    >>> for (c, s, e, v) in BedgBinnedSmooth(50)(data): print((c, s, e, '{0:0.3f}'.format(v)))
    ('chr1', 1000, 1050, '1.000')
    ('chr1', 1050, 1100, '1.000')
    ('chr1', 1100, 1150, '1.000')
    ('chr1', 1200, 1250, '1.000')
    ('chr1', 1950, 2000, '1.000')
    ('chr1', 2000, 2050, '1.000')
    ('chr3', 1000, 1050, '1.000')
    ('chr3', 1050, 1100, '1.000')
    ('chr2', 1000, 1050, '1.000')
    ('chr2', 1050, 1100, '1.000')

    >>> for (c, s, e, v) in BedgBinnedSmooth(50, left=1, right=2)(data): print((c, s, e, '{0:0.3f}'.format(v)))
    ('chr1', 900, 950, '0.250')
    ('chr1', 950, 1000, '0.500')
    ('chr1', 1000, 1050, '0.750')
    ('chr1', 1050, 1100, '0.750')
    ('chr1', 1100, 1150, '0.750')
    ('chr1', 1150, 1200, '0.500')
    ('chr1', 1200, 1250, '0.250')
    ('chr1', 1250, 1300, '0.250')
    ('chr1', 1850, 1900, '0.250')
    ('chr1', 1900, 1950, '0.500')
    ('chr1', 1950, 2000, '0.500')
    ('chr1', 2000, 2050, '0.500')
    ('chr1', 2050, 2100, '0.250')
    ('chr3', 900, 950, '0.250')
    ('chr3', 950, 1000, '0.500')
    ('chr3', 1000, 1050, '0.500')
    ('chr3', 1050, 1100, '0.500')
    ('chr3', 1100, 1150, '0.250')
    ('chr2', 900, 950, '0.250')
    ('chr2', 950, 1000, '0.500')
    ('chr2', 1000, 1050, '0.500')
    ('chr2', 1050, 1100, '0.500')
    ('chr2', 1100, 1150, '0.250')
    """

    def __init__(self, bin_size=50, left=0, right=0, validate=True):
        assert left >= 0 and right >= 0
        self.bin_size = bin_size
        self.left = left
        self.right = right
        self.bufsize = left + right
        self.denominator = 1. * (1 + left + right)
        self.validate = validate

    def __call__(self, intervals):
        denom = self.denominator
        bs = self.bin_size
        loffset = self.left  * bs
        roffset = self.right * bs
        if self.validate:
            intervals = with_validation(intervals, bs)

        for (chrom, ivs) in groupby(intervals, itemgetter(0)):  # group by chrom
            emit_start = 0
            buff = []
            sig = 0
            last_iv = BedgRec(chrom, float('inf'), float('inf'), 0.0)  # dummy record for avoid emission after ivs consumed

            for iv in chain(ivs, [last_iv]):
                while buff and emit_start + roffset < iv.start:
                    #sig = sum(v for (s, v) in buff)            # for speed, avoid summation
                    yield (chrom, emit_start, emit_start + bs, sig)
                    emit_start += bs
                    if buff[0][0] < emit_start - loffset:
                        sig -= buff[0][1]
                        buff.pop(0)

                val = iv.value / denom
                sig += val
                buff.append((iv.start, val))
                emit_start = max(emit_start, iv.start - roffset)


if __name__ == '__main__':
    exit(command.run())
