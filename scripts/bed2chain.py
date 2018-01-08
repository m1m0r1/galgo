#!/usr/bin/env python
from __future__ import print_function
from argtools import command, argument
import itertools
from operator import itemgetter, attrgetter
from collections import namedtuple

ChainHeader = namedtuple('ChainHeader', 'score tName tSize tStrand tStart tEnd qName qSize qStrand qStart qEnd id')


def bed_to_chains(bed, rlengths, inverse=False):
    """ Convert bed to chain file

    Args:
        bed: iterable of (chrom, start, end, length)
        rlenths : dictionary of {rname: length}
    Yields:
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

    >>> recs = [
    ...    LengthBed('chr1', 1000, 1001, 1),
    ...    LengthBed('chr1', 1002, 1002, 4),
    ...    LengthBed('chr1', 1010, 1010, 2),
    ...
    ...    LengthBed('chr2', 1000, 1010, 0),  # chrom changed
    ...    LengthBed('chr2', 1020, 1021, 1),
    ...    LengthBed('chr2', 1030, 1030, 8),
    ...    LengthBed('chr2', 1040, 1043, 8),
    ...    LengthBed('chr2', 2000, 3000, 3000),
    ...    LengthBed('chr2', 3500, 4000, 0),
    ... ]
    >>> rlengths = {'chr1': 10000, 'chr2': 30000, 'chr3': 10000}

    >>> chains = list(bed_to_chains(recs, rlengths))
    >>> attrgetter('score', 'tName', 'tSize', 'tStart', 'tEnd', 'qName', 'qSize', 'qStart', 'qEnd', 'id')(chains[0]['header'])
    (0, 'chr1', 10000, 0, 10000, 'chr1', 10006, 0, 10006, 1)
    >>> chains[0]['data']
    [(1002, 0, 4), (8, 0, 2), (8990,)]

    >>> attrgetter('score', 'tName', 'tSize', 'tStart', 'tEnd', 'qName', 'qSize', 'qStart', 'qEnd', 'id')(chains[1]['header'])
    (0, 'chr2', 30000, 0, 30000, 'chr2', 31503, 0, 31503, 2)
    >>> chains[1]['data']
    [(1000, 10, 0), (20, 0, 8), (10, 3, 8), (957, 1000, 3000), (500, 500, 0), (26000,)]

    >>> attrgetter('score', 'tName', 'tSize', 'tStart', 'tEnd', 'qName', 'qSize', 'qStart', 'qEnd', 'id')(chains[2]['header'])
    (0, 'chr3', 10000, 0, 10000, 'chr3', 10000, 0, 10000, 3)
    >>> chains[2]['data']
    [(10000,)]

    # inverse mode
    >>> chains = list(bed_to_chains(recs, rlengths, inverse=True))
    >>> attrgetter('score', 'tName', 'tSize', 'tStart', 'tEnd', 'qName', 'qSize', 'qStart', 'qEnd', 'id')(chains[0]['header'])
    (0, 'chr1', 10006, 0, 10006, 'chr1', 10000, 0, 10000, 1)
    >>> chains[0]['data']
    [(1002, 4, 0), (8, 2, 0), (8990,)]
    """
    rname_seen = set()
    chain_id = 1
    for (rname, recs) in itertools.groupby(bed, itemgetter(0)):
        assert rname not in rname_seen, "svbed are not sorted in rname! ({0} was appeared again)".format(rname)
        rname_seen.add(rname)
        tSize = qSize = rlengths[rname]
        data = []

        last_pos = 0
        for rec in recs:
            common_len = rec.start - last_pos
            ref_len = rec.end - rec.start
            alt_len = rec.length
            if ref_len == alt_len: # uncount the feature with same length
                continue
            elif ref_len > 0 or alt_len > 0:
                data.append((common_len, ref_len, alt_len))
                last_pos = rec.end
                qSize += (alt_len - ref_len)
            else:
                continue

        last_len = tSize - last_pos
        if inverse:
            tSize, qSize = qSize, tSize
            data = [(common_len, alt_len, ref_len) for common_len, ref_len, alt_len in data]

        data.append((last_len,))

        yield {'header': ChainHeader(0, rname, tSize, '+', 0, tSize,
                                        rname, qSize, '+', 0, qSize, chain_id),
               'data': data}

        chain_id += 1

    # remaining
    for rname in set(rlengths) - rname_seen:
        assert rname not in rname_seen
        rname_seen.add(rname)
        tSize = qSize = rlengths[rname]
        data = []
        last_pos = 0
        last_len = tSize - last_pos
        data.append((last_len,))
        yield {'header': ChainHeader(0, rname, tSize, '+', 0, tSize,
                                        rname, qSize, '+', 0, qSize, chain_id),
               'data': data}

        chain_id += 1


def iter_tabs(fname):
    with open(fname) as fp:
        for line in fp:
            row = line.rstrip('\r\n').split('\t')
            yield row

def read_genomes(fname):
    rlens = []
    for row in iter_tabs(fname):
        rlens.append((row[0], int(row[1])))
    return rlens

LengthBed = namedtuple('Bed', 'chrom start end length')
def iter_beds(fname, length_column=4):
    idx = length_column - 1
    for row in iter_tabs(fname):
        yield LengthBed(row[0], int(row[1]), int(row[2]), int(row[idx]))

@command
@argument('bed', help='bed file with "region replaced length"')
@argument('genomes', help='tsv file of rname and length as 1st and 2nd column, respectively (e.g. .fai file)')
@argument('-c', '--column', default=4, help='the column of "region replaced length" in bed')
@argument('-i', '--inverse', action='store_true')
def bed2chain(args):
    """ Convert bed file into chain file
    """
    rlens = read_genomes(args.genomes)

    # print chain file
    it = iter_beds(args.bed)
    for chains in bed_to_chains(it, dict(rlens), inverse=args.inverse):
        print('chain', *chains['header'], sep=' ')
        for row in chains['data']:
            print(*row, sep='\t')
        print()

if __name__ == '__main__':
    command.run()
