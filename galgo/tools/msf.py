from __future__ import print_function
from argtools import command, argument
from builtins import zip
from collections import namedtuple
from operator import attrgetter
import logging
from itertools import dropwhile
from ..utils import collect_while, isblank, blank_split, make_not
from ..bioseq import FastaContig


class MSF(object):
    """ Multiple Sequence Format

    ## example

    !!NA_MULTIPLE_ALIGNMENT

       MSF: 949  Type: N  Jan 15, 2015  12:21  Check: 0 ..

     Name: KIR3DP1*001      Len:   949  Check:  448  Weight:  1.00
     Name: KIR3DP1*002      Len:   949  Check:  324  Weight:  1.00
     Name: KIR3DP1*0030101  Len:   949  Check: 8583  Weight:  1.00
     Name: KIR3DP1*0030102  Len:   949  Check: 8583  Weight:  1.00
     Name: KIR3DP1*0030201  Len:   949  Check: 9501  Weight:  1.00
     ...
    //

        KIR3DP1*001  ATGTCGCTCA TGGTCGTCAG CATGGCGTGT GTTGGGTTCT TCTTGCTGCA
        KIR3DP1*002  ATGTCGCTCA TGGTCGTCAG CATGGCGTGT GTTGGGTTCT TCTTGCTGCA
    KIR3DP1*0030101  ATGTCGCTCA TGGTCGTCAG CATGGCGTGT GTTG...... ..........
    KIR3DP1*0030102  ATGTCGCTCA TGGTCGTCAG CATGGCGTGT GTTG...... ..........
    ...

    """
    def __init__(self, reader):
        self._reader = iter(reader)
        self._parsed = False

    def _parse(self):
        name_lines = []
        self._names = names = []
        #self._contigs = []
        self._contigs = {}
        it = self._reader
        is_name_line = lambda line: line.strip().startswith('Name:')
        it = dropwhile(make_not(is_name_line), it)
        name_lis, it = collect_while(is_name_line, it)

        for line in name_lis:
            name_line = line.strip()[len('Name: '):]   # name lines
            name_lines.append(name_line)
            row = blank_split(name_line)
            name = row[0]
            names.append(name)
            self._contigs[name] = FastaContig(name_line)

        it = dropwhile(lambda line: line.startswith('/') or line.strip() == '', it)
        while 1:
            try:
                seq_lis, it = collect_while(make_not(isblank), it)
            except StopIteration:
                break
            #if not self._contigs:
            #    for name_line in name_lines:
            #        contig = FastaContig(name_line)
            #        self._contigs.append(contig)

            for i, seq_line in enumerate(seq_lis):
                row = blank_split(seq_line)
                name = row[0]
                #print (name)
                #print (seq)
                seq = ''.join(row[1:])
                #self._contigs[i].append(seq)
                self._contigs[name].append(seq)

            it = dropwhile(isblank, it)
            #el = next(it, None)
            #if el is None:
            #    break

        self._parsed = True

    def __iter__(self):
        if not self._parsed:
            self._parse()
        #for contig in self._contigs:
        #    yield contig
        for name in self._names:
            yield self._contigs[name]


@command.add_sub
@argument('msf', help='all contigs should be the same lengths')
@argument('-m', '--missing-base', default='-', help='replace missing bases')
def msf2fa(args):
    """
    """
    if args.missing_base:
        assert len(args.missing_base) == 1, 'missing base should be a single character'

    with open(args.msf) as reader:
        msf = MSF(reader)
        for contig in msf:
            seq = contig.seq
            if args.missing_base:
                seq = seq.replace('.', args.missing_base)
            print ('>' + contig.name_line)
            print (seq)
