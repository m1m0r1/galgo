from builtins import filter, zip
from six import StringIO
import logging
from collections import namedtuple, defaultdict
from itertools import takewhile, dropwhile, chain, groupby
from .bioseq import dna_revcomp
from .samutil import aln2cigar
import pysam
import numpy as np
import pandas as pd

# Python 2 and 3:
#from io import BytesIO     # for handling byte strings
#from io import StringIO    # for handling unicode strings

class BlastTabFile:
    # Fields: query id, subject id, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score, gaps, query length, subject length, query seq, subject seq
    _field_keys = {
            'query id': 'qname',
            'subject id': 'sname',
            '% identity': 'identity',
            'alignment length': 'aln_len',
            'mismatches': 'mismatches',
            'gap opens': 'gap_opens',
            'q. start': 'qstart',
            'q. end': 'qend',
            's. start': 'sstart',
            's. end': 'send',
            'evalue': 'evalue',
            'bit score': 'bit_score',
            'gaps': 'gaps',
            'query length': 'qlen',
            'subject length': 'slen',
            'query seq': 'qseq',
            'subject seq': 'sseq',
            }

    def __init__(self, reader):
        # read header
        it = iter(dropwhile(lambda x: not x.startswith('# Fields'), reader))
        line = next(it)    # Fields
        self.fields = [f.strip() for f in line.split(':', 1)[1].split(',')]
        self.keys = [self._field_keys[f] for f in self.fields]  # blast key
        #self._convs = [self._types.get(k) for k in self.keys]
        BlastRecord = namedtuple('BlastRecord', self.keys)
        logging.info('blast record format: %s', BlastRecord)
        #class BlastRecord(namedtuple('BlastRecord', self.keys), BlastRecordMixin):
        #    pass

        self._Record = BlastRecord
        self._reader = it

    def iter_raw(self):
        for line in self._reader:
            if line.startswith('#'):
                continue
            rec = line.rstrip().split('\t')
            #rec = [rec if conv is None else conv(rec) for conv, rec in zip(self._convs, rec)]
            yield self._Record(*rec)

    def iter_query_tables(self, analyze=True):
        it = self._reader
        while 1:
            output = StringIO()
            it = takewhile(lambda x: not x.startswith('#'),
                    dropwhile(lambda x: x.startswith('#'), self._reader))
            for line in it:
                output.write(line)

            output.seek(0)
            tab = pd.read_table(output, names=self.keys, index_col=False)
            #logging.debug(tab)
            if not len(tab):
                break
            qname = tab['qname'][0]
            if analyze:
                self._analyze(tab)
            yield qname, tab
            output.close()

    @staticmethod
    def _analyze(tab):
        #logging.debug(tab)
        tab['edit_distance'] = tab['mismatches'] + tab['gaps']
        tab['qratio'] = 1. * ((tab['qend'] - tab['qstart']) + 1) / tab['qlen']
        tab['sratio'] = 1. * ((tab['send'] - tab['sstart']).abs() + 1) / tab['slen']
        tab['qlclip'] = tab['qstart'] - 1
        tab['qrclip'] = tab['qlen'] - tab['qend']
        tab['is_reverse'] = tab['send'] < tab['sstart']
        tab['slclip'] = np.where(tab['is_reverse'], tab['slen'] - tab['send'], tab['sstart'] - 1)
        tab['srclip'] = np.where(~tab['is_reverse'], tab['sstart'] - 1, tab['slen'] - tab['send'])
        tab['left_overhang'] = np.min([tab.qlclip, tab.slclip], axis=0)
        tab['right_overhang'] = np.min([tab.qrclip, tab.srclip], axis=0)

    def iter_sam_record(self, reference_ids, tags=None):
        mapq = 30
        bq = '+'
        it = self.iter_raw()

        # TODO secondary alignment flag
        for rec in it:
            is_reverse = int(rec.send) < int(rec.sstart)
            if is_reverse:
                qaln = dna_revcomp(rec.qseq)
                saln = dna_revcomp(rec.sseq)
                ref_start = int(rec.send) - 1
                lclip = int(rec.qstart) - 1
                rclip = int(rec.qend) - int(rec.qend)
            else:
                qaln = rec.qseq
                saln = rec.sseq
                ref_start = int(rec.sstart) - 1
                lclip = int(rec.qlen) - int(rec.qend)
                rclip = int(rec.qstart) - 1
            cigar = aln2cigar(qaln, saln)
            qseq = qaln.replace('-', '')
            if lclip > 0:
                cigar.prepend((cigar.H, lclip))
            if rclip > 0:
                cigar.append((cigar.H, rclip))

            try:
                edit = int(rec.mismatches) + int(rec.gaps)
            except:
                edit = -1

            a = pysam.AlignedSegment()
            a.query_name = rec.qname.encode('ascii')
            a.query_sequence = qseq.encode('ascii')
            a.reference_id = reference_ids[rec.sname]
            a.flag = (16 if is_reverse else 0)
            a.reference_start = ref_start
            a.mapping_quality = mapq
            a.cigar = cigar.values
            a.next_reference_id = -1
            a.next_reference_start = -1
            try:
                a.template_length = int(rec.slen)
            except:
                pass
            a.query_qualities = pysam.qualitystring_to_array(bq * len(qseq))
            a.tags = [
                    ("AS", float(rec.bit_score)), # alignment score
                    ("NM", edit),               # edit distance
                    ("ZE", float(rec.evalue)),  # E-value
            ]
            yield a
