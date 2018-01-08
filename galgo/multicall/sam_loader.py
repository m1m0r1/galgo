from __future__ import print_function
from builtins import filter
from . import model
from collections import defaultdict
from itertools import product, groupby
import logging
import numpy as np
from ..samutil import parse_region, Cigar


class SamDepthAnalyzer(object):
    def __init__(self, sam, regions=None):
        self._sam = sam
        if regions is None:
            self._regions = None
        else:
            self._regions = [parse_region(reg) for reg in regions]

    def _iter_pileups(self):
        if self._regions is None:
            for pcol in self._sam.pileup():
                yield pcol
        else:
            for ref, s, e in self._regions:
                for pcol in self._sam.pileup(reference=ref, start=s, end=e):
                    yield pcol

    def iter_records(self):
        for pcol in self._iter_pileups():
            yield SamPosInfo(pcol)

    def iter_ref_records(self):
        for tid, pcols in groupby(self._iter_pileups(), lambda x: x.tid):
            rname = self._sam.getrname(tid)
            yield tid, rname, (SamPosInfo(pcol) for pcol in pcols)

    #def get_tid_refnames(self):
    #    return {self._sam.get_tid(ref): ref for ref in self._sam.references}


class SamPosInfo(object):
    '''
    contig
    pos
    nsegs        number of mapped segments
    nubases      number of unique bases
    ubases       unique bases
    unsegs       number of mapped segmnets for each base
    uquals_mean  the mean of qualty values for each base
    wsegs        multimap weight of mapped segments
    uwsegs       multimap weight of mapped segments for each base
    # bases        full list of bases
    # quals        full list of bases
    '''
    def __init__(self, pcol):
        self._pcol = pcol
        self._pileups = pileups = list(filter(lambda x: not x.alignment.is_unmapped, pcol.pileups))   # filter unmapped
        self.tid = pcol.tid
        self.pos = pcol.pos
        self._base_qual_lis = base_qual_lis = [_get_base_qual(x) for x in pileups]
        self.nsegs = len(pileups)

        # group by base
        self.ubases = []
        base_quals = {}
        base_reads = {}
        base_weights = {}
        for i in range(self.nsegs):
            b, q = base_qual_lis[i]
            if b not in base_quals:
                self.ubases.append(b)
                base_quals[b] = []
                base_reads[b] = []
                base_weights[b] = []
            base_quals[b].append(q)
            base_reads[b].append(pileups[i])
            base_weights[b].append(_get_weight(pileups[i]))

        self.nubases = len(self.ubases)
        self.unsegs = [len(base_reads[b]) for b in self.ubases]
        self.uquals_mean = [np.mean(sum(base_quals[b], ()) or 0) for b in self.ubases]
        self.uwsegs = [sum(base_weights[b]) for b in self.ubases]
        self.wsegs = sum(self.uwsegs)

    def get_read_info(self):
        return [{'pcol': pcol, 'bases': bases, 'quals': quals} for pcol, (bases, quals) in zip(self._pileups, self._base_qual_lis)]

        # TODO
        # bases, quals

def _get_weight(pread, tag='YC'):
    """ tag is 'number of candidates mapped'
    """
    aln = pread.alignment
    try:
        return 1. / (aln.get_tag(tag) + 1)   # number of candidates
    except KeyError:
        return 1.   # default weight

def _get_base_qual(pread):
    """
    Returns: (seq, (quals))

    Used flags:
        indel
            indel length for the position following the current pileup site.
            This quantity peeks ahead to the next cigar operation in this
            alignment. If the next operation is an insertion, indel will
            be positive. If the next operation is a deletion, it will be
            negation. 0 if the next operation is not an indel.
        is_del
            1 iff the base on the padded read is a deletion
    """
    if pread.is_del:  # current is del
        return '-', ()
    qpos = pread.query_position
    aln = pread.alignment
    if pread.indel > 0:  # next is ins
        sl = slice(qpos, qpos + pread.indel + 1)
        return (aln.query_sequence[sl], tuple(aln.query_qualities[sl]))
    return (aln.query_sequence[qpos], (aln.query_qualities[qpos],))


from tqdm import tqdm
class SamModelBuilder(object):
    def __init__(self, sam, regions=None):
        #it = islice(it, 20000)
        self.model = dm = model.DataModel()
        da = SamDepthAnalyzer(sam, regions=regions)
        refs = da.iter_ref_records()
        _frags = {}  # {name: Fragment}
        _alnss = defaultdict(list)   # {name: [Alignment]}
        _aln_variants = defaultdict(list)   # {aln: [(Variant, bases, quals)]}
        tid_refs = {}

        # 1. build alignment graph
        #   variant => [alignment]
        #   alignment => [variant]
        for tid, refname, it in refs:
            ref = model.Reference(refname)
            tid_refs[tid] = ref
            dm.refs.append(ref)
            logging.info('Storing variants for %s', refname)
            it = tqdm(it)
            for pos_info in it:
                if pos_info.nubases <= 1:   #  at least two alleles were needed
                    continue
                if list(sorted(pos_info.unsegs, reverse=True))[1] > 1: # second major alleles > 1
                    v = model.Variant(ref, pos_info.pos)
                    read_info = pos_info.get_read_info()
                    alns = [info['pcol'].alignment for info in read_info]
                    for aln, info in zip(alns, read_info):
                        _alnss[aln.qname].append(aln)
                        _aln_variants[aln].append((v, info['bases'], info['quals']))
            logging.info('Storing variants is done')

        def get_uniq_alns(alns):
            """ #TODO 
            - choose most consistent ones ?
            - the same alignment object should be reduced in advance by attach ids ?
            """
            ualns = []
            uvars = set()
            for aln in alns:
                key = tuple(sorted([(v.pos, bs) for v, bs, qs in _aln_variants[aln]]))
                if key not in uvars:
                    ualns.append(aln)
                    uvars.add(key)
            return ualns

        def find_aln_pairs(alns):
            """
            """
            tid_alns = defaultdict(lambda : {1: [], 2: []})
            # print ([a.tid for a in alns])
            # print ([a.cigarstring for a in alns])
            # print ([a.pos for a in alns])
            # print ([a.is_read2 for a in alns])
            # print ([a.get_tags()['XC'] for a in alns])
            for aln in alns:
                part = 2 if aln.is_read2 else 1
                tid_alns[aln.tid][part].append(aln)

            pairs = []
            for tid in tid_alns.keys():
                parts = tid_alns[tid]
                parts1 = parts[1]
                parts2 = parts[2]
                # remove redundants
                parts1 = get_uniq_alns(parts1)
                parts2 = get_uniq_alns(parts2)
                # parts1f = list(filter(lambda a: not Cigar(a.cigartuples).has_clip(), parts1))
                # parts2f = list(filter(lambda a: not Cigar(a.cigartuples).has_clip(), parts2))
                # if no clip aln were exist, discard clipped alns
                # parts1 = parts1f or parts1
                # parts2 = parts2f or parts2

                MAX_PARTS = 100
                if len(parts1) > MAX_PARTS or len(parts2) > MAX_PARTS or len(parts1) * len(parts2) > MAX_PARTS:   # FIXME: better truncation algorithm
                    refname = tid_refs[tid].name
                    logging.warning('Too many alignments for read 1: %s or read 2: %s mapped on %s', len(parts1), len(parts2), refname)
                    continue
                if parts1 and parts2:
                    for pair in product(parts1, parts2):
                        pairs.append((tid, pair))
                elif parts1:
                    for aln in parts1:
                        pairs.append((tid, (aln, None)))
                elif parts2:
                    for aln in parts2:
                        pairs.append((tid, (None, aln)))

            return pairs

        def make_fragment(name, alns):
            pairs = find_aln_pairs(alns)
            if not pairs:
                return
            frag = model.Fragment(name)
            for tid, (aln1, aln2) in pairs:
                ref = tid_refs[tid]
                cand = model.PairMapCandidate(frag, ref)
                if aln1:
                    cand.set_read1(aln1)
                    for v, bs, qs in _aln_variants[aln1]:
                        cand.add_variant(v, bs, qs)
                if aln2:
                    cand.set_read2(aln2)
                    for v, bs, qs in _aln_variants[aln2]:
                        cand.add_variant(v, bs, qs)

                for v in cand.variants.values():
                    v.add_candidate(cand)
            return frag

        # 2. build fragment graph
        #   variant => [fragment_candidate]
        #   fragment_candidate => [variant]
        logging.info('Creating fragments')
        for name, alns in tqdm(_alnss.items()):
            frag = make_fragment(name, alns)
            if frag is None:
                continue
            _frags[name] = frag

        self.model.fragments.extend(_frags.values())
