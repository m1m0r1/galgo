import logging
import itertools
import numpy as np
from scipy.special import gammaln
from tqdm import tqdm
from collections import defaultdict
from ..utils import memoize, Counter

class Reference(object):
    """
    - name
    - variants: [Variant]
    - candidates: [MapCandidate]
    """
    def __init__(self, name):
        self.name = name
        self.variants = []
        self.candidates = []


class Variant(object):
    """
    - pos
    # - depth
    # - weights
    - candidates: [MapCandidate]   # candidate fragemnts
    """
    def __init__(self, ref, pos):
        self.ref = ref
        self.pos = pos
        self.candidates = []
        ref.variants.append(self)

    def add_candidate(self, cand):
        assert self.ref == cand.ref
        self.candidates.append(cand)

    def get_bases(self):
        """ [Bases]
        """
        return [cand.bases[self.pos] for cand in self.candidates]

    def get_quals(self):
        """ [Quals]
        """
        return [cand.quals[self.pos] for cand in self.candidates]


class MapCandidate(object):
    """
    - ref: Reference
    - refname: Reference.name
    - pos: [pos]
    - fragment: Fragment
    - variants: {pos: Variant}
    - bases: {pos: Bases}
    - quals: {pos: Quals}
    """
    is_paired = False

    def __init__(self, fragment, ref):
        self.ref = ref
        self.fragment = fragment
        self.pos = []
        self.variants = {}
        self.bases = {}
        self.quals = {}
        ref.candidates.append(self)
        fragment.candidates.append(self)

    def add_variant(self, variant, bases, quals):
        pos = variant.pos
        self.pos.append(pos)
        self.variants[pos] = variant
        self.bases[pos] = bases
        self.quals[pos] = quals


class PairMapCandidate(MapCandidate):
    """
    - is_paired: Bool
    - pos1
    - pos2
    - end1
    - end2
    - orientation   # FF | FR | RF | RR
    """
    has_read1 = False
    has_read2 = False

    @property
    def is_paired(self):
        return self.has_read1 and self.has_read2

    @property
    def orientation(self):
        return 'FR'[self.is_rev1] + 'FR'[self.is_rev2]

    def set_read1(self, aln):
        self.has_read1 = True
        self.pos1 = aln.pos
        self.end1 = aln.aend
        self.is_rev1 = aln.is_reverse

    def set_read2(self, aln):
        self.has_read2 = True
        self.pos2 = aln.pos
        self.end2 = aln.aend
        self.is_rev2 = aln.is_reverse


class Fragment(object):
    """
    - name     # read qname
    - candidates: [MapCandidate]
    """
    def __init__(self, name):
        self.name = name
        self.candidates = []

class DataModel(object):
    """
    - refs: [Reference]
    - fragments: [Fragment]
    """
    def __init__(self):
        self.refs = []
        self.fragments = []


def log2probs(log_vals):
    """
    >>> [round(v, 2) for v in log2probs(np.array([1., 1., 1., 1.]))]
    [0.25, 0.25, 0.25, 0.25]
    """
    vals = np.exp(log_vals - log_vals.max())
    return vals / vals.sum()


QC = .1 * np.log(10)
DEL_QUAL = 20
# simple
@memoize
def get_aln_score(bases, quals, true_bases):
    """
    get_aln_score('A', (10,), 'A')  # match
    get_aln_score('AT', (10, 20), '--')  # mismatch
    get_aln_score('-', (), 'AT')  # mismatch
    """
    if bases == '-':
        if bases == true_bases:
            return np.log(1 - 10 ** (- .1 * DEL_QUAL))
        else:
            return - QC * DEL_QUAL

    s = 0.
    l = min(len(bases), len(true_bases))
    s -= max(0, len(bases) - l) * QC * DEL_QUAL  # add score of missing bases for trial bases
    for i in xrange(l):
        if bases[i] == true_bases[i]:
            s += np.log(1 - 10 ** (- .1 * quals[i]))
        else:
            s += - QC * quals[i]
    return s


import networkx as nx
class VariantGraph(object):
    def __init__(self, refname):
        self.refname = refname
        self._graph = nx.MultiGraph()
        self._aln_vars = defaultdict(list)   # {aln_rec: [variant]}

    def add_variant(self, v):
        assert v.ref.name == self.refname
        self._graph.add_node(v.pos, variant=v)

    def add_edge(self, cand, v1, v2):
        assert v1.ref.name == self.refname
        assert v2.ref.name == self.refname
        self._graph.add_edge(v1.pos, v2.pos, **{'variants': (v1, v2), 'cand': cand})

    def get_graph(self):
        return self._graph

    def get_subgraphs(self):
        return list(nx.connected_component_subgraphs(self._graph))



class InferModel(object):
    """
    - data: DataModel
    - ref_phases: {refname: RefPhase}
    # - fragments: data.fragments
    - cands: [MapCandidate]          # sampling
    # - cand_phases: {cand: phase_id}  # phase_id

    - Z = {frag_id: [z_0, ..., z_Cn]}
    - H = {(ref_id, pos): [h_0, ..., ]}

    Make graph to traverse connected variants
    """
    def __init__(self, data_model, hap_depth, ploidies=None):
        """
        Args:
            data_model: DataModel
            mean_depth: Float
            ploidies:   {ref: ploidy}
        """
        self._dm = data_model
        self._ploidies = {ref.name: 2 for ref in self._dm.refs}
        if ploidies:
            self._ploidies.update(ploidies)

        self._part_phases = []
        self._base_cands = {}      # {var: [bases]}    # possible bases
        self._var_probs = defaultdict()    # {v: [[hap_base_prob]]}    # Q(h_gx)
        self._var_alloc_lscores = defaultdict(lambda : defaultdict(dict)) # {var: {cand: [lscore]}}        # scores for possible bases and read allocations

        self._alloc_var_depths = defaultdict(lambda : defaultdict(float)) # {alloc_code: {pos: mean_depth}}   # cached depth values
        self._alloc_probs = {}   # {frag.name: [prob.]}    # Q(z_n)
        self._alloc_codes = {}   # {frag.name: [0, 1, 2, 4, 5, 1, 2, ...]}   # value is a list of hap codes
        self._alloc_cands = {}   # {frag.name: [None, cand1, cand1, cand2, cand2, ...]}

        self._graphs = {}   # {refname: VariantGraph}

        self._ref_hap_codes = {}   # {ref_name: [hap_code1, ...]}
        i = 0
        for ref in self._dm.refs:
            codes = []
            for _ in range(self._ploidies[ref.name]):
                i += 1
                codes.append(i)
            self._ref_hap_codes[ref.name] = codes

        self._code_refs = list(self._dm.refs)                                          # {code: refname}
        self._hap_depth = hap_depth
        self._init_base_cands()
        self._init_var_alloc_scores()
        self._init_dist()
        self._init_het_graph()

    def _init_base_cands(self):
        # self._base_combis = {}      # {var: [(bases1, bases2, ..)]}    # possible genotypes
        logging.info('Init base combinations')
        for ref in self._dm.refs:
            rname = ref.name
            for v in tqdm(ref.variants):
                bases = v.get_bases()
                base_counts = Counter(bases)
                ubases = tuple(sorted(base_counts.keys()))
                self._base_cands[v] = ubases

    def _init_var_alloc_scores(self):
        logging.info('Init variant allocation scores')
        for ref in self._dm.refs:
            rname = ref.name
            for cand in tqdm(ref.candidates):
                for v in cand.variants.values():
                    pos = v.pos
                    hap_bases = self._base_cands[v]  # variant candidate bases
                    #print (cand.bases[pos], cand.quals[pos], hap_bases)
                    score = [get_aln_score(cand.bases[pos], cand.quals[pos], bases) for bases in hap_bases]
                    self._var_alloc_lscores[v][cand] = np.array(score)

    def _init_dist(self):
        # init all fragment assignments
        logging.info('Init all fragment assignments: Q(Z)')
        for frag in self._dm.fragments:
            fname = frag.name
            hap_codes = [0]   # null hap
            self._alloc_cands[fname] = [None]
            for c in frag.candidates:
                hap_codes.extend(self._ref_hap_codes[c.ref.name])
                self._alloc_cands[fname].extend([c] * self._ploidies[c.ref.name])
            self._alloc_codes[fname] = hap_codes
            probs = np.ones(len(hap_codes)) / len(hap_codes)
            self._set_frag_probs(fname, probs)

        logging.info('Setting initial variant probs: Q(H)')
        for ref in self._dm.refs:
            for v in ref.variants:
                hap_probs = self._eval_hap_probs(v)
                self._var_probs[v] = hap_probs

    def _list_variants(self, ref):
        vinfos = []
        for v in ref.variants:
            gen_scores = self.get_init_gen_scores(v)
            bases = [gs['bases'] for gs in gen_scores]
            probs = [gs['prob'] for gs in gen_scores]
            idx = np.argmax(probs)
            gt = bases[idx]
            sinfo = self.get_site_info(v)
            is_hom = len(set(gt)) == 1
            vinfos.append({'variant': v, 'is_hom': is_hom, 'gt': gt, 'site_score': sinfo['site_score'], 'prob': probs[idx]})
        return vinfos

    def _init_het_graph(self):
        logging.info('Creating HET graph')
        for ref in tqdm(self._dm.refs):
            rname = ref.name
            vg = VariantGraph(rname)
            self._graphs[rname] = vg
            var_infos = self._list_variants(ref)
            het_infos = [vi for vi in var_infos if not vi['is_hom']]
            het_set = set([vi['variant'] for vi in het_infos])
            for c in ref.candidates:
                c_hets = [c.variants[pos] for pos in c.pos if c.variants[pos] in het_set]
                for v in c_hets:
                    vg.add_variant(v)
                for v1, v2 in zip(c_hets, c_hets[1:]):
                    vg.add_edge(c, v1, v2)

    def _iter_overlap_allocs(self, v):
        """
        Yields: (Candidate, hap_idx)
        """
        for c in v.candidates:
            for idx in xrange(self._ploidies[c.ref.name]):
                yield (c, idx)

    def _eval_hap_probs(self, v):
        rname = v.ref.name
        log_probs = [np.zeros(len(self._base_cands[v])) for idx in range(self._ploidies[rname])]  # probs for base combis
        for cand, hap_idx in self._iter_overlap_allocs(v):  # no need to use this?
            frag = cand.fragment
            prob = self._alloc_probs[frag.name][hap_idx]
            lscore = self._var_alloc_lscores[v][cand]      # weights for base combis
            log_probs[hap_idx] += prob * lscore
        #return [log2probs(lp) for lp in log_probs]
        return log_probs

    def _set_frag_probs(self, name, probs):
        probs0 = self._alloc_probs.get(name)
        if probs0 is None:
            diffs = probs
        else:
            diffs = probs - probs0
        self._alloc_probs[name] = probs

        # update depth cache
        for i, code in enumerate(self._alloc_codes[name]):
            cand = self._alloc_cands[name][i]
            if cand is None:
                continue
            for pos in cand.pos:
                self._alloc_var_depths[code][pos] += diffs[code]

    @property
    def refs(self):
        return self._dm.refs

    @property
    def variants(self):
        vs = []
        for ref in self._dm.refs:
            vs.extend(ref.variants)
        return vs

    def get_site_info(self, v):
        """
        Args:
            v: Variant
        Returns:
            {hap_code: depth}
        """
        codes = self._ref_hap_codes[v.ref.name]
        ncands = len(v.candidates)
        site_score = ncands * np.log(self._hap_depth * 2) - gammaln(ncands + 1) - self._hap_depth * 2

        return {
            'ncands': len(v.candidates),
            'site_score': site_score,
            'base_counts': Counter(v.get_bases()),
            'mean_depths': {code: self._alloc_var_depths[code][v.pos] for code in codes},
        }

    def get_variant_probs(self, v):
        """
        Args:
            v: Variant
        Returns:
            [((b1, b2, ..), prob)]
        """
        rname = v.ref.name
        bases = tuple(self._base_cands[v])
        probs = self._var_probs[v]
        probs = [log2probs(ll) for ll in probs]

        return [{'bases': bases, 'prob': prob} for prob in probs]

    def get_het_blocks(self, refname):
        g = self._graphs[refname]
        subs = g.get_subgraphs()

        blocks = [] # {'pos': left_most_pos, 'variants': [Variant]}
        for sub in subs:
            nodes = sub.nodes(data=True)
            vs = [dat['variant'] for pos, dat in nodes]
            pos = min(v.pos for v in vs)
            blocks.append({'pos': pos, 'variants': vs})

        return blocks

    def get_init_gen_scores(self, v):
        """
        Args:
            v: Variant
        Returns:
            [{'bases': (b1, b2, ..), 'prob': prob}]
        """
        rname = v.ref.name
        # (0, 0), (0, 1), (1, 1)  # bases index
        basess = []
        lls = []

        # TODO non diploid samples
        assert self._ploidies[rname] == 2
        genotypes = []
        for i in xrange(len(self._base_cands[v])):
            for j in xrange(i+1):
                genotypes.append((i, j))

        for base_idxs in genotypes:
            bases = tuple(self._base_cands[v][bi] for bi in base_idxs)
            basess.append(bases)

            ll = 0
            for c in v.candidates:
                scores = self._var_alloc_lscores[v][c]
                ll += np.log(np.sum([np.exp(scores[bi]) for bi in base_idxs]))
            lls.append(ll)

        probs = log2probs(np.array(lls))
        ncands = len(v.candidates)
        site_score = ncands * np.log(self._hap_depth * 2) - gammaln(ncands + 1) - self._hap_depth * 2

        return [{'bases': bases, 'prob': prob, 'site_score': site_score, 'score': np.log(prob) + site_score} for bases, prob in zip(basess, probs)]

    def get_total_depth(self, v):
        """
        Args:
            v: Variant
        Returns:
            depth
        """
        return len(v.candidates)
