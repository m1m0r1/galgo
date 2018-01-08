from __future__ import print_function
from builtins import zip, range
import logging
import itertools
from ..utils import cached_property


class RefDB(object):
    def __init__(self, ref_db):
        self._db = ref_db

    def get_haps(self, rname):
        return []

    def get_contigs(self, refname):
        # reuqire hashmap of contigs for refname
        fasta = self._refs[refname]
        with open(fasta) as fp:
            fa = Fasta(fp)
            return fa.contigs


USE_PHASE_PROBS = True

class RefVarCall(object):
    def __init__(self, tab):
        self.copy_number = len(tab.iloc[0]['hap_bases'].split('|'))
        self._blocks = []   # [{'block_id': block_id, pos: [pos], phase_bases: [[base]]}]
        self._pos_infos = {}  # {pos: {'block_id': block_id, 'phase_bases': [[base]]}}

        for block_id, tab1 in tab.groupby('block_id'):
            phase_bases = [[] for _ in range(self.copy_number)]
            phase_probs = []
            zygosities = []
            if 'end' in tab1.columns:
                end_data = tab1['end']
                pos_list = list(pos1 for pos, end in zip(tab1['pos'], end_data) for pos1 in range(pos, end))
            else:
                end_data = tab1['pos'] + 1
                pos_list = list(tab1['pos'])

            if 'phase_prob' in tab1.columns:
                org_phase_probs = tab1['phase_prob']
            else:
                org_phase_probs = list(None for _ in range(len(tab1)))

            for gt, phase_prob in zip(tab1['hap_bases'], org_phase_probs):
                pos_phase_bases = [[] for _ in range(self.copy_number)]
                col_bases = gt.split('|')

                # scan rows
                for i, bases in enumerate(col_bases):
                    for base in bases.split('.'):
                        pos_phase_bases[i].append(base)

                # scan columns
                first_het = True
                for col_bases in zip(*pos_phase_bases):
                    bset = set(col_bases)
                    zygosities.append(len(bset))
                    if 'N' in bset:
                        assert len(bset) == 1, 'if N is in a column, all the bases at the column should be N'

                    if len(bset) >= 2:  # HET
                        if first_het:
                            if USE_PHASE_PROBS and phase_prob is not None:
                                phase_prob = float(phase_prob)

                            #phase_prob = min(phase_prob, .95)  # adhoc
                            #phase_prob = min(phase_prob, .5)  # adhoc
                            phase_probs.append(phase_prob)
                            first_het = False
                        else:
                            phase_probs.append(1.)
                    else:
                        phase_probs.append(None)

                for i, bases in enumerate(pos_phase_bases):
                    phase_bases[i].extend(bases)

            pos_bases = list(zip(*phase_bases))   # transpose phase_bases

            block_id = None if block_id == 'None' else block_id  # TODO fix input value
            logging.info('Number of positions for block %s: %s', block_id, len(pos_list))
            for i, hap in enumerate(phase_bases):
                logging.info('Number of bases in haplotype %s: %s', i, len(hap))
            assert all(len(pos_list) == len(hap) for hap in phase_bases)

            self._blocks.append({'block_id': block_id, 'pos': pos_list, 'pos_bases': pos_bases, 'phase_bases': phase_bases, 'phase_probs': phase_probs, 'zygosities': zygosities})
            for pos, end, gt in zip(tab1['pos'], end_data, tab1['hap_bases']):
                #has_call = block_id is None or last_step > 0  # TODO fix input value?
                col_bases = gt.split('|')
                hap_bases = [[] for _ in col_bases]
                # scan rows
                for i, bases in enumerate(col_bases):
                    for base in bases.split('.'):
                        hap_bases[i].append(base)

                for j, pos in enumerate(range(pos, end)):   # iterate per position
                    bases = [bs[j] for bs in hap_bases]
                    self._pos_infos[pos] = {'block_id': block_id, 'bases': bases} #, 'has_call': has_call}

        self._binfo = BlockInfo(self._blocks)
        self._validate()

    def _validate(self):
        hetd = self._binfo.get_het_trellis()
        homd = self._binfo.get_hom_data()
        assert len(self._pos_infos) == len(hetd) + len(homd)
        for d in hetd:
            pos = d['pos']
            gt1 = self._pos_infos[pos]['bases']
            gt2 = [d['uniq_bases'][idx] for idx in d['states'][0]]
            assert gt1 == gt2, (gt1, gt2)
        for d in homd:
            pos = d['pos']
            gt1 = self._pos_infos[pos]['bases']
            gt2 = [d['uniq_bases'][0] for _ in range(self.copy_number)]
            assert gt1 == gt2, (gt1, gt2)

    def get_pos_info(self, pos):
        return self._pos_infos.get(pos)

    def iter_blocks(self):
        for rec in self._blocks:
            yield rec

    def iter_poss(self):
        return sorted(pos for b in self._blocks for pos in b['pos'])

    def match_haps(self, haps, score='edit'):
        """
        Args:
            haps: [Contig]   # selected combinations of contigs
        """
        if score == 'hmm':
            matcher = MatchInfoHMM(self._binfo, self.copy_number, haps)
        else:
            matcher = MatchInfo(self._blocks, self.copy_number, haps, score=score)  # using `blocks` is confusing
        return matcher


class BlockInfo(object):
    def __init__(self, blocks):
        """
        Args:
            blocks: {'block_id': block_id, 'pos': pos_list, 'phase_bases': phase_bases, 'phase_probs': probs}  # TODO assign phase_conf

        Properties:
        - het_trellis: [{
                pos: 155,
                states: [(0, 0, 1), (0, 1, 0), (1, 0, 0)]   # if not multiple, this is homozygous call
                bases: ['AC', 'AT']   # bases for base_idx
                trans: [[.9, .05, .05], [.05, .9, .05], [.05, .05, .9]]   # matrix
            }]
        - hom_data: [{
                pos: 156,
                bases: ['AC']   # bases for base_idx
            }]
        """
        hetd = []
        self._homd = []
        for bl in blocks:
            pos_list = bl['pos']
            pos_bases = bl['pos_bases']
            phase_bases = bl['phase_bases']
            zygosities = bl['zygosities']
            phase_probs = bl['phase_probs']
            assert len(pos_list) == len(pos_bases) == len(zygosities) == len(phase_probs)
            for pos, bases, z, phase_prob in zip(pos_list, pos_bases, zygosities, phase_probs):
                uniq_bases = list(sorted(set(bases)))
                if z == 1:   # hom
                    self._homd.append({
                        'pos': pos,
                        'uniq_bases': uniq_bases,
                    })
                else:
                    first_state = [uniq_bases.index(bs) for bs in bases]
                    states = list(itertools.permutations(first_state)) # length is copy_number! and
                    lower_bound = 1. / len(states)
                    if phase_prob is None:
                        phase_prob = lower_bound
                    else:
                        phase_prob = float(phase_prob)
                    assert phase_prob >= lower_bound, 'phase_prob: {} should be larger than lowerbound: {}'.format(phase_prob, lower_bound)
                    hetd.append({
                        'pos': pos,
                        'uniq_bases': uniq_bases,
                        'states': states,
                        'phase_prob': phase_prob,
                    })

        self._hetd = list(sorted(hetd, key=lambda d: d['pos']))
        for d1, d2 in zip(self._hetd, self._hetd[1:]):
            dim1 = len(d1['states'])
            dim2 = len(d2['states'])
            states1 = d1['states']
            states2 = d2['states']
            phase_prob = d2['phase_prob']
            rem_len = len(states2) - 1
            rem_prob = (1 - phase_prob) / rem_len
            assert dim1 == dim2
            trans = np.diagflat([phase_prob - rem_prob] * dim1) + rem_prob * np.ones((dim1, dim2))
            assert trans.shape == (dim1, dim2)
            d2['trans'] = trans

    def get_het_trellis(self):
        return self._hetd

    def get_hom_data(self):
        return self._homd


class MatchInfo(object):
    """
        'score': sum of scores
        'edit': sum of edits
        'match': sum of matches
        'blocks': [{
            'block_id': block_id,
            'indexes': (index,),  # selected order e.g. (0,) or (1,) for one hap, (0, 1) or (1, 0) for two haps
            'match': match,  # sum of matches
            'edit': edit, # sum of edits
            'matches': [match], # match for each (phased, hap)
            'edits': [edit], # edit for each (phased, hap)
            'edit_details': [[(pos, ref_base, alt_base)]],   # each edit
            'score': # total score,
            'score_mode': # total score,
        }]
    """
    def __init__(self, input_blocks, copy_number, haps, score='edit'):
        if score == 'edit':
            def get_score(match_info):
                return - match_info['edit']
        elif score == 'match':
            def get_score(match_info):
                return match_info['match']
        else:
            raise NotImplementedError

        # TODO HMM score using phase confidence
        self._matches = matches = []
        indexes_cands = list(itertools.permutations(tuple(range(copy_number)), len(haps)))   # selection of blocks for each hap
        for bl in input_blocks:
            match_info_cands = [self._get_match_info(haps, bl, indexes) for indexes in indexes_cands]
            match_info = getmax(match_info_cands, key=get_score)
            match_info['score'] = get_score(match_info)
            matches.append(match_info)

        self.score = sum(b['score'] for b in matches)
        self.edit = sum(b['edit'] for b in matches)
        self.match = sum(b['match'] for b in matches)

    def get_block_matches(self):
        """
        Returns:
            [{
                'block_id': block_id,
                'indexes': indexes,
                'match': sum(matches),
                'edit': sum(edits),
                'matches': matches,
                'edits': edits,
                'edit_details': [(pos, ref, alt)],
            }]
        """
        return self._matches

    @staticmethod
    def _get_match_info(haps, block, indexes):
        #assert len(haps) == len(indexes)
        # these data are added in order of haps
        edits = []
        matches = []
        edit_details = []
        for hap, index in zip(haps, indexes):
            hap_bases = [hap[pos] for pos in block['pos']]
            bases = block['phase_bases'][index]
            # TODO calculate edit distance
            is_edits = [b1 != 'N' and b2 != 'N' and b1 != b2 for b1, b2 in zip(hap_bases, bases)]
            match = sum([b1 != 'N' and b2 != 'N' and b1 == b2 for b1, b2 in zip(hap_bases, bases)])
            details = []
            edit = 0
            for i, is_edit in enumerate(is_edits):
                if is_edit:
                    pos = block['pos'][i]
                    ref = hap[pos]
                    alt = bases[i]
                    details.append((pos, ref, alt))
                    edit += 1
            edits.append(edit)
            matches.append(match)
            edit_details.append(details)

        return {
            'block_id': block['block_id'],
            'match': sum(matches),
            'edit': sum(edits),
            'matches': matches,
            'edits': edits,
            'edit_details': edit_details,
        }


MISMATCH_PROB = .001
MISMATCH_PROB = .0001
#MISMATCH_PROB = .1

class MatchInfoHMM(object):
    def __init__(self, binfo, copy_number, haps, delta=MISMATCH_PROB, block_threshold=.95):
        self.delta = delta
        self._haps = haps
        self._hmm_calc = None
        self._copy_number = copy_number
        self._block_threshold = block_threshold

        LD0 = np.log(delta)
        LD1 = np.log(1 - delta)
        LDN = np.log(1./4)   # this will result in selecting falsy candidate
        LDN = LD1
        def get_lemit(bases, hap_bases):  # len(bases) >= len(hap_bases)
            #return sum([LD1 if (b == 'N' or hb == 'N' or b == hb) else LD0 for b, hb in zip(bases, hap_bases)]) # TODO single N is ok?
            return sum([LDN if (b == 'N' or hb == 'N') else LD1 if (b == hb) else LD0 for b, hb in zip(bases, hap_bases)]) # TODO single N is ok?

        # hom_emissions = #
        self._binfo = binfo
        self._hetd = hetd = self._binfo.get_het_trellis()    # linit, ltrans
        if hetd:
            l = len(hetd[0]['states'])
            linit = np.log(np.ones(l) / l)
            dims = []
            het_lemit = []
            ltrans = []
            for d in hetd:
                pos = d['pos']
                uniq_bases = d['uniq_bases']
                hap_bases = [hap[pos] for hap in haps]
                het_lemit.append(np.array([get_lemit([uniq_bases[idx] for idx in idxs], hap_bases) for idxs in d['states']]))
                trans = d.get('trans')
                if trans is not None:
                    trans = np.maximum(trans, HMMCalc.min_prob)  # TODO adhoc?
                    trans = np.minimum(trans, 1 - HMMCalc.min_prob)  # TODO adhoc?
                    ltrans.append(lnormal(np.log(trans), axis=1))
                dims.append(len(d['states']))
            #ltrans # head is not used

            self._hmm_calc = HMMCalc(dims, linit=linit, ltrans=ltrans, lemit=het_lemit)

        self._homd = homd = self._binfo.get_hom_data()
        hom_lemit = []
        for d in homd:
            pos = d['pos']
            uniq_bases = d['uniq_bases']
            hap_bases = [hap[pos] for hap in haps]
            lemit = get_lemit([uniq_bases[0] for _ in range(copy_number)], hap_bases)
            hom_lemit.append(lemit)

        if self._hmm_calc:
            self.score = self._hmm_calc.loglikelihood + sum(hom_lemit, 0)
        else:
            self.score = sum(hom_lemit, 0)

    @property
    def edit(self):
        return sum([bm['edit'] for bm in self._block_matches])

    @property
    def match(self):
        return sum([bm['match'] for bm in self._block_matches])

    def get_block_matches(self, prob_threshold=.9):
        # TODO
        """
        Returns:
            [{
                'block_id': block_id,
                'match': sum(matches),
                'edit': sum(edits),
                'matches': matches,
                'edits': edits,
                'edit_details': [(pos, ref, alt)],
            }]
        """
        return self._block_matches
        #if self._hmm_calc:
        #    states = self._hmm_calc.get_best_path() # TODO

    @cached_property
    def _block_matches(self):
        logging.info('Search block matches')
        def get_result(poss, all_hap_bases, all_typed_bases):
            assert [len(poss)] * self._copy_number == [len(x) for x in all_typed_bases] == [len(x) for x in all_hap_bases]
            edits = []
            matches = []
            edit_details = []
            for hap_bases, bases in zip(all_hap_bases, all_typed_bases):
                is_edits = [b1 != 'N' and b2 != 'N' and b1 != b2 for b1, b2 in zip(hap_bases, bases)]
                match = sum([b1 != 'N' and b2 != 'N' and b1 == b2 for b1, b2 in zip(hap_bases, bases)])
                details = []
                edit = 0
                for i, is_edit in enumerate(is_edits):
                    if is_edit:
                        pos = poss[i]
                        ref = hap_bases[i]
                        alt = bases[i]
                        assert ref != alt
                        details.append((pos, ref, alt))
                        edit += 1
                edits.append(edit)
                matches.append(match)
                edit_details.append(details)
            block_id = min(poss) if poss else None
            return {
                'block_id': block_id,
                'match': sum(matches),
                'edit': sum(edits),
                'matches': matches,
                'edits': edits,
                'edit_details': edit_details,
                'raw': {
                    'pos_list': poss,
                    'all_typed_bases': all_typed_bases,
                    'all_hap_bases': all_hap_bases,
                }
            }

        def get_hom_block():
            poss = []
            all_typed_bases = [[] for _ in range(self._copy_number)]
            block_id = None
            for d in self._homd:
                pos = d['pos']
                uniq_bases = d['uniq_bases']
                poss.append(pos)
                for j in range(self._copy_number):
                    all_typed_bases[j].append(uniq_bases[0])  # hom
            all_hap_bases = [[hap[pos] for pos in poss] for hap in self._haps]
            return get_result(poss, all_hap_bases, all_typed_bases)

        # Get single block for all the HET variants
        def get_het_block():
            poss = []
            all_typed_bases = [[] for _ in range(self._copy_number)]
            het_states = self._hmm_calc.best_states   # TODO it had better to use likely path
            assert len(het_states) == len(self._hetd)
            for i, d in enumerate(self._hetd):
                pos = d['pos']
                uniq_bases = d['uniq_bases']
                state_idx = het_states[i]
                state = d['states'][state_idx]
                poss.append(pos)
                #logging.info((i, pos, uniq_bases, state_idx, state, self._hmm_calc.state_lls[i], self._hmm_calc._lemit[i]))
                for j, s in enumerate(state):
                    all_typed_bases[j].append(uniq_bases[s])
            all_hap_bases = [[hap[pos] for pos in poss] for hap in self._haps]
            return get_result(poss, all_hap_bases, all_typed_bases)

        def iter_het_blocks():
            #vres = self._hmm_calc.viterbi
            vres = self._hmm_calc.viterbi_marginal
            next_probs = vres['state_probs'][1:] + [1.]
            it = iter(zip(vres['path'], next_probs, self._hetd))
            while 1:
                poss = []
                all_typed_bases = [[] for _ in range(self._copy_number)]
                for state_idx, next_prob, d in it:
                    pos = d['pos']
                    uniq_bases = d['uniq_bases']
                    state = d['states'][state_idx]
                    for j, s in enumerate(state):
                        all_typed_bases[j].append(uniq_bases[s])

                    poss.append(pos)
                    if next_prob < self._block_threshold:
                        break

                if not poss:
                    break
                logging.info('Phase block %s - %s (sites: %s, next_prob: %s)', poss[0], poss[-1], len(poss), next_prob)
                all_hap_bases = [[hap[pos] for pos in poss] for hap in self._haps]
                yield get_result(poss, all_hap_bases, all_typed_bases)

        def iter_het_blocks2():
            vb = self._hmm_calc.get_viterbi_blocks()
            offset = 0
            for block in vb['blocks']:
                path = block['path']
                all_typed_bases = [[] for _ in range(self._copy_number)]
                poss = []
                for state_idx, d in zip(path, self._hetd[offset:offset + len(path)]):
                    poss.append(d['pos'])
                    uniq_bases = d['uniq_bases']
                    state = d['states'][state_idx]
                    for j, s in enumerate(state):
                        all_typed_bases[j].append(uniq_bases[s])
                all_hap_bases = [[hap[pos] for pos in poss] for hap in self._haps]
                yield get_result(poss, all_hap_bases, all_typed_bases)
                offset += len(path)

        matches = []
        matches.append(get_hom_block())
        if self._hmm_calc:
            #matches.append(get_het_block())
            #matches.extend(iter_het_blocks())  # does not work well
            matches.extend(iter_het_blocks2())  # does not work well
        return matches


class VarCallInfo(object):
    def __init__(self, vc_info):
        self._tab = pd.read_table(vc_info)
        tab = pd.read_table(vc_info)
        self._refs = {}
        self.refnames = []
        for refname, tab1 in tab.groupby('refname'):
            self.refnames.append(refname)
            vc = RefVarCall(tab1)
            self._refs[refname] = vc

    def get_var_call(self, refname):
        return self._refs.get(refname)

    def rank_haps(self, refname, contigs, threshold=5, score=None):
        """
        Args:
            contigs : [FastaContig]
        Returns:
            {
                'order': [hap_idxs]
                'score': [score]
                'match_info': [match_info]
            }
        """
        copy_number = self._refs[refname].copy_number
        if copy_number == 0:
            return {'order': [], 'score': [], 'match_info': [], 'copy_number': 0}

        #if self._refs[refname].copy_number == 1:
        #    assert self._refs[refname].copy_number == 1
        #    return {}
        #assert self._refs[refname].copy_number == 2  # TODO polyploid

        mask_seqs = [(c.name, c.get_clip_mask_seq(mask_char='N')) for c in contigs]  # [(name, seq)]
        vc = self._refs[refname]
        if score is None:
            score = 'edit'
        def score_fn(hap_idxs):  # TODO optimize
            match_info = vc.match_haps([mask_seqs[i][1] for i in hap_idxs], score=score)
            return match_info.score

        if score == 'hmm':
            threshold = threshold * - np.log(MISMATCH_PROB)   # adjust score

        # 'A*01:01:01:01', 'A*11:01:01:01'
        # test_idxs = [idx for idx, c in enumerate(contigs) if c.name in ('HLA:HLA00001', 'HLA:HLA00043')]
        # logging.info((test_idxs, score_fn(test_idxs)))
        # test_idxs = [idx for idx, c in enumerate(contigs) if c.name in ('HLA:HLA14798',)]
        # logging.info((test_idxs, score_fn(test_idxs)))
        # test_idxs = [idx for idx, c in enumerate(contigs) if c.name in ('HLA:HLA00001',)]
        # logging.info((test_idxs, score_fn(test_idxs)))
        # test_idxs = [idx for idx, c in enumerate(contigs) if c.name in ('HLA:HLA00043',)]  # actually low score
        # logging.info((test_idxs, score_fn(test_idxs)))

        if copy_number == 1:
            result = rank_combinations1(len(mask_seqs), score_fn=score_fn, threshold=threshold)
        elif copy_number == 2:
            result = rank_combinations2(len(mask_seqs), score_fn=score_fn, threshold=threshold)
        else:
            raise NotImplementedError()
            result = rank_combinations(copy_number, len(mask_seqs), score_fn=score_fn, threshold=threshold)

        match_infos = []
        for hap_idxs in result['order']:
            #logging.info(hap_idxs)
            match_info = vc.match_haps([mask_seqs[i][1] for i in hap_idxs], score=score)
            #logging.info(idxs)
            #logging.info(match_info)
            match_infos.append(match_info)
        result['copy_number'] = vc.copy_number
        result['match_info'] = match_infos
        return result


def rank_combinations1(n, score_fn, threshold=0):
    L = -float('inf')
    scores1 = [score_fn((j,)) for j in range(n)]
    O1 = [j for j, score in sorted(enumerate(scores1), key=lambda x: x[1], reverse=True)]
    S1 = [scores1[j] for j in O1]
    ranks = []
    rank = 0
    last = float('inf')
    for s in S1:
        if s < last:
            rank += 1
        ranks.append(rank)
        last = s

    order = [(j,) for j in O1]
    return {'order': order, 'score': S1, 'rank': ranks}

# for diploid
def rank_combinations2(n, score_fn, threshold=0):
    """
    Assumptions for score_fn : S(els)
    # definition
      els = [el]
      el in {0, ..., n-1}
      0 <= len(els) <= dim
      S(permute(els)) == S(els) # exchangable
      # S(els1 v els2) + S(els1 ^ els2) <= S(els1) + S(els2) # submodular ?
      S(els1 v els2) <= S(els1) + S(els2)

    >>> values = {
    ...     (0,): 5,
    ...     (1,): 2,
    ...     (2,): 4,
    ...     (3,): 5,
    ...     (0, 0): 8, (0, 1): 7, (0, 2): 9, (0, 3): 10,
    ...     (1, 1): 3, (1, 2): 5, (1, 3): 7,
    ...     (2, 2): 8, (2, 3): 9,
    ...     (3, 3): 8,
    ... }
    >>> assert all(values[i, j] <= values[i,] + values[j,] for i in range(4) for j in range(i, 4))  # validation of values
    >>> result = rank_combinations2(4, score_fn=(lambda v: values[tuple(sorted(v))]))
    >>> result['score']
    [10]
    >>> result['order']
    [(0, 3)]
    >>> result['rank']
    [1]

    >>> result = rank_combinations2(4, score_fn=(lambda v: values[tuple(sorted(v))]), threshold=1)
    >>> result['score']
    [10, 9, 9]
    >>> result['order']
    [(0, 3), (0, 2), (3, 2)]
    >>> result['rank']
    [1, 2, 2]

    >>> result = rank_combinations2(4, score_fn=(lambda v: values[tuple(sorted(v))]), threshold=2)
    >>> result['score']
    [10, 9, 9, 8, 8, 8]
    >>> result['order']
    [(0, 3), (0, 2), (3, 2), (0, 0), (3, 3), (2, 2)]
    >>> result['rank']
    [1, 2, 2, 3, 3, 3]
    """
    L = -float('inf')
    scores1 = [score_fn((j,)) for j in range(n)]
    O1 = [j for j, score in sorted(enumerate(scores1), key=lambda x: x[1], reverse=True)]
    S1 = [scores1[j] for j in O1]
    U1 = []  # init
    for i in range(n):
        U1.append(S1[i])
    O2 = []
    S2 = []

    for i1 in range(n):
        j1 = O1[i1]
        if S1[i1] + U1[i1] < L:  # pruning
            continue
        for i2 in range(i1, n):
            if S1[i1] + S1[i2] < L: # pruning
                break
            j2 = O1[i2]
            score = score_fn((j1, j2))
            if score >= L:
                O2.append((j1, j2))
                S2.append(score)
                os2 = list(sorted(zip(O2, S2), key=lambda x: x[1], reverse=True))  # TODO optimize
                O2 = [o for o, s in os2]
                S2 = [s for o, s in os2]
                # update lower bound
                L = S2[0] - threshold
                while 1:
                    s = S2[-1]
                    if s < L:
                        O2.pop()
                        S2.pop()
                    else:
                        break

    ranks = []
    rank = 0
    last = float('inf')
    for s in S2:
        if s < last:
            rank += 1
        ranks.append(rank)
        last = s

    return {'order': O2, 'score': S2, 'rank': ranks}


def rank_combinations(n, dim, score_fn, threshold=0):
    """
    Assumptions for score_fn : S(els)
    # definition
        els = [el]
        el in {0, ..., n-1}
        0 <= len(els) <= dim
        S([]) = 0
        S(permute(els)) == S(els) # exchangable
        # S(els1 v els2) + S(els1 ^ els2) <= S(els1) + S(els2) # submodular ?
        S(els1 v els2) <= S(els1) + S(els2)
    """
    L = -float('inf')
    scores1 = [score_fn(i) for i in range(n)]
    ordered = [i for i, score in sorted(enumerate(scores1), key=lambda x: x[1], reverse=True)]
    S = []
    S1 = [scores1[i] for i in ordered]
    S.append(S1)
    U = []
    U1 = []  # init
    for i in ordered:
        U1.append(S[0][i])
    U.append(U1)
    # TODO


class HapMatch(object):
    def __init__(self, vc_info, ref_db):
        self._vc_info = vc_info
        self._ref_db = ref_db

    @property
    def refs(self):
        return args._vc_info.refs

    def refnames(self):
        return self._vc_info.refnames

    def get_bests(self, refname, threshold=5, score=None):
        """
            result: {
                'order': [indexes]
                'score': [score]
                'rank': [rank]
                'match_info': [match_info]
                'hrecs': [{
                    'score': sum of scores
                    'edit': sum of edits
                    'match': sum of matches
                    'blocks': [{
                        'block_id': block_id,
                        'indexes': (index,),  # selected order e.g. (0,) or (1,) for one hap, (0, 1) or (1, 0) for two haps
                        'match': match,  # sum of matches
                        'edit': edit, # sum of edits
                        'matches': [match], # match for each (phased, hap)
                        'edits': [edit], # edit for each (phased, hap)
                        'edit_details': [[(pos, ref_base, alt_base)]],   # each edit
                        'score': # total score,
                        'score_mode': # total score,
                    }]
                }]
            }
        """
        contigs = self._ref_db.get_contigs(refname)   # [(name, seq)]
        logging.info('number of contigs: %s', len(contigs))
        res = self._vc_info.rank_haps(refname, contigs, threshold=threshold, score=score)
        def get_hrec(blocks, j, hap_name):
            match = edit = 0
            edit_details = []
            for block in blocks:
                #k = block['indexes'].index(j)
                #k = block['indexes'][j]
                #match += block['matches'][k]
                #edit += block['edits'][k]
                #edit_details.extend(block['edit_details'][k])
                match += block['matches'][j]
                edit += block['edits'][j]
                edit_details.extend(block['edit_details'][j])

            return {
                'name': hap_name,
                'match': match,
                'edit': edit,
                'edit_details': edit_details,
            }

        cn = res['copy_number']
        recs = []
        if cn == 0:
            return recs
        for i, hap_idxs in enumerate(res['order']):
            hap_names = [contigs[idx].name for idx in hap_idxs]
            match_info = res['match_info'][i]
            block_matches = match_info.get_block_matches()
            hrecs = [get_hrec(block_matches, j, hap_names[j]) for j in range(cn)]
            rec = {
                    'refname': refname,
                    'ref_cn': cn,
                    'score': match_info.score,
                    'edit': match_info.edit,
                    'match': match_info.match,
                    'rank': res['rank'][i],
                    'haps': hrecs,
                    'nblocks': len(block_matches),
                    'block_ids': [m['block_id'] for m in block_matches],
            }
            recs.append(rec)
        return recs

    def refine_variant_table(self, refname):
        contigs = self._ref_db.get_contigs(refname)   # [(name, seq)]
        logging.info('number of contigs: %s', len(contigs))
        res = self._vc_info.rank_haps(refname, contigs, threshold=0, score='hmm')
        hap_idxs = res['order'][0]
        block_matches = res['match_info'][0].get_block_matches()

    def get_bests(self, refname):
        haps = self._ref_db.get_haps(refname)
        self._vc_info.get_score(refname, haps)
