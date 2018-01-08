from __future__ import print_function
import logging
from builtins import range
import numpy as np
from ..utils import cached_property


def logsumexp(vals, min_prob=0, axis=None):
    """
    logsumexp(vals)
    = log(sum(exp(vals - m) * exp(m)))
    = log(sum(exp(vals - m)) * exp(m))
    = log(sum(exp(vals - m))) + m

    >>> np.isclose(logsumexp(np.array([0, 0, 0])), np.log(np.sum(np.exp(np.array([0, 0, 0])))))
    True
    >>> np.isclose(logsumexp(np.array([100, 50, 20])), np.log(np.sum(np.exp(np.array([100, 50, 20])))))
    True
    >>> np.isclose(logsumexp(np.array([[1, 2], [3, 4]])), np.log(np.sum(np.exp(np.array([[1, 2], [3, 4]])))))
    True
    >>> np.isclose(logsumexp(np.array([[1, 2], [3, 4]]), axis=1), np.log(np.sum(np.exp(np.array([[1, 2], [3, 4]])), axis=1))).all()
    True
    >>> val = [[[1, 2], [3, 4]], [[5, 6], [7, 8]]]
    >>> np.isclose(logsumexp(np.array(val)), np.log(np.sum(np.exp(np.array(val)))))
    True
    >>> np.isclose(logsumexp(np.array(val), axis=0), np.log(np.sum(np.exp(np.array(val)), axis=0))).all()
    True
    >>> np.isclose(logsumexp(np.array(val), axis=1), np.log(np.sum(np.exp(np.array(val)), axis=1))).all()
    True
    >>> np.isclose(logsumexp(np.array(val), axis=2), np.log(np.sum(np.exp(np.array(val)), axis=2))).all()
    True
    """
    m = vals.max(axis=axis)
    if axis is not None:
        diffs = vals - np.expand_dims(m, axis=axis)
    else:
        diffs = vals - m
    return np.log(np.maximum(np.exp(diffs), min_prob).sum(axis=axis)) + m


def lnormal(val, axis=None):
    """
    >>> na = np.array
    >>> np.isclose(1, np.sum(np.exp(lnormal(na([0, 0, 0])))))
    True
    >>> list(np.isclose(1, np.sum(np.exp(lnormal(na([[1, 2, 1], [2, 1, 3]]), axis=1)), axis=1)))
    [True, True]
    """
    t = logsumexp(val, axis=axis)
    if axis is None:
        return val - t
    else:
        return val - np.expand_dims(t, axis=axis)


class HMMCalc(object):
    """
    >>> dims = [2, 3]
    >>> na = np.array
    >>> linit = lnormal(na([0, 1]))
    >>> ltrans = [lnormal(na([[0, 0, 0], [1, 0, 0]]), axis=1)]
    >>> lemit = [lnormal(na([0, 0])), lnormal(na([0, 0, 0]))]
    >>> h = HMMCalc(dims, linit, ltrans, lemit)
    >>> h.best_states
    [1, 0]
    >>> ltrans = [lnormal(na([[0, 0, 0], [-1, 1, 2]]), axis=1)]
    >>> lemit = [lnormal(na([0, 0])), lnormal(na([0, 0, 0]))]
    >>> h = HMMCalc(dims, linit, ltrans, lemit)
    >>> h.best_states
    [1, 2]

    # State rotation
    >>> linit1 = lnormal(na([0, 1]))
    >>> ltrans1 = [lnormal(na([[0, 0, 0], [1, 0, 0]]), axis=1)]
    >>> lemit1 = [lnormal(na([0, 1])), lnormal(na([1, 0, 2]))]
    >>> h1 = HMMCalc(dims, linit1, ltrans1, lemit1)
    >>> linit2 = lnormal(na([1, 0]))
    >>> ltrans2 = [lnormal(na([[0, 0, 1], [0, 0, 0]]), axis=1)]
    >>> lemit2 = [lnormal(na([1, 0])), lnormal(na([0, 2, 1]))]
    >>> h2 = HMMCalc(dims, linit2, ltrans2, lemit2)
    >>> h1.loglikelihood == h2.loglikelihood
    True
    >>> h.loglikelihood != h1.loglikelihood
    True
    """
    min_prob = 1e-100

    def __init__(self, dims, linit, ltrans, lemit):
        """
        dims: dimensions of states
        linit: log of initail state probabilities  (len(dims[0]))
        ltrans: list of log of transition probabilities  (len(dims) - 1)
        lemit: list of log of emission probabilities  (len(dims))
        """
        self.length = len(dims)
        self._dims = dims # number of states for each state
        self._linit = linit
        self._ltrans = ltrans  # first element
        self._lemit = lemit

    @cached_property
    def loglikelihood(self):
        M = len(self._dims)
        la = logsumexp(self._linit + self._lemit[0], min_prob=self.min_prob)  # first element
        ll = logsumexp(self._lbs[0] + la, min_prob=self.min_prob)
        return ll

    @cached_property
    def _las(self):
        M = len(self._dims)
        las = [] # forward msgs : alpha_ij = P(y[1:i], x[i] = j)
        la = self._linit + self._lemit[0]
        las.append(la)
        for i in range(1, M):
            dim = self._dims[i]
            term = logsumexp(la[:, None] + self._ltrans[i-1], axis=0, min_prob=self.min_prob)
            la = term + self._lemit[i]
            las.append(la)

        return las

    @cached_property
    def _lbs(self):
        M = len(self._dims)
        rev_lbs = [] # backward msgs : beta_ij = P(y[i:M] | x[i] = j)
        lb = self._lemit[M - 1]
        rev_lbs.append(lb)
        for i in range(M-1 - 1, -1, -1):
            dim = self._dims[i]
            term = logsumexp(lb[None, :] + self._ltrans[i], axis=1, min_prob=self.min_prob)
            lb = term + self._lemit[i]
            rev_lbs.append(lb)

        return list(reversed(rev_lbs))

    @cached_property
    def viterbi(self):
        """
        path = argmax_x P(x | y)
        Returns:
            'path': [state_index]   # dim: M
            'state_probs': [prob]   # dim: M
            'gammas': []    # gamma_{i+1}(x_{i+1}) ~ P(y_{i+1} | x_{i+1}) max_{x_i} [ gamma_i(x_i) P(x_{i+1} | x_i) ]
            'pointers': []  # pointers for backtrack

        gamma should be normalized
        """
        M = len(self._dims)
        gammas = []
        pointers = []   # pointers from (i+1)th to ith node
        gamma = self._linit + self._lemit[0]   # this is already normalized
        gammas.append(gamma)
        for i in range(1, M):
            dim = self._dims[i]
            term = gamma[:, None] + self._ltrans[i-1]
            ps = np.argmax(term, axis=0)
            gamma = term.max(axis=0) + self._lemit[i]
            gamma -= logsumexp(gamma, min_prob=self.min_prob) # normalize
            pointers.append(ps)
            gammas.append(gamma)
        assert len(gammas) == M
        assert len(pointers) == M - 1

        # backtracking
        rev_path = []
        sidx = np.argmax(gammas[-1])
        rev_path.append(sidx)
        for ps in pointers:
            sidx = ps[sidx]
            rev_path.append(sidx)
        path = list(reversed(rev_path))
        assert len(path) == M

        state_probs = [np.exp(gammas[i][sidx]) for i, sidx in enumerate(path)]
        return {'path': path, 'state_probs': state_probs, 'gammas': gammas, 'pointers': pointers}

    @cached_property
    def viterbi_marginal(self):
        """
        path = argmax_x P(x1, .., xl | y)
        Returns:
            'path': [state_index]   # dim: M
            'state_probs': [prob]   # dim: M
            'deltas': []    # gamma_{i+1}(x_{i+1}) ~ P(y_{i+1} | x_{i+1}) max_{x_i} [ gamma_i(x_i) P(x_{i+1} | x_i) ]
            'pointers': []  # pointers for backtrack

        gamma should be normalized
        """
        M = len(self._dims)
        v = self.viterbi
        deltas = []
        pointers = []
        delta = self._linit + self._lbs[0]
        deltas.append(delta)
        for i in range(1, M):
            dim = self._dims[i]
            gamma = v['gammas'][i-1]
            term = (self._ltrans[i-1] + gamma[:, None])
            delta = term.max(axis=0) + self._lbs[i]
            delta -= logsumexp(delta, min_prob=self.min_prob) # normalize
            deltas.append(delta)
            ps = np.argmax(term, axis=0)
            pointers.append(ps)
        assert len(deltas) == M
        assert len(pointers) == M - 1

        # backtracking
        rev_path = []
        sidx = np.argmax(deltas[-1])
        rev_path.append(sidx)
        for ps in pointers:
            sidx = ps[sidx]
            rev_path.append(sidx)
        path = list(reversed(rev_path))
        assert len(path) == M

        state_probs = [np.exp(deltas[i][sidx]) for i, sidx in enumerate(path)]
        return {'path': path, 'state_probs': state_probs, 'deltas': deltas, 'pointers': pointers}

    def get_viterbi_blocks(self, threshold=.9):
        """
        path = argmax_x P(x | y)
        Returns: {
                'blocks': [{
                    'path': [state_index]   # dim: M
                    'state_probs': [prob]   # dim: M
                }]
                'gammas': []    # gamma_{i+1}(x_{i+1}) ~ P(y_{i+1} | x_{i+1}) max_{x_i} [ gamma_i(x_i) P(x_{i+1} | x_i) ]
                'pointers': []  # pointers for backtrack
            }

        gamma should be normalized
        """
        gammas = self.viterbi['gammas']
        pointers = self.viterbi['pointers']
        M = len(self._dims)
        assert len(gammas) == M
        assert len(pointers) == M - 1
        path_revs = []
        # backtracking
        it = enumerate(reversed(pointers + [None]), 1)
        try:
            rev_path = []
            while 1:
                i, ps = next(it)
                egamma = np.exp(gammas[-i])
                if ps is None:
                    sidx = np.argmax(egamma)
                    rev_path.append(sidx)
                else:
                    sidx = ps[sidx]
                    if egamma[sidx] < threshold:
                        path = list(reversed(rev_path))
                        path_revs.append(path)
                        sidx = np.argmax(egamma)  # replace sidx
                        rev_path = [sidx]
                        continue
                    else:
                        rev_path.append(sidx)
        except StopIteration:
            if rev_path:
                path = list(reversed(rev_path))
                path_revs.append(path)

        blocks = []
        offset = 0
        for path in reversed(path_revs):
            state_probs = [np.exp(gammas[i][sidx]) for i, sidx in enumerate(path, offset)]
            blocks.append({
                'path': path,
                'state_probs': state_probs,
            })
            offset += len(path)
            #logging.info(path)

        #logging.info((sum(len(b['path']) for b in blocks), M))
        assert sum(len(b['path']) for b in blocks) == M
        return {'blocks': blocks, 'gammas': gammas, 'pointers': pointers}

    @cached_property
    def state_lls(self):   # log likelihood of each state
        M = len(self._dims)
        lls = []
        ll = self._linit + self._lbs[0]
        lls.append(ll)
        for i in range(1, M):
            dim = self._dims[i]
            #ll = np.array([logsumexp(self._las[i-1] + self._ltrans[i-1][:, k], min_prob=self.min_prob) for k in range(dim)]) + self._lbs[i]
            ll = logsumexp(self._las[i-1][:, None] + self._ltrans[i-1], axis=0, min_prob=self.min_prob) + self._lbs[i]
            lls.append(ll)
            assert ll.shape == (self._dims[i],), (i, M, self._las[i-1], self._ltrans[i-1], self._lbs[i])
        return lls

    @cached_property
    def best_states(self):   # log likelihood of each state
        """
        Returns: e.g. [0, 1, 1, ..]
        """
        assert self._dims == [len(lls) for lls in self.state_lls]
        return [np.argmax(lls) for lls in self.state_lls]

    def get_best_path(self):  # TODO blocked path with threshold
        """
        Returns: [states]
        """  # TODO
        return

    def get_opt_blocks(self, threshold=.9):  # TODO blocked path with threshold
        """
        Returns: [{'block': [], 'states': [],  'prob': }]
        """
        blocks = []
        M = len(self._dims)
        rev_lbs = []
        for i in range(1, M):
            dim = self._dims[i]
            term = logsumexp(lb[None, :] + self._ltrans[i], axis=1, min_prob=self.min_prob)
            lb = term + self._lemit[i]
            rev_lbs.append(lb)

        return blocks

