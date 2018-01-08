from __future__ import print_function
import logging
import numpy as np
from collections import namedtuple
from ..utils import Counter
from scipy.stats import norm
#from depth_norm import em
cimport numpy as np
ctypedef np.float64_t DTYPE_t

Result = namedtuple('Result', 'a b theta step ll')

cdef DTYPE_t LOG_DUMMY_ZERO = -60
cdef DTYPE_t DUMMY_ZERO = np.exp(LOG_DUMMY_ZERO)

cdef class Normalizer(object):
    cdef:
        public int max_copy, step, step_limit
        DTYPE_t a, b, la, lb, max_diff, x_sum, mu_sum
        int Z_dim, Z2_dim, N, total_ploidy
        np.ndarray x, mu, theta, ploidies, Z, Z2, z1_mean, z2_mean
        np.ndarray theta_mat, Z_sum_mat, z1_probs, z2_probs # 2 dim
        str name
        list z_probs

    def __init__(self, counts, exp_counts, ploidies, max_copy=3, a=0., b=1., la=1., lb=1., theta=None, max_diff=1e-6, step_limit=0, name='default'):
        self.x = counts
        self.mu = exp_counts
        self.x_sum = self.x.sum()
        self.mu_sum = self.mu.sum()
        self.ploidies = np.array(ploidies)
        self.total_ploidy = self.ploidies.sum()
        assert self.total_ploidy > 0, 'total ploidy should be > 0'
        self.max_copy = max_copy
        self.a = a
        self.b = b
        self.la = la
        self.lb = lb
        self.max_diff = max_diff
        self.step_limit = step_limit
        self.step = 0
        self.Z_dim = max_copy + 1
        self.Z2_dim = max_copy * 2 + 1
        self.Z = np.arange(self.Z_dim)
        self.Z2 = np.arange(self.Z2_dim)
        self.N = len(counts)
        self.theta = theta
        self.name = name
        self.z_probs = []

        if self.theta is None:
            self.theta = 1. * np.zeros(self.Z_dim)
            # setting initial theta values
            #logging.info(self.x)
            #logging.info(self.mu)
            #logging.info(self.x / self.mu)
            tmp_hap_cns = [int(k) for k in np.round((self.x - a) / self.mu / b / 2.)]
            counter = Counter(tmp_hap_cns)
            less = sum([counter[k] for k in counter.keys() if k < 0])
            more = sum([counter[k] for k in counter.keys() if k > max_copy])
            counter_trim = Counter(dict((k, counter[k]) for k in xrange(max_copy + 1)))
            counter_trim[0] += less
            counter_trim[max_copy] += more
            theta0 = 1. # smoothing
            self.theta = theta = 1. * np.array([counter_trim[k] + theta0 for k in xrange(max_copy + 1)]) / (len(self.x) + theta0 * (max_copy + 1))  # initial value for theta
        logging.info('[%s] Initial values: a=%s, b=%s, theta=%s', name, a, b, theta)

        self.theta_mat = self.theta[:, None] * self.theta[None, :]
        self.Z_sum_mat = np.array(self.Z[:, None]) + np.array(self.Z[None, :])

    cdef DTYPE_t _get_next_a(self):
        cdef DTYPE_t a = self.a
        cdef DTYPE_t a1 = float('inf')
        cdef int loop = 1
        cdef DTYPE_t mu_sum_b = self.mu_sum * self.b
        #cdef DTYPE_t x0_term = ((self.x - 1) * np.array([z_prob[0, 0] if p == 2 else z_prob[0] for z_prob, p in zip(self.z_probs, self.ploidies)])).sum()
        cdef DTYPE_t x0_term = (self.x * self.z2_probs[:, 0]).sum()
        cdef DTYPE_t h_term

        while abs(a1 - a) > 1e-3:
            a1 = a  # prev
            #h_term = np.array([(a * z_prob / (DUMMY_ZERO + self.Z_sum_mat + a)).sum() if p == 2 else (a * self.Z / (DUMMY_ZERO + z_prob + a)).sum() for z_prob, p in zip(self.z_probs, self.ploidies)]).sum()
            h_term = (self.x * a * (self.z2_probs[:, 1:] / (self.Z2[None, 1:] + a)).sum(axis=1)).sum()
            #logging.info('a: %s -> %s, mu_sum_b: %s, la: %s, h_term: %s', a1, a, mu_sum_b, self.la, h_term)
            a = (x0_term + h_term) / (self.la + mu_sum_b)
            loop += 1
            if loop % 100 == 0:
                logging.warning('[%s] loop: %s, a: %s -> %s', self.name, loop, a1, a)
        #logging.info('[%s] loop: %s, a: %s -> %s, (%s)', self.name, loop, a1, a, (x0_term, h_term, self.la, mu_sum_b))
        return a

    cdef DTYPE_t _get_next_b(self):
        cdef DTYPE_t A = ((self.a + self.z2_mean) * self.mu).sum()
        cdef DTYPE_t b = self.b
        cdef DTYPE_t log_b1 = float('inf')
        cdef DTYPE_t log_b = np.log(b)
        cdef int loop = 1

        while abs(log_b1 - log_b) > 1e-3:
            log_b1 = log_b  # prev
            log_b = log_b + (self.x_sum - 1 - A * b - self.lb * log_b) / (A * b + self.lb)
            b = np.exp(log_b)
            loop += 1
            if loop % 100 == 0:
                logging.warning('[%s] loop: %s, log_b: %s -> %s', self.name, loop, log_b1, log_b)
        return b

    def _iter_z_probs(self):
        cdef DTYPE_t x1, mu1
        cdef int ploidy
        cdef np.ndarray likes, log_likes1, y1
        cdef int i, j
        cdef np.ndarray probs2

        for x1, mu1, ploidy in zip(self.x, self.mu, self.ploidies):
            if ploidy == 0:
                yield (np.array(0),
                        (0,) * self.Z_dim,
                        (0,) * self.Z2_dim)
            elif ploidy == 1:
                y1 = (self.a + self.Z) * mu1 * self.b
                log_likes = x1 * np.log(y1 + DUMMY_ZERO) - y1 + np.log(self.theta + DUMMY_ZERO)
                log_likes -= np.max(log_likes)
                #log_likes = np.where(np.isnan(log_likes) | (log_likes < LOG_DUMMY_ZERO), LOG_DUMMY_ZERO, log_likes)
                likes = np.exp(log_likes)
                likes = np.where(np.isnan(likes), 0, likes)
                prob = likes / likes.sum()
                assert not np.any(np.isnan(prob)), 'prob should not be nan'
                yield (prob,
                       tuple(prob),
                       tuple(prob) + (0,) * (self.Z2_dim - self.Z_dim))
            elif ploidy == 2:
                y1 = (self.a + self.Z_sum_mat) * mu1 * self.b
                log_likes = x1 * np.log(y1 + DUMMY_ZERO) - y1 + np.log(self.theta_mat + DUMMY_ZERO)
                log_likes -= np.max(log_likes)
                #log_likes = np.where(np.isnan(log_likes) | (log_likes < LOG_DUMMY_ZERO), LOG_DUMMY_ZERO, log_likes)
                likes = np.exp(log_likes)
                likes = np.where(np.isnan(likes), 0, likes)
                #assert not np.any(np.isnan(likes))
                prob = likes / likes.sum()
                assert not np.any(np.isnan(prob)), 'prob should not be nan'
                probs2 = np.zeros(self.Z2_dim)
                for i in xrange(self.Z_dim):
                    for j in xrange(self.Z_dim):
                        probs2[i + j] += prob[i, j]
                yield (prob,
                       tuple(prob.sum(axis=1)),
                       tuple(probs2))

    def iter_genotypes(self):
        dip_copy = self.max_copy * 2

        for p, (zp, z1p) in zip(self.ploidies, self._iter_z_probs()):
            if p == 0:
                yield (1,) + (0,) * dip_copy
            elif p == 1:
                remain = dip_copy + 1 - len(zp)
                yield tuple(zp) + (0,) * remain
            else:
                dip_probs = [0] * (dip_copy + 1)
                for i in xrange(self.Z_dim):
                    for j in xrange(self.Z_dim):
                        dip_probs[i + j] += zp[i, j]
                yield tuple(dip_probs)

    cdef DTYPE_t _calc_ll(self):
        log_b = np.log(self.b)
        ll_next = (self.x[:, None] * ((np.log(self.a + self.Z2 + DUMMY_ZERO)[None, :] * self.z2_probs).sum(axis=1) + log_b)).sum()
        ll_next += - ((self.a + self.z2_mean) * self.mu).sum() * self.b
        ll_next += (self.ploidies[:, None] * np.log(self.theta + DUMMY_ZERO)[None, :] * self.z1_probs).sum()
        ll_next += - np.sum([(z_prob * np.log(z_prob + DUMMY_ZERO)).sum() for z_prob in self.z_probs])   # entropy
        ll_next += - self.la * self.a
        ll_next += - log_b - (self.lb / 2.) * log_b ** 2
        return ll_next

    def iteration(self):
        cdef list _z1_probs, _z2_probs
        cdef DTYPE_t p
        cdef DTYPE_t ll_next, term1
        cdef DTYPE_t diff = 1
        cdef DTYPE_t ll = -float('inf')
        cdef DTYPE_t log_b
        cdef np.ndarray z_prob
        cdef np.ndarray[long, ndim=2] Z_sum_mat
        cdef np.ndarray zp
        cdef tuple z1p, z2p
        cdef int ll_check_iv = 100
        cdef DTYPE_t prev_a = self.a
        cdef DTYPE_t prev_b = self.b
        cdef np.ndarray[DTYPE_t] prev_theta = self.theta
        Z_sum_mat = self.Z_sum_mat

        while diff > self.max_diff:
            if self.step_limit and self.step > self.step_limit:
                raise AssertionError('The number of steps exceeds the limit')
            self.step += 1

            self.z_probs = []
            _z1_probs = []
            _z2_probs = []
            for zp, z1p, z2p in self._iter_z_probs():
                self.z_probs.append(zp)
                _z1_probs.append(z1p)
                _z2_probs.append(z2p)
            self.z1_probs = np.array(_z1_probs)  # {x, Z_dim}
            self.z2_probs = np.array(_z2_probs)  # {x, Z2_dim}
            self.z1_mean = (self.Z[None, :] * self.z1_probs).sum(axis=1)
            self.z2_mean = (self.Z2[None, :] * self.z2_probs).sum(axis=1)

            self.theta = np.array([list(z_prob.sum(axis=0) + z_prob.sum(axis=1)) if p == 2 else list(z_prob)
                    for z_prob, p in zip(self.z_probs, self.ploidies) if p > 0]).sum(axis=0) / self.total_ploidy
            self.theta_mat = self.theta[:, None] * self.theta[None, :]
            self.a = self._get_next_a()
            self.b = self._get_next_b()
            assert not np.isnan(self.a), 'a is nan'
            assert not np.isnan(self.b), 'b is nan'

            diff = max(np.abs(self.a - prev_a), np.abs(self.b - prev_b), np.max(np.abs(self.theta - prev_theta)))
            prev_a = self.a
            prev_b = self.b
            prev_theta = self.theta.copy()
            #d = (self.x - self.a) / self.b
            if self.step % ll_check_iv == 0:
                ll_next = self._calc_ll()
                ll_diff = ll_next - ll
                #logging.info(ll_next)
                #assert ll_diff >= -1e-2, 'The lower bound of the log-likelihood should be increased (ll_diff: {0}, diff: {1})'.format(ll_diff, diff)
                logging.warning('The lower bound of the log-likelihood decreased (ll_diff: {0}, diff: {1})'.format(ll_diff, diff))
                ll = ll_next
                logging.info('[%s] ll: %s (%s), %s', self.name, ll, diff, Result(self.a, self.b, self.theta, self.step, ll))
            # y = a + b * z2

            ll = self._calc_ll()
        return Result(self.a, self.b, self.theta, self.step, ll)


def normalize_em(*args, **kwds):
    normalizer = Normalizer(*args, **kwds)
    return normalizer.iteration()
