#!/usr/bin/env python
from __future__ import print_function
import logging
import numpy as np
from collections import namedtuple
from ..utils import Counter
from scipy.stats import norm
#from depth_norm import em
cimport numpy as np
ctypedef np.float64_t DTYPE_t

Result = namedtuple('Result', 'a b c d theta step ll')

cdef class Normalizer(object):
    cdef:
        #np.float_t[:] x, theta
        DTYPE_t a, b, c, la, lb, lc, max_diff, alpha, omega
        int max_copy, step_limit, Z_dim, N, total_ploidy
        np.ndarray x, d_mean, d2_mean, x_a_d_mean, theta, ploidies, Z, Z2
        np.ndarray theta_mat, Z_sum_mat, Z2_sum_mat   # 2 dim
        str name

    def __init__(self, x, ploidies, max_copy=3, a=0., b=1., c=10., la=0., lb=.3, lc=.3, theta=None, max_diff=1e-6, step_limit=0, name='default'):
        self.x = x
        self.ploidies = np.array(ploidies)
        self.total_ploidy = self.ploidies.sum()
        assert self.total_ploidy > 0, 'total ploidy should be > 0'
        self.max_copy = max_copy
        self.a = a
        self.b = b
        self.c = c
        self.la = la
        self.lb = lb
        self.lc = lc
        self.max_diff = max_diff
        self.step_limit = step_limit
        self.Z_dim = max_copy + 1
        self.Z = np.arange(self.Z_dim)
        self.N = len(x)
        self.theta = theta
        self.name = name

        if self.theta is None:
            self.theta = 1. * np.zeros(self.Z_dim)
            # setting initial theta values
            tmp_hap_cns = [int(k) for k in np.round((x - a) / b / 2.)]
            counter = Counter(tmp_hap_cns)
            less = sum([counter[k] for k in counter.keys() if k < 0])
            more = sum([counter[k] for k in counter.keys() if k > max_copy])
            counter_trim = Counter(dict((k, counter[k]) for k in xrange(max_copy + 1)))
            counter_trim[0] += less
            counter_trim[max_copy] += more
            theta0 = 1. # smoothing
            self.theta = theta = 1. * np.array([counter_trim[k] + theta0 for k in xrange(max_copy + 1)]) / (len(x) + theta0 * (max_copy + 1))  # initial value for theta
        logging.info('[%s] Initial values: a=%s, b=%s, c=%s, theta=%s', name, a, b, c, theta)

        self.theta_mat = self.theta[:, None] * self.theta[None, :]
        self.Z_sum_mat = np.array(self.Z[:, None]) + np.array(self.Z[None, :])
        self.Z2 = np.array(self.Z) ** 2
        self.Z2_sum_mat = np.array(self.Z_sum_mat) ** 2

        #alpha = 20 # skew parameter of a
        self.alpha = 1. # skew parameter of a
        #omega = 0.2   # learning rate of a
        self.omega = 1.   # learning rate of a

    cdef DTYPE_t _get_next_a(self):
        cdef DTYPE_t a = self.a
        cdef DTYPE_t a1 = float('inf')
        cdef DTYPE_t rate
        cdef int loop = 1

        while abs(a1 - a) > 1e-3:
            rate = norm.pdf(self.alpha * a) / norm.cdf(self.alpha * a)
            a1 = a  # prev
            a = a * (1. - self.omega) + self.omega * (self.c * (self.x.sum() - self.b * self.d_mean.sum()) + self.alpha * rate) / (self.c * self.N + self.la)
            loop += 1
            if loop % 100 == 0:
                logging.warning('[%s] loop: %s, a: %s -> %s', self.name, loop, a1, a)
        return a

    cdef DTYPE_t _get_next_b(self):
        cdef DTYPE_t A = self.x_a_d_mean.sum()
        cdef DTYPE_t B = self.d2_mean.sum()
        cdef DTYPE_t C = self.lb/self.c
        cdef DTYPE_t b = self.b
        cdef DTYPE_t log_b1 = float('inf')
        cdef DTYPE_t log_b = np.log(b)
        cdef int loop = 1
        while abs(log_b1 - log_b) > 1e-3:
            log_b1 = log_b  # prev
            log_b = log_b - (A*b - B*b**2 - 1/self.c - C*log_b) / (A*b - 2*B*b**2 - C)
            b = np.exp(log_b)
            loop += 1
            if loop % 100 == 0:
                logging.warning('[%s] loop: %s, log_b: %s -> %s', self.name, loop, log_b1, log_b)
            #logging.info('updating b: %s', b)
        return b

    #cdef likelihood(self, x):
    #    return np.exp(- c/2 * (x - (a + b * Z_sum_mat)) ** 2) * theta_mat  # ~ z1, z2

    def _iter_z_probs(self):
        cdef DTYPE_t x1
        cdef int ploidy
        cdef np.ndarray likes
        #cdef np.ndarray[np.float_t, ndim=2] Z_sum_mat
        Z_sum_mat = self.Z_sum_mat

        for x1, ploidy in zip(self.x, self.ploidies):
            if ploidy == 0:
                yield np.array(0), (0,) * self.Z_dim
            elif ploidy == 1:
                likes = np.exp(- self.c/2 * (x1 - (self.a + self.b * self.Z)) ** 2) * self.theta # ~ z1, z2
                prob = likes / likes.sum()
                yield prob, tuple(prob)
            elif ploidy == 2:
                likes = np.exp(- self.c/2 * (x1 - (self.a + self.b * Z_sum_mat)) ** 2) * self.theta_mat  # ~ z1, z2
                prob = likes / likes.sum()
                yield prob, tuple(prob.sum(axis=1))

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

    def iteration(self):
        cdef int step = 1
        cdef DTYPE_t dummy_zero = 1e-120
        cdef list z_probs, z1_probs
        cdef DTYPE_t p
        cdef DTYPE_t ll_next, term1
        cdef DTYPE_t diff = 1
        cdef DTYPE_t ll = -float('inf')
        cdef np.ndarray z_prob
        cdef np.ndarray np_z1_probs
        cdef np.ndarray[long, ndim=2] Z_sum_mat, Z2_sum_mat
        Z_sum_mat = self.Z_sum_mat
        Z2_sum_mat = self.Z2_sum_mat

        while diff > self.max_diff:
            if self.step_limit and step > self.step_limit:
                raise AssertionError('The number of steps exceeds the limit')

            #z_likes = np.array([likelihood(x[n]) for n in enumerate(ploidies)])   # {x, z1, z2}
            #z_probs = np.array([like / like.sum() for like in z_likes])  # {x, z1, z2}
            #z1_probs = z_probs.sum(axis=2)   # {x, z1}
            #z1_exp = (Z[None, :] * z1_probs).sum(axis=1)    # {x}
            #d_mean = 2 * z1_exp  # {x}
            #d2_mean = ((Z[None, None, :] + Z[None, :, None]) ** 2 * z_probs).sum(axis=1).sum(axis=1)  # {x}
            z_probs = []
            z1_probs = []
            for zp, z1p in self._iter_z_probs():
                z_probs.append(zp)
                z1_probs.append(z1p)
            nd_z1_probs = np.array(z1_probs)  # {x, z1}
            self.d_mean = np.array([(Z_sum_mat * z_prob).sum() if p == 2 else (self.Z * z_prob).sum() for z_prob, p in zip(z_probs, self.ploidies)])
            self.d2_mean = np.array([(Z2_sum_mat * z_prob).sum() if p == 2 else (self.Z2 * z_prob).sum() for z_prob, p in zip(z_probs, self.ploidies)])

            self.theta = np.array([list(z_prob.sum(axis=0) + z_prob.sum(axis=1)) if p == 2 else list(z_prob)
                    for z_prob, p in zip(z_probs, self.ploidies) if p > 0]).sum(axis=0) / self.total_ploidy
            self.theta_mat = self.theta[:, None] * self.theta[None, :]
            #theta1 = z1_probs.sum(axis=0) / N
            #a = (x.sum() - b * d_mean.sum() - la0) / (N + la/c)
            #b = (((x - a) * d_mean).sum() + lb/c * b0) / (d2_mean.sum() + lb/c)
            self.a = self._get_next_a()
            x_a = self.x - self.a
            self.x_a_d_mean = x_a * self.d_mean
            self.b = self._get_next_b()

            term1 = (x_a**2 - 2 * self.b * self.x_a_d_mean + self.b**2 * self.d2_mean).sum()
            self.c = self.N / (term1 + 2 * self.lc)
            #theta = theta1

            ll_next = self.N/2. * np.log(self.c)
            #ll_next += ((- c/2. * (x[:, None, None] - (a + b * Z_sum_mat[None, :, :])) ** 2 + 2 * np.log(theta + dummy_zero)[None, :, None]) * z_probs).sum()
            #ll_next += - (z_probs * np.log(z_probs + dummy_zero)).sum()
            ll_next += - self.c/2. * term1
            ll_next += (self.ploidies[:, None] * np.log(self.theta + dummy_zero)[None, :] * nd_z1_probs).sum()
            ll_next += - sum((z_prob * np.log(z_prob + dummy_zero)).sum() for z_prob in z_probs)   # entropy
            #ll_next += - la/2. * a**2 - la0 * a
            ll_next += - self.la/2. * self.a**2 + np.log(norm.cdf(self.alpha * self.a))   # skew norm
            #ll_next += - lb/2. * (b - b0)**2
            ll_next += - np.log(self.b) - self.lb/2. * np.log(self.b)**2
            ll_next += - self.lc * self.c
            diff = ll_next - ll
            assert diff >= -1e-2, 'The lower bound of the log-likelihood should be increased (diff: {0})'.format(diff)
            ll = ll_next
            d = np.array(self.x - self.a) / self.b
            #logging.info('ll: %s (%s, %s, %s), %s', ll, diff, self.max_diff, diff>self.max_diff, Result(self.a, self.b, self.c, d, self.theta, step, ll))
            if step % 100 == 0:
                d = np.array(self.x - self.a) / self.b
                logging.info('[%s] ll: %s (%s), %s', self.name, ll, diff, Result(self.a, self.b, self.c, None, self.theta, step, ll))

            step += 1
            # x = a + b * d
            d = (self.x - self.a) / self.b

        return Result(self.a, self.b, self.c, d, self.theta, step, ll)


def normalize_em(*args, **kwds):
    normalizer = Normalizer(*args, **kwds)
    return normalizer.iteration()

