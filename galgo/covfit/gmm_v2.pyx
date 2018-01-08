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
        public int step
    cdef:
        #np.float_t[:] x, theta
        DTYPE_t a, b, c, la, lb, lc, max_diff, alpha, delta, omega
        int max_copy, step_limit, Z_dim, N, total_ploidy
        np.ndarray x, d_mean, d2_mean, x_a_d_mean, theta, ploidies, Z, Z2
        np.ndarray vars
        np.ndarray rc_mean, rcz_mean, rcz2_mean, x_a_rcz_mean
        np.ndarray theta_mat, Z_sum_mat, Z2_sum_mat, var_sum_mat   # 2 dim
        str name

    def __init__(self, x, ploidies, max_copy=3, a=0., b=1., c=10., la=0., lb=.3, lc=.3, delta=.3, theta=None, max_diff=1e-6, step_limit=0, name='default'):
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
        self.delta = delta
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

        self.vars = self.Z.copy().astype('float')
        self.vars[0] = self.delta
        self.var_sum_mat = self.Z_sum_mat.copy().astype('float')
        self.var_sum_mat[0, 0] = delta
        #logging.info(self.vars)
        #logging.info(self.var_sum_mat)

        #alpha = 20 # skew parameter of a
        #self.alpha = 1. # skew parameter of a
        #omega = 0.2   # learning rate of a
        #self.omega = 1.   # learning rate of a

    # OK
    cdef inline DTYPE_t _get_next_a(self):
        return ((self.rc_mean * self.x).sum() - self.b * self.rcz_mean.sum()) / (self.rc_mean.sum() + self.la / self.c)

    cdef DTYPE_t _get_next_b(self):
        cdef DTYPE_t A = self.x_a_rcz_mean.sum()
        cdef DTYPE_t B = self.rcz2_mean.sum()
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
        cdef np.ndarray likes, probs
        #cdef np.ndarray[np.float_t, ndim=2] Z_sum_mat
        Z_sum_mat = self.Z_sum_mat

        for x1, ploidy in zip(self.x, self.ploidies):
            if ploidy == 0:
                yield np.array(0), (0,) * self.Z_dim, (0,) * self.Z_dim
            elif ploidy == 1:
                likes = 1./np.sqrt(self.vars) * np.exp(- self.c/2./self.vars * (x1 - (self.a + self.b * self.Z)) ** 2) # ~ z1, z2
                probs = likes * self.theta
                probs /= probs.sum()
                yield probs, tuple(probs), likes
            elif ploidy == 2:
                likes = 1./np.sqrt(self.var_sum_mat) * np.exp(- self.c/2./self.var_sum_mat * (x1 - (self.a + self.b * Z_sum_mat)) ** 2) # ~ z1, z2
                probs = likes * self.theta_mat
                probs /= probs.sum()
                yield probs, tuple(probs.sum(axis=1)), likes

    def iter_cns(self):
        dip_copy = self.max_copy * 2

        for p, (zp, z1p, _) in zip(self.ploidies, self._iter_z_probs()):
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

    def iter_genotypes(self):
        dip_copy = self.max_copy * 2

        for p, (zp, z1p, likes) in zip(self.ploidies, self._iter_z_probs()):
            if p == 0:
                yield {'0': (1.0, 1.0)}
            elif p == 1:
                remain = dip_copy + 1 - len(zp)
                yield dict((str(n), (zp[n], likes[n])) for n in xrange(len(zp)))
            else:
                probs = {}
                for i in xrange(self.Z_dim):
                    probs['{0}/{1}'.format(i, i)] = (zp[i, i], likes[i, i])
                    for j in xrange(i):
                        probs['{0}/{1}'.format(i, j)] = (zp[i, j] + zp[j, i], likes[i, j])
                yield probs

    def iteration(self):
        cdef DTYPE_t dummy_zero = 1e-120
        cdef list z_probs, z1_probs
        cdef DTYPE_t p
        cdef DTYPE_t ll_next, term1
        cdef DTYPE_t diff = 1
        cdef DTYPE_t ll = -float('inf')
        cdef np.ndarray ln_rc_mean
        cdef np.ndarray z_prob
        cdef np.ndarray np_z1_probs
        cdef np.ndarray[long, ndim=2] Z_sum_mat, Z2_sum_mat,
        cdef np.ndarray[DTYPE_t, ndim=2] var_sum_mat
        Z_sum_mat = self.Z_sum_mat
        Z2_sum_mat = self.Z2_sum_mat
        var_sum_mat = self.var_sum_mat

        self.step = 1
        while diff > self.max_diff:
            if self.step_limit and self.step > self.step_limit:
                raise AssertionError('The number of steps exceeds the limit')

            #z_likes = np.array([likelihood(x[n]) for n in enumerate(ploidies)])   # {x, z1, z2}
            #z_probs = np.array([like / like.sum() for like in z_likes])  # {x, z1, z2}
            #z1_probs = z_probs.sum(axis=2)   
            #z1_exp = (Z[None, :] * z1_probs).sum(axis=1)    # {x}
            #d_mean = 2 * z1_exp  # {x}
            #d2_mean = ((Z[None, None, :] + Z[None, :, None]) ** 2 * z_probs).sum(axis=1).sum(axis=1)  # {x}
            z_probs = []   # {x, z1, z2}
            z1_probs = []  # {x, z1}
            for zp, z1p, _ in self._iter_z_probs():
                z_probs.append(zp)
                z1_probs.append(z1p)
            nd_z1_probs = np.array(z1_probs)  # {x, z1}
            self.rc_mean = np.array([(1./var_sum_mat * z_prob).sum() if p == 2 else (1./self.vars * z_prob).sum() for z_prob, p in zip(z_probs, self.ploidies)])
            self.rcz_mean = np.array([(1./var_sum_mat * Z_sum_mat * z_prob).sum() if p == 2 else (1./self.vars * self.Z * z_prob).sum() for z_prob, p in zip(z_probs, self.ploidies)])
            self.rcz2_mean = np.array([(1./var_sum_mat * Z2_sum_mat * z_prob).sum() if p == 2 else (1./self.vars * self.Z2 * z_prob).sum() for z_prob, p in zip(z_probs, self.ploidies)])
            ln_rc_mean = np.array([(- np.log(var_sum_mat) * z_prob).sum() if p == 2 else (- np.log(self.vars) * z_prob).sum() for z_prob, p in zip(z_probs, self.ploidies)])

            self.theta = np.array([list(z_prob.sum(axis=0) + z_prob.sum(axis=1)) if p == 2 else list(z_prob)
                    for z_prob, p in zip(z_probs, self.ploidies) if p > 0]).sum(axis=0) / self.total_ploidy
            self.theta_mat = self.theta[:, None] * self.theta[None, :]
            self.a = self._get_next_a()
            x_a = self.x - self.a
            self.x_a_rcz_mean = x_a * self.rcz_mean
            self.b = self._get_next_b()

            term1 = (x_a**2 * self.rc_mean - 2 * self.b * self.x_a_rcz_mean + self.b**2 * self.rcz2_mean).sum()
            self.c = self.N / (term1 + 2 * self.lc)

            ll_next = self.N/2. * np.log(self.c)
            ll_next += ln_rc_mean.sum() / 2.
            ll_next += - self.c/2. * term1
            ll_next += (self.ploidies[:, None] * np.log(self.theta + dummy_zero)[None, :] * nd_z1_probs).sum()
            ll_next += - sum((z_prob * np.log(z_prob + dummy_zero)).sum() for z_prob in z_probs)   # entropy
            ll_next += - self.la/2. * self.a**2 #+ np.log(norm.cdf(self.alpha * self.a))   # skew norm
            ll_next += - np.log(self.b) - self.lb/2. * np.log(self.b)**2
            ll_next += - self.lc * self.c
            diff = ll_next - ll
            assert diff >= -1e-2, 'The lower bound of the log-likelihood should be increased (diff: {0})'.format(diff)
            ll = ll_next
            d = np.array(self.x - self.a) / self.b
            #logging.info('ll: %s (%s, %s, %s), %s', ll, diff, self.max_diff, diff>self.max_diff, Result(self.a, self.b, self.c, d, self.theta, step, ll))
            if self.step % 100 == 0:
                d = np.array(self.x - self.a) / self.b
                logging.info('[%s] ll: %s (%s), %s', self.name, ll, diff, Result(self.a, self.b, self.c, None, self.theta, self.step, ll))

            self.step += 1
            # x = a + b * d
            d = (self.x - self.a) / self.b

        return Result(self.a, self.b, self.c, d, self.theta, self.step, ll)


def normalize_em(*args, **kwds):
    normalizer = Normalizer(*args, **kwds)
    return normalizer.iteration()
