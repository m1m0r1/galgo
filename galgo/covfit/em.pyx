#!/usr/bin/env python
from __future__ import print_function
import logging
import numpy as np
from collections import namedtuple
from counter import Counter
from scipy.stats import norm
#from depth_norm import em
cimport numpy as np
ctypedef np.float64_t DTYPE_t
#ctypedef float DTYPE_t

Result = namedtuple('Result', 'a b c d theta step ll')

#class Normalizer(object):
#    def __init__(self, *args, **kwds):
#        self._core = _Normalizer(*args, **kwds)
#
#    def iteration(self):
#        return self._core.iteration()
#
def normalize_em(x, ploidies, max_copy=3, a=0., b=1., c=10., b0=1., la=.3, lb=.3, lc=.3, max_diff=1e-6, step_limit=None):
    Result = namedtuple('Result', 'a b c d theta step ll')

    Z_dim = max_copy + 1
    Z = np.arange(Z_dim)
    N = len(x)
    theta = 1. * np.zeros(Z_dim)
    ploidies = np.array(ploidies)
    total_ploidy = ploidies.sum()

    # setting initial theta values
    tmp_hap_cns = [int(k) for k in np.round((x - a) / b / 2.)]
    counter = Counter(tmp_hap_cns)
    less = sum([counter[k] for k in counter.keys() if k < 0])
    more = sum([counter[k] for k in counter.keys() if k > max_copy])
    counter_trim = Counter(dict((k, counter[k]) for k in xrange(max_copy + 1)))
    counter_trim[0] += less
    counter_trim[max_copy] += more
    theta0 = 1. # smoothing
    theta = 1. * np.array([counter_trim[k] + theta0 for k in xrange(max_copy + 1)]) / (len(x) + theta0 * (max_copy + 1))  # initial value for theta
    logging.info('Initial values: a=%s, b=%s, c=%s, theta=%s', a, b, c, theta)

    #theta[1] = 1.  # default is 1-copy for haploid
    #z_sum_mat = np.array([[z1 + z2 for z2 in Z] for z1 in Z]) # [(z1 + z2)]
    Z_sum_mat = Z[:, None] + Z[None, :]
    Z2 = Z ** 2
    Z2_sum_mat = Z_sum_mat ** 2

    # testing
    la0 = 0
    #la = 0
    alpha = 20 # skew parameter of a
    alpha = 1. # skew parameter of a
    #omega = 0.2   # learning rate of a
    omega = 1.   # learning rate of a

    step = 1
    diff = 1
    ll = -float('inf')
    def get_next_a(a):
        a1 = float('inf')
        while abs(a1 - a) > 1e-3:
            rate = norm.pdf(alpha * a) / norm.cdf(alpha * a)
            a1 = a  # prev
            a = a * (1. - omega) + omega * (c * (x.sum() - b * d_mean.sum()) + alpha * rate) / (c * N + la)
            #logging.info('updating a: %s, rate: %s', a, rate)
        return a

    def get_next_b(b):
        A = x_a_d_mean.sum()
        B = d2_mean.sum()
        C = lb/c
        #b1 = float('inf')
        #while abs(b1 - b) > 1e-3:
        #    b1 = b  # prev
        #    b = b - b * (A*b - B*b**2 - 1/c - C*np.log(b)) / (A*b - 2*B*b**2 - C)
        #    #logging.info('updating b: %s', b)

        log_b1 = float('inf')
        log_b = np.log(b)
        while abs(log_b1 - log_b) > 1e-3:
            log_b1 = log_b  # prev
            log_b = log_b - (A*b - B*b**2 - 1/c - C*log_b) / (A*b - 2*B*b**2 - C)
            b = np.exp(log_b)
            #logging.info('updating b: %s', b)
        return b

    def likelihood(x):
        return np.exp(- c/2 * (x - (a + b * Z_sum_mat)) ** 2) * theta_mat  # ~ z1, z2

    def iter_z_probs(x, ploidies):
        for x1, ploidy in zip(x, ploidies):
            if ploidy == 0:
                yield np.array(0), (0,) * Z_dim
            elif ploidy == 1:
                likes = np.exp(- c/2 * (x1 - (a + b * Z)) ** 2) * theta # ~ z1, z2
                prob = likes / likes.sum()
                yield prob, tuple(prob)
            elif ploidy == 2:
                likes = np.exp(- c/2 * (x1 - (a + b * Z_sum_mat)) ** 2) * theta_mat  # ~ z1, z2
                prob = likes / likes.sum()
                yield prob, tuple(prob.sum(axis=1))

    while diff > max_diff:
        if step_limit and step > step_limit:
            raise AssertionError('The number of step should not exceed limit')

        theta_mat = theta[:, None] * theta[None, :]
        #z_likes = np.array([likelihood(x[n]) for n in enumerate(ploidies)])   # {x, z1, z2}
        #z_probs = np.array([like / like.sum() for like in z_likes])  # {x, z1, z2}
        #z1_probs = z_probs.sum(axis=2)   # {x, z1}
        #z1_exp = (Z[None, :] * z1_probs).sum(axis=1)    # {x}
        #d_mean = 2 * z1_exp  # {x}
        #d2_mean = ((Z[None, None, :] + Z[None, :, None]) ** 2 * z_probs).sum(axis=1).sum(axis=1)  # {x}
        z_probs = []
        z1_probs = []
        for zp, z1p in iter_z_probs(x, ploidies):
            z_probs.append(zp)
            z1_probs.append(z1p)
        z1_probs = np.array(z1_probs)  # {x, z1}
        d_mean = np.array([(Z_sum_mat * z_prob).sum() if p == 2 else (Z * z_prob).sum() for z_prob, p in zip(z_probs, ploidies)])
        d2_mean = np.array([(Z2_sum_mat * z_prob).sum() if p == 2 else (Z2 * z_prob).sum() for z_prob, p in zip(z_probs, ploidies)])
        theta = np.array([list(z_prob.sum(axis=0) + z_prob.sum(axis=1)) if p == 2 else list(z_prob)
                for z_prob, p in zip(z_probs, ploidies) if p > 0]).sum(axis=0) / total_ploidy
        #theta1 = z1_probs.sum(axis=0) / N
        #a = (x.sum() - b * d_mean.sum() - la0) / (N + la/c)
        #b = (((x - a) * d_mean).sum() + lb/c * b0) / (d2_mean.sum() + lb/c)
        a = get_next_a(a)
        x_a = x - a
        x_a_d_mean = x_a * d_mean
        b = get_next_b(b)

        term1 = (x_a**2 - 2 * b * x_a_d_mean + b**2 * d2_mean).sum()
        c = N / (term1 + 2 * lc)
        #theta = theta1

        dummy_zero = 1e-120
        ll_next = N/2. * np.log(c)
        #ll_next += ((- c/2. * (x[:, None, None] - (a + b * Z_sum_mat[None, :, :])) ** 2 + 2 * np.log(theta + dummy_zero)[None, :, None]) * z_probs).sum()
        #ll_next += - (z_probs * np.log(z_probs + dummy_zero)).sum()
        ll_next += - c/2. * term1
        ll_next += (ploidies[:, None] * np.log(theta + dummy_zero)[None, :] * z1_probs).sum()
        ll_next += - sum((z_prob * np.log(z_prob + dummy_zero)).sum() for z_prob in z_probs)   # entropy
        #ll_next += - la/2. * a**2 - la0 * a
        ll_next += - la/2. * a**2 + np.log(norm.cdf(alpha * a))   # skew norm
        #ll_next += - lb/2. * (b - b0)**2
        ll_next += - np.log(b) - lb/2. * np.log(b)**2
        ll_next += - lc * c
        diff = ll_next - ll
        assert diff >= -1e-2, 'the lower bound of the log-likelihood should be increased (diff: {0})'.format(diff)
        ll = ll_next
        step += 1
        if step % 100 == 0:
            d = (x - a) / b
            logging.info('ll: %s, %s', ll, Result(a, b, c, d, theta, step, ll))

        #print (theta)
        #print (a)
        #print (b)
        #print (c)
        #print (ll)
        #print (diff)
        #if step > 1000:
        #    break

    # x = a + b * d
    d = (x - a) / b

    return Result(a, b, c, d, theta, step, ll)
