import numpy as np
cimport numpy as np

cdef extern from "logsumexp.h":
    float flogsumexp(const float* buf, int N) nogil
    double dlogsumexp(const double* buf, int N) nogil   # redundant?


cdef logsumexpf(np.ndarray[dtype=np.float32_t] a):
    """
    Compute the log of the sum of exponentials of input elements.
    Parameters
    ----------
    a : np.ndarray
        Input data. Must be contiguous.
    """
    if not (a.flags['C_CONTIGUOUS'] or a.flags['F_CONTIGUOUS']):
        raise TypeError('a must be contiguous')

    return flogsumexp(&a[0], a.size)


cdef logsumexpd(np.ndarray[dtype=np.float64_t] a):
    """
    Compute the log of the sum of exponentials of input elements.
    Parameters
    ----------
    a : np.ndarray
        Input data. Must be contiguous.
    """
    if not (a.flags['C_CONTIGUOUS'] or a.flags['F_CONTIGUOUS']):
        raise TypeError('a must be contiguous')

    return dlogsumexp(&a[0], a.size)


def logsumexp(np.ndarray a):
    if a.dtype == np.float32:
        return logsumexpf(a)
    if a.dtype == np.float64:
        return logsumexpd(a)
    raise Exception('Unknown')
