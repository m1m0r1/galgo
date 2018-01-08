from builtins import zip
import numpy as np


def get_banded_matrix(band1, band2, m, n, dtype=float):
    """
    >>> b = get_banded_matrix(1, 1, 2, 3, dtype=int)
    >>> b[1, 2] = 1
    >>> b[1, 2]
    1
    >>> b = get_banded_matrix(1, 1, 3, 2, dtype=int)
    >>> b[2, 1] = 1
    >>> b[2, 1]
    1
    >>> b[1, 2]
    0
    """
    assert band1 > 0 and band2 > 0
    data = np.zeros((band1 + band2 + 1, max(m, n)), dtype=dtype)
    return sparse_banded_matrix(data, -band1, m, n)


# similar interface to scipy.sparse.spdiags
class sparse_diag_matrix:
    """
    >>> data = np.array([[1, 2, 3, 4], [1, 2, 3, 4], [1, 2, 3, 4]])
    >>> diags = np.array([0, -1, 2])
    >>> sp = sparse_diag_matrix(data, diags, 4, 4)
    >>> sp.toarray()
    array([[1, 0, 3, 0],
           [1, 2, 0, 4],
           [0, 2, 3, 0],
           [0, 0, 3, 4]])
    >>> sp[0, 0]
    1
    >>> sp[2, 1]
    2
    >>> sp[3, 1]
    0
    """
    def __init__(self, data, diags, m, n):
        self._diags = dict((j1, vals) for j1, vals in zip(diags, data))
        self.shape = (m, n)

    def __getitem__(self, key):
        if isinstance(key, tuple):
            i1, i2 = key
            if i1 < 0 or i2 < 0:
                raise NotImplementedError('Cannot handle negative index')
            j1 = i2 - i1   # diag index
            j2 = i2
            d = self._diags.get(j1)
            return d[j2] if d is not None else 0
        raise NotImplementedError

    def __setitem__(self, key, value):
        if isinstance(key, tuple):
            i1, i2 = key
            if i1 < 0 or i2 < 0:
                raise NotImplementedError('Cannot handle negative index')
            j1 = i2 - i1   # diag index
            j2 = i2
            self._diags[j1][j2] = value   # cannot set values to non-initialized positions
            return
        raise NotImplementedError

    def toarray(self):  # adhoc
        m, n = self.shape
        from scipy.sparse import spdiags
        data = []
        diags = []
        for j1, dat in self._diags.items():
            data.append(dat)
            diags.append(j1)
        return spdiags(data, diags, m, n).toarray()


class sparse_banded_matrix:
    """
    >>> data = np.array([[1, 2, 3, 4], [1, 2, 3, 4], [1, 2, 3, 4]])
    >>> offset = -1
    >>> sp = sparse_banded_matrix(data, offset, 4, 4)
    >>> sp.toarray()
    array([[1, 2, 0, 0],
           [1, 2, 3, 0],
           [0, 2, 3, 4],
           [0, 0, 3, 4]])
    >>> sp[0, 0]
    1
    >>> sp[2, 1]
    2
    >>> sp[3, 1]
    0
    """
    def __init__(self, data, offset, m, n):
        self.shape = (m, n)
        self._offset = offset
        self._data = data
        assert data.shape[1] >= max(m, n), 'data shape should be longer than matrix length'
        self._width = data.shape[0]

    def __getitem__(self, key):   # FIXME not so fast
        if isinstance(key, tuple):
            i1, i2 = key
            if i1 < 0 or i2 < 0:
                raise NotImplementedError('Cannot handle negative index')
            j1 = i2 - i1 - self._offset   # diag index
            j2 = i2
            return self._data[j1, j2] if 0 <= j1 < self._width else 0
        raise NotImplementedError

    def __setitem__(self, key, value):
        if isinstance(key, tuple):
            i1, i2 = key
            if i1 < 0 or i2 < 0:
                raise NotImplementedError('Cannot handle negative index')
            j1 = i2 - i1 - self._offset  # diag index
            j2 = i2
            self._data[j1, j2] = value   # cannot set values to non-initialized positions
            return
        raise NotImplementedError

    def toarray(self):  # adhoc
        m, n = self.shape
        from scipy.sparse import spdiags
        data = self._data
        diags = list(range(self._offset, self._offset + self._width))
        return spdiags(data, diags, m, n).toarray()
