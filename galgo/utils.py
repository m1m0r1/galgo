from __future__ import print_function
try:
    from collections import OrderedDict
except ImportError:
    from ordereddict import OrderedDict
try:
    from collections import Counter
except ImportError:
    from counter import Counter
try:
    from UserDict import DictMixin
    class _attrdict(object, DictMixin):
        pass
except ImportError:
    from collections import MutableMapping
    class _attrdict(MutableMapping):
        pass
from builtins import zip, range
from six.moves import zip_longest
from itertools import chain
import os
import re
import tempfile
import logging
import warnings
import numpy as np
import functools
from .sselogsumexp import logsumexp

_term_colors = {
    'clear': '\033[0m',
    'black': '\033[30m',
    'red': '\033[31m',
    'green': '\033[32m',
    'yellow': '\033[33m',
    'blue': '\033[34m',
    'purple': '\033[35m',
    'cyan': '\033[36m',
    'white': '\033[37m'
}


class attrdict(_attrdict):
    """ Dictionary class whose properties are accessible as attribute

    #TODO implement merge method if needed

    >>> d0 = attrdict()  # empty
    >>> d0
    {}

    >>> d1 = attrdict([('a', 1), ('b', {'key': 'value'}), ('c', 3)], __recursive__=True)  # init with (key, value) pairs
    >>> d1.a, d1.b, d1.b.key, d1.c                                    # access with attribute
    (1, {'key': 'value'}, 'value', 3)
    >>> d1['a']
    1

    >>> d2 = attrdict({'a': 1, 'b': 2, 'd': ['e', {'f': 2}]}, __recursive__=True)         # init with dict
    >>> d3 = attrdict(a=1, b=2, d=['e', {'f': 2}], __recursive__=True)                    # init with kwds style
    >>> d2 == d3
    True
    >>> d3.d[1]
    {'f': 2}

    >>> d1.update(d2)   # override dict
    >>> d1.a, d1.b, d1.c, d1.d
    (1, 2, 3, ['e', {'f': 2}])
    """
    def __init__(self, *args, **kwds):
        if args:
            self.__dict__ = dict(args[0])
        if kwds:
            self.__dict__.update(**kwds)

        recursive = kwds.pop('__recursive__', False)
        if recursive:
            for (k, v) in self.items():
                if isinstance(v, dict):
                    self[k] = attrdict(v)
                elif isinstance(v, (list, tuple, set, frozenset)):
                    self[k] = v.__class__((attrdict(el) if isinstance(el, dict) else el) for el in v)

    def __iter__(self):
        return iter(self.__dict__)

    def __getitem__(self, key):
        return getattr(self, key)

    def __setitem__(self, key, value):
        return setattr(self, key, value)

    def __delitem__(self, key):
        return delattr(self, key)

    def __len__(self):
        return len(self.__dict__)

    def __repr__(self):
        return repr(dict(self.items()))


class defaultattrdict(attrdict):
    """ Dictionary class whose properties are accessible as attribute

    #TODO implement merge method if needed

    >>> d0 = defaultattrdict(dict)  # empty
    >>> d0
    {}

    >>> d1 = defaultattrdict(list)
    >>> d1.a.append(1)
    >>> d1.a.append(2)
    >>> d1['a']
    [1, 2]
    >>> d1.a
    [1, 2]
    >>> d1.b
    []
    """
    def __init__(self, default):
        self._default = default
        super(defaultattrdict, self).__init__()

    def __getitem__(self, key):
        return getattr(self, key)

    def __getattr__(self, key):
        if key.startswith('_'):
            return self.__dict__[key]
        try:
            return super(defaultattrdict, self).__getattr__(key)
        except AttributeError:
            setattr(self, key, self._default())
            return self[key]

    def __iter__(self):
        return (k for k in iter(self.__dict__) if not k.startswith('_'))


def color_term(term, color):
    return ''.join((_term_colors[color], term, _term_colors['clear']))


def fill_text(text, length, char=' '):
    """
    >>> fill_text('abc', 5, char='+')
    'abc++'
    >>> fill_text('abc', 2, char='+')
    'abc'
    """
    return (text + char * length)[:max(len(text), length)]

rfill_text = fill_text

def lfill_text(text, length, char=' '):
    """
    >>> lfill_text('abc', 5, char='+')
    '++abc'
    >>> lfill_text('abc', 2, char='+')
    'abc'
    """
    return (char * length + text)[-max(len(text), length):]


def remove_suffix(text, suffix):
    """
    >>> remove_suffix('NA12878.sorted.bam', '.bam')
    'NA12878.sorted'
    >>> remove_suffix('NA12878.sorted.bam', '.sam')
    'NA12878.sorted.bam'
    >>> remove_suffix('NA12878.sorted.bam', '')
    'NA12878.sorted.bam'
    """
    if suffix and text.endswith(suffix):
        return text[:-len(suffix)]
    else:
        return text


def iter_tabs(it, sep='\t', skip_if=None):
    if skip_if is None:
        skip_if = lambda line : False

    for line in it:
        if skip_if(line):
            continue
        yield line.rstrip('\r\n').split(sep)

def with_header(rows, container=OrderedDict, use_header=False):
    """
    >>> tab = iter_tabs(['a,b,c', '1,2,3', '4,5,6,7'], sep=',')
    >>> recs = list(with_header(tab))
    >>> recs[0]['a']
    '1'
    >>> list(recs[1].values())
    ['4', '5', '6']
    """
    rows = iter(rows)
    keys = next(rows)
    if use_header:
        yield tuple(keys)
    for vals in rows:
        yield container(zip(keys, vals))


def _safetype(type):
    def f(v, default=None):
        try:
            return type(v)
        except Exception as e:
            return default
    return f

safeint = _safetype(int)
safefloat = _safetype(float)

def safediv(a, b, error=float('nan')):
    try:
        return a / b
    except ZeroDivisionError:
        return error


def weighted_choice(weights, size=1):
    """ weights is an instance of numpy.array
    >>> np.random.seed(0)
    >>> dict(Counter(weighted_choice(np.array([1, 0, 0]), 10000))) == {0: 10000}
    True
    >>> np.random.seed(0)
    >>> dict(Counter(weighted_choice(np.array([1.5, 0, 3.0, 0.5]), 10000))) == {0: 3055, 2: 5943, 3: 1002}
    True
    """
    t = weights.cumsum()
    s = weights.sum()
    return np.searchsorted(t, np.random.rand(size) * s)


def getmax(iterable, key=lambda x: x):
    """
    >>> getmax([0, 1, 2])
    2
    >>> getmax([0, 1, 2], key=lambda x: -x)
    0
    """
    el_scores = [(el, key(el)) for el in iterable]
    el_score = max(el_scores, key=lambda x: x[1])
    return el_score[0]


_missing = object()
def chunked(seq, size, fillvalue=_missing):
    """
    >>> [tuple(lis) for lis in chunked([1, 2, 3, 4, 5, 6, 7], 3)]
    [(1, 2, 3), (4, 5, 6)]
    >>> [tuple(lis) for lis in chunked([1, 2, 3, 4, 5, 6, 7], 3, fillvalue=None)]
    [(1, 2, 3), (4, 5, 6), (7, None, None)]
    """
    if fillvalue is _missing:
        return zip(*([iter(seq)] * size))
    else:
        return zip_longest(fillvalue=fillvalue, *([iter(seq)] * size))


def bucket(it, size):
    """
    >>> [tuple(lis) for lis in bucket([1, 2, 3, 4, 5, 6, 7], 3)]
    [(1, 2, 3), (4, 5, 6), (7,)]
    """
    assert size > 0
    i = 0
    stock = []
    for el in it:
        stock.append(el)
        i += 1
        if i == size:
            yield stock
            stock = []
            i = 0
    if stock:
        yield stock


def unzip(it, size):
    """
    >>> a, b = unzip(zip((1, 2, 3, 4, 5), ('a', 'b', 'c')), 2)
    >>> a
    [1, 2, 3]
    >>> b
    ['a', 'b', 'c']
    """
    els = tuple([] for _ in range(size))
    for el in it:
        for i in range(size):
            els[i].append(el[i])
    return els


def iter_pairs(seq, last=None):
    """
    >>> list(iter_pairs((1, 2, 3, 4, 5)))
    [(1, 2), (2, 3), (3, 4), (4, 5), (5, None)]
    >>> list(iter_pairs((1,)))
    [(1, None)]
    """
    it = iter(seq)
    prev = next(it)
    for el in it:
        yield (prev, el)
        prev = el
    yield (prev, last)


def iter_counts(seq, same=lambda a, b: a==b):
    """
    >>> list(iter_counts(('a', 'a', 'b', 'A', 'a', None, 'c')))
    [('a', 2), ('b', 1), ('A', 1), ('a', 1), (None, 1), ('c', 1)]
    >>> list(iter_counts(('a', 'a', 'b', 'A', 'a', 'c'), same=lambda a, b: a.upper()==b.upper()))
    [('a', 2), ('b', 1), ('A', 2), ('c', 1)]
    """
    it = iter(seq)
    try:
        prev = next(it)
    except StopIteration:
        return

    count = 1
    for el in it:
        if same(prev, el):
            count += 1
        else:
            yield (prev, count)
            prev = el
            count = 1
    yield (prev, count)


def iter_uniques(iterable):
    """ Creates an ordered set from a list of tuples or other hashable items

    >>> ''.join(iter_uniques('abracadabra'))
    'abrcd'
    >>> ''.join(iter_uniques('simsalabim'))
    'simalb'
    """
    seen = set()

    for item in iterable:
        #Save unique items in input order
        if item not in seen:
            seen.add(item)
            yield item


def collect_while(pred, it):
    """
    >>> cols, it = collect_while(lambda x: x=='A', 'AAABAAAA')
    >>> ''.join(cols)
    'AAA'
    >>> next(it)
    'B'
    >>> cols, it = collect_while(lambda x: x=='A', it)
    >>> ''.join(cols)
    'AAAA'
    >>> next(it, 'end')
    'end'
    """
    it = iter(it)
    lis = []
    try:
        while 1:
            el = next(it)
            if not pred(el):
                break
            lis.append(el)
        return lis, chain((el,), it)
    except StopIteration as e:
        if lis:
            return lis, it
        else:
            raise e

def skip_until(pred, it):
    """
    >>> el, it = skip_until(lambda x: x=='B', 'AAABAAAA')
    >>> el
    'B'
    >>> el, it = skip_until(lambda x: x=='A', it)
    >>> el
    'A'
    >>> ''.join(it)
    'AAA'
    """
    cols, it = collect_while(make_not(pred), it)
    el = next(it)
    return el, it

_blank_pat1 = re.compile('^\s*$')
_blank_pat2 = re.compile('\s+')

def isblank(line):
    '''
    >>> all(map(isblank, [' ', '', '  ', '\\n']))
    True
    >>> any(map(isblank, [' A ', 'A', ' ABC', 'ABCD \\n']))
    False
    '''
    return bool(_blank_pat1.match(line))

def make_not(fn):
    """ Invert return value of input function
    >>> notblank = make_not(isblank)
    >>> any(map(notblank, [' ', '', '  ', '\\n']))
    False
    >>> all(map(notblank, [' A ', 'A', ' ABC', 'ABCD \\n']))
    True
    """
    return lambda *args, **kwds: not fn(*args, **kwds)

def blank_split(string, maxsplit=0):
    """ 
    >>> blank_split('A   B  C')
    ['A', 'B', 'C']
    >>> blank_split('   A  B C   ')
    ['A', 'B', 'C']
    """
    return _blank_pat2.split(string.strip(), maxsplit=maxsplit)


# decorator
def deprecated(fn):
    def g(*args, **kwds):
        warnings.warn("deprecated", DeprecationWarning)
        return fn(*args, **kwds)
    return g


def memoize(fn, *args):
    """ A simple memoize function
    """
    _d = {}
    @functools.wraps(fn)
    def f(*args):
        if args not in _d:
            _d[args] = fn(*args)
        return _d[args]
    return f


# taken from https://github.com/mitsuhiko/werkzeug
class cached_property(object):
    """A decorator that converts a function into a lazy property.  The
    function wrapped is called the first time to retrieve the result
    and then that calculated result is used the next time you access
    the value::

        class Foo(object):

            @cached_property
            def foo(self):
                # calculate something important here
                return 42

    The class has to have a `__dict__` in order for this property to
    work.
    """ 

    def __init__(self, func, name=None, doc=None):
        self.__name__ = name or func.__name__
        self.__module__ = func.__module__
        self.__doc__ = doc or func.__doc__
        self.func = func

    def __get__(self, obj, type=None):
        if obj is None:
            return self
        value = obj.__dict__.get(self.__name__, _missing)
        if value is _missing:
            value = self.func(obj)
            obj.__dict__[self.__name__] = value
        return value


def get_empty_port(default=None, port_range=None):
    """ get port
    """
    import psutil
    conns = psutil.net_connections()
    used_ports = set([conn.laddr[1] for conn in conns])
    if default and default not in used_ports:
        return default
    if port_range is None:
        port_range = list(range(49152, 65535+1))  # ephemeral ports
    for port in port_range:
        if port not in used_ports:
            return port


def file_split(input, nlines=1000, prefix='', suffix='', dir=None, header=None):
    """
    >>> data = range(50)
    >>> d = tempfile.mkdtemp()
    >>> results = list(file_split(data, nlines=25, prefix='abc.', suffix='.num', dir=d))
    >>> [os.path.basename(f) for f in results]
    ['abc.1.num', 'abc.2.num']
    >>> open(results[0]).read() == ''.join(map(str, range(25)))
    True
    >>> open(results[1]).read() == ''.join(map(str, range(25, 50)))
    True
    >>> os.remove(results[0])
    >>> os.remove(results[1])
    >>> os.removedirs(d)
    """
    assert nlines > 0
    if dir is None:
        dir = tempfile.mkdtemp()
    logging.info('Write split files in {0}'.format(dir))

    def get_next(n):
        name = os.path.join(dir, '{0}{1}{2}'.format(prefix, n, suffix))
        writer = open(name, 'w+')
        if header is not None:
            print (header, file=writer)
        return name, writer

    name = writer = None
    try:
        n = 1
        l = 0
        name, writer = get_next(n)
        for line in input:
            if l == nlines:
                writer.close()
                yield name
                n += 1
                name, writer = get_next(n)
                l = 0
            print (line, end='', file=writer)
            l += 1
        writer.close()
        writer = None
        yield name
    finally:
        if writer is not None:
            writer.close()


def format(st, *args, **kwds):
    """ Format with string interpolation (currently, only local variables is usable)

    >>> a = 1
    >>> b = 2
    >>> format('{a} + {a} = {b}')
    '1 + 1 = 2'
    >>> format('{a} + {a} != {b}', b=3)
    '1 + 1 != 3'
    >>> format('{a} + {0} = {b}', 2, b=3)
    '1 + 2 = 3'

    # TODO impossible?
    # >>> def nested1():
    # ...     a = 2
    # ...     def nested2():
    # ...         return format('{a} = {a} = {b}')
    # ...     return nested2
    # >>> nested1()()
    # '2 = 2 = 2'
    """
    import inspect
    from itertools import islice
    frame = inspect.currentframe()
    try:
        d = {}
        #d.update(frame.f_back.f_globals)
        d.update(frame.f_back.f_locals)
        d.update(kwds)
        return st.format(*args, **d)
    finally:
        del frame


def get_process_pool(njobs, init=None):
    from multiprocess import Pool
    import signal
    def init_pool():
        logging.info('Init pool')
        signal.signal(signal.SIGINT, signal.SIG_IGN)
    if init is None:
        init = init_pool
    pool = Pool(njobs, init_pool)
    return pool
