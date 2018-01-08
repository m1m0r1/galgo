from __future__ import absolute_import
import io

def open(file, *args, **kwds):
    if file == '-':
        file = '/dev/stdin'
    return io.open(file, *args, **kwds)
