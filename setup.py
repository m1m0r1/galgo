import os
import sys
import glob
from setuptools import setup, find_packages
from setuptools.extension import Extension
#from distutils.core import setup

REQUIRES = ['argtools>=0.1.5', 'numpy', 'matplotlib', 'pandas>=0.19', 'sklearn', 'statsmodels', 'multiprocess', 'joblib', 'future', 'xlsxwriter', 'tqdm']
REQUIRES.extend(['cython>=0.22',  'pysam>=0.11', 'networkx'])
if sys.version_info >= (2, 7):
    REQUIRES.append('seaborn')

try:
    from collections import OrderedDdict
except ImportError:
    REQUIRES.append('ordereddict')
try:
    from collections import Counter
except ImportError:
    REQUIRES.append('counter')

README = os.path.join(os.path.dirname(__file__), 'README.md')
if os.path.exists(README):
    long_description = open(README).read() + '\n\n'
else:
    long_description = 'A collection of algorithms and tools for genome analysis'
from galgo import __VERSION__

import numpy as np
if 'build_ext' in sys.argv:
    from Cython.Distutils import build_ext
    from Cython.Compiler import Options
    from Cython.Build import cythonize
    Options.fast_fail = True
    exts = (
            cythonize('galgo/covfit/*.pyx')
         + cythonize('galgo/alignment/*.pyx')
         + cythonize('galgo/csamutil.pyx')
         + [Extension('galgo.sselogsumexp', sources=['galgo/sselogsumexp.pyx', 'src/logsumexp.c'], depends=glob.glob('src/*.h'), extra_compile_args=['-msse2'])]
    )
    cmd_class = {'build_ext': build_ext}
else:
    exts = [
            Extension('galgo.covfit.pmm', sources=['galgo/covfit/pmm.c']),
            Extension('galgo.covfit.gmm_v2', sources=['galgo/covfit/gmm_v2.c']),
            Extension('galgo.alignment.basics', sources=['galgo/alignment/basics.c']),
            Extension('galgo.csamutil', sources=['galgo/csamutil.c']),
            Extension('galgo.sselogsumexp', sources=['galgo/sselogsumexp.c', 'src/logsumexp.c'], depends=glob.glob('src/*.h'), extra_compile_args=['-msse2'])
            ]
    cmd_class = {}

exclude = ['galgo.__main__']
setup(name='galgo',
      version=__VERSION__,
      install_requires=REQUIRES,
      description='A collection of algorithms and tools for genome analysis',
      long_description=long_description,
      classifiers=['Programming Language :: Python',
                   'Topic :: Scientific/Engineering :: Bio-Informatics'],
      keywords='genome bio-informatics fasta',
      author='Takahiro Mimori',
      author_email='takahiro.mimori@gmail.com',
      #py_modules=['galgo'],
      package_data = {'': ['README.md']},
      #package_dir = {'galgo': 'galgo'},
      packages=find_packages(exclude=exclude),
      #scripts = ['bin/galgo'],
      scripts = ['scripts/bed2chain.py'],
      entry_points={'console_scripts': ['galgo = galgo.__main__:main']},
      url='https://github.com/m1m0r1/galgo',
      cmdclass = cmd_class,
      include_dirs = [np.get_include(), 'src'],
      ext_modules = exts,
)
