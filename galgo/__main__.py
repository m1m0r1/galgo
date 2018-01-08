#!/usr/bin/env python
from __future__ import absolute_import
import os
import sys
if __package__ == '':   # hack for relative import
     (search_path, __package__) = os.path.split(os.path.split(os.path.abspath(__file__))[0])
     sys.path.insert(0, search_path)
     __import__(__package__)

from . import __VERSION__
import logging
from argtools import command
from .tools import fasta
from .tools import seq_sim
from .tools import msf
from .tools import msa
from .tools import bam
from .tools import bam_profile
from .tools import bam_summary
from .tools import bam_summary2
from .tools import blast
from .tools import vcf
from .tools import assoc
from .tools import feature
from .tools import imgt
from .tools import segdup
#from .tools import bam_aln
from .tools import gc_profile
#from .tools import covfit_gmm
from .tools import covfit_gmm_v2
#from .tools import covfit_pmm
from .tools import igv
from .tools import misc

command.before_run(lambda self: logging.info('galgo version: %s', __VERSION__))
def main():
    command.run()

if __name__ == '__main__':
    main()
