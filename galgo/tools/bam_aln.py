from __future__ import print_function
from argtools import command, argument
from pysam import Samfile


@command.add_sub
@argument('vcf')
@argument('bam')
@argument('--fasta')
def bam_variant_aln(args):
    samfile = Samfile(args.bam)
    for rec in args.vcf:
        fp = open(vcf)
        reader = pyvcf.Reader(fp)
        self.positions = []

    for rec in samfile.fetch(vcf):
        samfile.getrname(rec.tid)
        rec
