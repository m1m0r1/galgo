from __future__ import print_function
from argtools import command, argument
import logging
import pysam


def iter_tabs(fp):
    for line in fp:
        if line.startswith('#'):
            continue
        yield line.rstrip('\r\n').split('\t')


@command.add_sub
@argument('region_bed')
@argument('bam')
@argument('--use-1start', action='store_true', default=False, help='regard region bed as 1 start and 1 end')
@argument('-m', '--mass', choices=['middle', 'start', 'end'], default='middle', help='mass point of read to determine the read is inside of the region')
def bam_coverage(args):
    """
    The following reads will be omitted
    - secondary alignment
    - supplementary alignment
    - unmapped

    Append following columns to input bed

    coverage of
    - all reads
    - pair-mapped
    - MQ0
    - MQ0 pair-mapped
    """
    offset = 1 if args.use_1start else 0

    get_mass = {
            'start': lambda rec: rec.pos,
            'end':   lambda rec: rec.aend - 1,   # end position of the alignment (note that aend points one past the last aligned residue)
            'middle': lambda rec: (rec.pos + rec.aend - 1) / 2.,   # middle position of the alignment
    }[args.mass]

    with pysam.Samfile(args.bam, 'rb') as samfile:
        for bed in iter_tabs(open(args.region_bed)):
            ref = bed[0]
            start = int(bed[1]) - offset
            end = int(bed[2])
            logging.info('fetching %s:%s-%s', ref, start, end)

            bam_recs = samfile.fetch(reference=ref, start=start, end=end)
            total = 0
            paired = 0
            total0 = 0
            paired0 = 0
            for rec in bam_recs:
                if rec.flag & 0x904:  # skip unmapped, secondary and supplementary alignment
                    continue
                mass = get_mass(rec)
                if not (start <= mass < end):  # skip if mass point is not included
                    continue

                is_paired = \
                    not rec.is_unmapped \
                    and not rec.mate_is_unmapped \
                    and rec.is_proper_pair

                total += 1
                if is_paired:
                    paired += 1

                if rec.mapq == 0:
                    total0 += 1
                    if is_paired:
                        paired0 += 1

            covs = (total, paired, total0, paired0)
            print (*(tuple(bed) + covs), sep='\t')
