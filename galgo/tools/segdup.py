from __future__ import print_function
import os
import glob
from argtools import command, argument
from .. import sh
from ..utils import iter_tabs, file_split
from ..blast import BlastTabFile
import logging
import numpy as np
from ..job_queue import UGEQueue


def random_rgb():
    return (np.random.random_integers(256),
            np.random.random_integers(256),
            np.random.random_integers(256))


#@command.add_sub
@argument('segdup_tab')
def segdup2trackbed(args):
    """
    * set unique name to pairs (use uid)
    * randomize color by name
    """
    with open(args.segdup_tab) as fp:
        rows = iter_tabs(fp)
        header = with_header(rows)
        for row in rows:
            chrom = row['chrom']
            start = row['chromStart']
            end = row['chromEnd']
            name = row['name']
            score = row['score']
            strand = row['strand']
            other_chrom = row['otherChrom']
            other_start = row['otherStart']
            other_end = row['otherEnd']


def require(file):
    if not os.path.exists(file):
        raise Exception('File does not exist: {0}'.format(file))

def is_valid_file(file, min_size=1):
    if not os.path.exists(file):
        return False
    if os.path.getsize(file) < min_size:
        return False
    return True

def remove_all(*targets):
    for target in targets:
        if os.path.exists(target):
            sh.call(['rm', '-r'] + [target])

def makedirs(*targets):
    sh.call(['mkdir', '-p'] + list(targets))


@command.add_sub
@argument('region_bed', help='with header')
@argument('fasta', help='reference fasta')
@argument('-db', required=True, help='blast db')
@argument('--col-id', default='name', help='column name of region id')
@argument('-o', '--out-dir', default='refdup_analyze')
@argument('-j', '--njobs', type=int, default=8)
@argument('-f', '--force-rerun', action='store_true')
@argument('-R', '--rerun-from', type=int, default=float('inf'))
@argument('-m', '--min-size', type=int, default=500)
@argument('-M', '--max-size', type=int, default=300000)
def bed_refdup_analyze(args):
    """
    Outputs:
        {outdir}/region.fa
        {outdir}/region.ref.blastn
        {outdir}/region.summary.txt
    """
    step = 0
    out_dir = args.out_dir
    region_bed = '{0}/region.bed'.format(out_dir)
    split_dir = '{0}/region.split'.format(out_dir)
    region_summary = '{0}/region.summary.txt'.format(out_dir)

    sh.call(['mkdir', '-p', out_dir])

    logging.info('Rerun from step: %s', args.rerun_from)
    step += 1
    logging.info('step %s: create copy of input bed', step)
    if args.rerun_from <= step or not is_valid_file(region_bed):
        sh.call('''
        cat {0} | awk 'NR==1; NR>1 && $3-$2 >= {min_size} && $3-$2 <= {max_size}'
        '''.format(
            args.region_bed, min_size=args.min_size, max_size=args.max_size),
            stdout=region_bed)

    queue = UGEQueue()

    step += 1
    logging.info('step %s: clean up and split files', step)
    if args.rerun_from <= step or not is_valid_file('{dir}/region.1.bed'.format(dir=split_dir)):  # TODO use more strict check?
        require(region_bed)
        remove_all(split_dir)
        makedirs(split_dir)
        with open(region_bed) as fp:
            header = next(fp).rstrip()
            _sp = list(file_split(fp, nlines=10, prefix='region.', suffix='.bed', dir=split_dir, header=header))
            logging.info('%s was split into %s files', region_bed, len(_sp))


    step += 1
    logging.info('step %s: make region fasta and apply blastn', step)
    split_files = glob.glob('{dir}/region.*.bed'.format(dir=split_dir))
    if args.rerun_from <= step or not all(is_valid_file('{dir}/region.{task_id}.blastn'.format(dir=split_dir, task_id=i + 1)) for i in range(len(split_files))):
        for f in split_files:
            require(f)
        ntask = len(split_files)
        cmd = '''
            bedtools getfasta -name -bed <(csvcut -t -c chrom,start,end,{id} {bed} | csvformat -T | sed 1d) -fi {fasta} -fo {region_fasta}
            blastn -db {db} -num_threads {njobs} -perc_identity 95 -max_target_seqs 30 -best_hit_overhang 0.1 -query {region_fasta} -evalue 0.1 -outfmt '7 std gaps qlen slen qseq sseq' > {blastn_out}
        '''.format(
            id=args.col_id,
            bed='{dir}/region.$SGE_TASK_ID.bed'.format(dir=split_dir),
            fasta=args.fasta,
            region_fasta='{dir}/region.$SGE_TASK_ID.fa'.format(dir=split_dir),
            db=args.db,
            blastn_out='{dir}/region.$SGE_TASK_ID.blastn'.format(dir=split_dir),
            njobs=args.njobs)
        queue.array_call(cmd, name='region_blastn', memory=60, slot=args.njobs, ntask=ntask)

        #sh.call(""" rm -r {0}/task/* """.format(out_dir)
        #bedtools getfasta -name -bed <(csvcut -t -c chrom,start,end,{id} {bed} | csvformat -T | sed 1d) -fi {fasta} -fo {out}
        #""".format(id=args.col_id, bed=region_bed, fasta=1))

    step += 1
    logging.info('step {0}: filter blastn', step)
    if args.rerun_from <= step or not is_valid_file(region):
        _filter_blastn(region_blastn, region_summary)


#TODO make command
def _filter_blastn(blastn, summary):
    min_qlen_ratio = .9
    min_matched = 98
    with open(blastn) as reader: #, open(summary, 'w+') as writer:
        tab = BlastTabFile(reader)
        for query, tab in tab.iter_query_tables():
            aln_qlen = tab.qseq.map(lambda x: len(x.replace('-', '')))
            cond = (aln_qlen >= qlen * min_qlen_ratio) & (identity >= min_matched)
            print (query, tab[cond])
