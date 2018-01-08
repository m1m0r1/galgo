#!/usr/bin/env python
from __future__ import print_function
from argtools import command, argument
import logging
import numpy as np
import itertools
from collections import namedtuple
from itertools import groupby
from multiprocess import Pool
from scipy.stats import norm
from ..covfit import load_gender_info, get_ploidies
from ..covfit import pmm
from ..utils import iter_tabs, Counter, safediv


@command.add_sub
@argument('coverage_bed')
@argument('--col-region-id')
@argument('--col-sample-id', default='sampleid')
@argument('--col-count', default='valid_count')
@argument('--col-exp-count', default='valid_exp_count')
@argument('-d', '--diploid-value', type=float, default=1.0)
@argument('-m', '--max-copy', default=10)
@argument('-la', type=float, default=0.3, help='penalty of offset')
@argument('-lb', type=float, default=1.0, help='penalty of ratio')
@argument('-a', type=float, default=0., help='penalty of offset')
@argument('-b', type=float, default=1.0, help='penalty of ratio')
@argument('-f', '--max-diff', type=float, default=1e-6, help='max diff of log likelihood')
@argument('-p', '--process', type=int, default=8)
@argument('--step-limit', type=int, help='the number of steps until inference stopped')
@argument('--gender-info', help='col1: sample_id, col2: {1:male, 2:female}')
def covfit_pmm(args):
    """

    1. chrom
    2. start
    3. end
    4. region id
    """
    sample_genders = load_gender_info(args.gender_info)
    zoom = 2. / args.diploid_value

    def run():
        param_keys = 'status top_cn top_theta a b c theta ia ib ic ll nsamples step'.split(' ')
        pool = Pool(args.process)
        #logging.info(sample_genders)

        task_inputs = gen_inputs()

        if args.process == 1:
            results = itertools.imap(task, task_inputs)
        else:
            results = pool.imap(task, task_inputs, 100)

        data_header = ['chrom', 'start', 'end', 'region_id', 'nsamples']
        param_header = ['status', 'step', 'ia', 'ib', 'a', 'b', 'top_cn', 'top_theta', 'theta', 'll']
        out_header = data_header + param_header
        print (*out_header, sep='\t')
        for data, params in results:
            row = [data[c] for c in data_header] \
                + [params[c] for c in param_header]
            print (*row, sep='\t')


    def gen_inputs():
        with open(args.coverage_bed) as fp:
            it = iter_tabs(fp)
            header = next(it)
            idx_region_id = header.index(args.col_region_id) if args.col_region_id else 3
            idx_sample_id = header.index(args.col_sample_id)
            idx_count = header.index(args.col_count)
            idx_exp_count = header.index(args.col_exp_count)

            for region_id, rows in groupby(it, lambda x: x[idx_region_id]):
                rows = list(rows)
                chrom = rows[0][0]

                samples = []
                genders = []
                counts = []
                exp_counts = []

                for row in rows:
                    sample_id = row[idx_sample_id]
                    count = float(row[idx_count])
                    exp_count = float(row[idx_exp_count]) / zoom
                    raw_cn = safediv(count, exp_count)
                    if np.isnan(raw_cn):
                        continue
                    if raw_cn > args.max_copy * 1.1:   # regard this as outlier
                        logging.warning('Omit outlier (data: %s)', zip(header, row))
                        continue
                    samples.append(sample_id)
                    genders.append(sample_genders[sample_id])
                    counts.append(count)
                    exp_counts.append(exp_count)

                ploidies = get_ploidies(chrom, genders)

                data = dict(
                        chrom = rows[0][0],
                        start = int(rows[0][1]),
                        end = int(rows[0][2]),
                        region_id = region_id,
                        counts = counts,
                        exp_counts = exp_counts,
                        la = args.la,
                        lb = args.lb,
                        a = args.a,
                        b = args.b,
                        max_diff = args.max_diff,
                        step_limit = args.step_limit,
                        ploidies = ploidies,
                        max_copy = args.max_copy,
                        nsamples = len(ploidies),
                )

                yield data

    def task(data):
        chrom = data['chrom']
        start = data['start']
        end = data['end']
        name = data['region_id']
        ploidies = data['ploidies']
        init_a = data['a']
        init_b = data['b']
        la = data['la']
        lb = data['lb']
        counts = np.array(data['counts'])
        exp_counts = np.array(data['exp_counts'])
        max_diff = data['max_diff']
        step_limit = data['step_limit']
        max_copy = data['max_copy']

        n_callable = len(counts)
        logging.info('[%s] n_callable: %s', name, n_callable)
        results = []
        norm = None
        try:
            assert n_callable > 0, 'The number of callable samples ({0}) is less than the minimum value'.format(n_callable)
            norm = pmm.Normalizer(counts=counts, exp_counts=exp_counts, ploidies=ploidies, max_copy=max_copy, a=init_a, b=init_b, la=la, lb=lb, max_diff=max_diff, step_limit=step_limit, name=name)
            result = norm.iteration()
            top_cn = np.argmax(result.theta)
            params = {
                    'status': 1,
                    'ia': init_a,
                    'ib': init_b,
                    'a': result.a,
                    'b': result.b,
                    'theta': ','.join(map(str, result.theta)),
                    'step': result.step,
                    'top_cn': top_cn,
                    'top_theta': result.theta[top_cn],
                    'll': result.ll,
            }
        except AssertionError as e:
            logging.warning('Assert error: %s', e)
            logging.warning('Cannot optimize model for %s', data)
            step = 0 if norm is None else norm.step
            params = {
                    'status': 0,
                    'ia': init_a,
                    'ib': init_b,
                    'a': init_a,
                    'b': init_b,
                    'theta': '',
                    'step': step,
                    'top_cn': 0,
                    'top_theta': 0,
                    'll': 0,
            }
        results.append(params)

        logging.info('Emitting %s:%s-%s, %s, %s', chrom, start, end, name, params)
        return (data, params)

    run()


@command.add_sub
@argument('fit_result')
@argument('--min-call-nsample', type=int, default=1)
@argument('-d', '--diploid-value', type=float, default=1.0)
@argument('-t', '--threshold', type=float, default=.8)
@argument('--sex-info', help='list file of 1:male, 2:female')
def covfit_pmm_genotype(args):
    """
    """
    nsample = args.nsample
    zoom = 2. / args.diploid_value
    gender_info = [int(line.rstrip()) for line in open(gender_info)]

    for row in iter_tabs(open(args.fit_result)):
        chrom = row[0]
        start = int(row[1])
        end = int(row[2])
        name = row[3]

        depths = np.array(row[-nsample:], dtype=float)
        callable_idxs = np.where(~ np.isnan(depths))[0]
        n_callable = len(callable_idxs)

        mod_depths = (depths * zoom)[callable_idxs]

        #mod_depths = depths * zoom
        keys = row[-args.nsample-2].split(':')
        params = dict(zip(keys, row[-args.nsample-1].split(':')))
        a = float(params['a'])
        b = float(params['b'])
        c = float(params['c'])
        try:
            assert n_callable >= args.min_call_nsample, 'The number of callable samples ({0}) is less than the minimum value'.format(n_callable)
            theta = np.array(map(float, params['theta'].split(',')))
            max_copy = len(theta) - 1

            ploidies = tuple(np.array(get_ploidies(chrom))[callable_idxs])

            #ploidies = get_ploidies(chrom)
            #ploidies = np.array(ploidies)

            norm = em.Normalizer(mod_depths, ploidies=ploidies, max_copy=max_copy, a=a, b=b, c=c, theta=theta, name=name)
            gts = list(norm.iter_genotypes())
            gts1 = []   # genotyped copy numbers
            gts2 = []   # genotype probability of 0/1/2
            for prob in gts:
                max_cn = np.argmax(prob)
                if prob[max_cn] >= args.threshold:
                    gts1.append(max_cn)
                else:
                    gts1.append('.')
                gts2.append('{0:.1f}:{1:.1f}:{2:.1f}'.format(prob[0], prob[1], prob[2]))
        except (ValueError, AssertionError):
            gts1 = ['.'] * nsample
            gts2 = ['.:.:.'] * nsample

        out_row = list(row[:-args.nsample-2])
        out_row.append(a)
        out_row.append(b)
        out_row.append(c)
        out_row.append(params['theta'])
        out_row.append(','.join(map(str, gts1)))
        out_row.append(','.join(map(str, gts2)))
        print (*out_row, sep='\t')
