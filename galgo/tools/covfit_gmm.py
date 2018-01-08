#!/usr/bin/env python
from __future__ import print_function
from argtools import command, argument
import logging
import numpy as np
import itertools
from collections import namedtuple, defaultdict
from itertools import groupby
from builtins import zip
from multiprocess import Pool
import signal
from scipy.stats import norm
from ..covfit import load_gender_info, get_ploidies
#from ..covfit import gmm
from ..covfit import gmm_v2 as gmm
from ..utils import iter_tabs, Counter, with_header


@command.add_sub
@argument('coverage_bed')
@argument('--col-region-id')
@argument('--col-sample-id', default='sampleid')
@argument('--col-cov', default='mean_cov')
@argument('-d', '--diploid-value', type=float, default=1.0)
@argument('-m', '--max-copy', type=int, default=10)
@argument('-la', type=float, default=0.3, help='penalty of offset')
@argument('-lb', type=float, default=1.0, help='penalty of ratio')
@argument('-lc', type=float, default=0.3, help='penalty of precision')
@argument('-a', nargs='*', type=float, default=[0., 0.4, 0.8], help='penalty of offset')
@argument('-b', nargs='*', type=float, default=[1.0], help='penalty of ratio')
@argument('-c', type=float, default=1.0, help='penalty of precision')
@argument('-f', '--max-diff', type=float, default=1e-6, help='max diff of log likelihood')
@argument('-p', '--process', type=int, default=8)
@argument('--delta', type=float, default=0.3, help='variance for zero copy')
@argument('--step-limit', type=int, help='the number of steps until inference stopped')
@argument('--gender-info', help='col1: sample_id, col2: {1:male, 2:female}')
def covfit_gmm(args):
    """

    1. chrom
    2. start
    3. end
    4. region id
    """
    if args.gender_info:
        sample_genders = load_gender_info(args.gender_info)
    else:
        sample_genders = defaultdict(lambda : 1)
    zoom = 2. / args.diploid_value

    # not required
    def imap_wrap(pool, task, task_inputs, chunksize=1):  # test
        result = pool.imap(task, task_inputs, chunksize=1)
        while 1:
            yield result.next(timeout=999999)


    def run():
        def init_pool():
            logging.info('Init pool')
            signal.signal(signal.SIGINT, signal.SIG_IGN)

        pool = Pool(args.process, init_pool)
        try:
            task_inputs = gen_inputs()

            if args.process == 1:
                results = itertools.imap(task, task_inputs)
            else:
                results = pool.imap(task, task_inputs, 1)
                #results = imap_wrap(pool, task, task_inputs, 100)

            data_header = ['chrom', 'start', 'end', 'region_id', 'nsamples']
            param_header = ['status', 'ia', 'ib', 'ic', 'a', 'b', 'c', 'top_cn', 'top_theta', 'theta', 'll', 'step']
            out_header = data_header + param_header
            print (*out_header, sep='\t')
            for data, params in results:
                row = [data[c] for c in data_header] \
                    + [params[c] for c in param_header]
                print (*row, sep='\t')
            pool.close()
        except Exception as e:
            pool.terminate()
            raise e
        finally:
            pool.join()

    def gen_inputs():
        with open(args.coverage_bed) as fp:
            it = iter_tabs(fp)
            header = next(it)
            idx_region_id = header.index(args.col_region_id) if args.col_region_id else 3
            idx_sample_id = header.index(args.col_sample_id)
            idx_cov = header.index(args.col_cov)

            for region_id, rows in groupby(it, lambda x: x[idx_region_id]):
                rows = list(rows)
                chrom = rows[0][0]

                samples = []
                genders = []
                covs = []

                for row in rows:
                    sample_id = row[idx_sample_id]
                    cov = float(row[idx_cov])
                    if np.isnan(cov):
                        continue
                    if cov * zoom > args.max_copy:  # coverage larger than the max copy is not used for parameter estimation (adhoc treatment)
                        continue
                    samples.append(sample_id)
                    genders.append(sample_genders[sample_id])
                    covs.append(cov)

                ploidies = get_ploidies(chrom, genders)
                ploidy_counts = Counter(ploidies)

                data = dict(
                        chrom = rows[0][0],
                        start = int(rows[0][1]),
                        end = int(rows[0][2]),
                        region_id = region_id,
                        covs = covs,
                        zoom = zoom,
                        la = args.la,
                        lb = args.lb,
                        lc = args.lc,
                        a = args.a,
                        b = args.b,
                        c = args.c,
                        delta = args.delta,
                        max_diff = args.max_diff,
                        step_limit = args.step_limit,
                        ploidies = ploidies,
                        max_copy = args.max_copy,
                        nsamples = ','.join([str(ploidy_counts[cn]) for cn in (0, 1, 2)]),
                )

                yield data

    def task(data):
        #signal.signal(signal.SIGINT, signal.SIG_IGN)
        chrom = data['chrom']
        start = data['start']
        end = data['end']
        name = data['region_id']
        ploidies = data['ploidies']
        init_as = data['a']
        init_bs = data['b']
        init_abs = [(a, b) for a in init_as for b in init_bs]
        init_c = data['c']
        la = data['la']
        lb = data['lb']
        lc = data['lc']
        covs = np.array(data['covs'])
        zoom = data['zoom']
        max_diff = data['max_diff']
        step_limit = data['step_limit']
        max_copy = data['max_copy']

        mod_depths = (covs * zoom)
        #ploidy_counts = Counter(ploidies)  # 0, 1, 2
        #nsamples_str =  ','.join(str(ploidy_counts[i]) for i in (0, 1, 2))

        n_callable = len(covs)
        logging.info('[%s] n_callable: %s', name, n_callable)
        #init_as = (0.,)
        results = []
        for init_a, init_b in init_abs:
            norm = None
            try:
                assert n_callable > 0, 'The number of callable samples ({0}) is less than the minimum value'.format(n_callable)
                norm = gmm.Normalizer(mod_depths, ploidies=ploidies, max_copy=max_copy, a=init_a, b=init_b, c=init_c, la=la, lb=lb, lc=lc, max_diff=max_diff, step_limit=step_limit, name=name)
                #init_a, init_b, _, _ = get_best_fit(mod_depths)
                result = norm.iteration()
                top_cn = np.argmax(result.theta)
                params = {
                        'status': 1,
                        'ia': init_a,
                        'ib': init_b,
                        'ic': init_c,
                        'a': result.a,
                        'b': result.b,
                        'c': result.c,
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
                        'ic': init_c,
                        'a': init_a,
                        'b': init_b,
                        'c': init_c,
                        'theta': '',
                        'step': step,
                        'top_cn': 0,
                        'top_theta': 0,
                        'll': 0,
                }
            results.append(params)
        best_params = max(results, key=lambda v: (v['status'], v['ll']))

        logging.info('Emitting %s:%s-%s, %s, %s', chrom, start, end, name, best_params)
        return (data, best_params)

    run()


@command.add_sub
@argument('fit_result')
@argument('--min-call-nsample', type=int, default=1)
@argument('-d', '--diploid-value', type=float, default=1.0)
@argument('--delta', type=float, default=0.3, help='variance for zero copy')    # TODO incorporate into the result
@argument('-t', '--threshold', type=float, default=.8)
@argument('--sex-info', help='list file of 1:male, 2:female')
def covfit_gmm_genotype(args):
    """
    """
    nsample = args.nsample
    zoom = 2. / args.diploid_value
    gender_info = [int(line.rstrip()) for line in open(gender_info)]
    get_ploidies = get_ploidy_factory(nsample, gender_info)

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

            norm = gmm.Normalizer(mod_depths, ploidies=ploidies, max_copy=max_copy, a=a, b=b, c=c, theta=theta, name=name, delta=delta)
            gts = list(norm.iter_cns())
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


@command.add_sub
@argument('coverage')
@argument('-p', '--fit-params', required=True)
@argument('--min-valid-ratio', type=float, default=.5)
@argument('--min-valid-len', type=int, default=100)
#@argument('--col-region-id')
@argument('--col-cov', default='mean_cov', help='column name of normalized coverage')
@argument('--delta', type=float, default=0.3, help='variance for zero copy')    # TODO incorporate into the result
@argument('-d', '--diploid-value', type=float, default=1.0)
@argument('-t', '--threshold', type=float, default=.8)
@argument('--sex', type=int, choices=[1, 2], required=True)
def covfit_gmm_genotype1(args):
    """
    gt        # filtered by specified condition
    map_gt    # 
    likelihoods
    """
    zoom = 2. / args.diploid_value
    cov_it = with_header(iter_tabs(open(args.coverage)), use_header=True)
    fit_it = with_header(iter_tabs(open(args.fit_params)), use_header=True)
    cov_header = next(cov_it)
    fit_header = next(fit_it)
    #col_region_id = args.col_region_id if args.col_region_id else 'name'
    header = list(cov_header) + ['copy_number', 'top_cn', 'top_gt', 'top_gt_prob', 'gts', 'cn_prob0', 'cn_prob1', 'cn_prob2', 'cn_probs', 'cn_probs']
    print (*header, sep='\t')

    NA = float('nan')
    for cov_row, params in zip(cov_it, fit_it):
        assert all(cov_row[col] == params[col] for col in ('chrom', 'start', 'end')), 'chrom, start, and end should be aligned for the two input files'
        chrom = params['chrom']
        start = int(params['start'])
        end = int(params['end'])
        name = params['region_id']

        depth = float(cov_row[args.col_cov])
        mod_depth = depth * zoom
        ploidy = get_ploidies(chrom, [args.sex])[0]

        try:
            a = float(params['a'])
            b = float(params['b'])
            c = float(params['c'])
            theta = np.array(map(float, params['theta'].split(',')))
        except ValueError:
            theta = NA

        cn_probs = (NA,)
        low_cn_probs = (NA,) * 3
        gt_keys = ['.']
        gt_vals = [1.]
        if ploidy == 0:
            top_cn = 0
            cn = 0
            gt_keys = ['0']
            gt_vals = [1.]
            low_cn_probs = (1., 0., 0.)
        elif not np.isnan(depth) and not np.any(np.isnan(theta)):
            max_copy = len(theta) - 1
            norm = gmm.Normalizer(np.array([mod_depth]), ploidies=[ploidy], max_copy=max_copy, a=a, b=b, c=c, theta=theta, name=name, delta=args.delta)
            cn_probs = list(norm.iter_cns())[0]
            top_cn = np.argmax(cn_probs)
            cn = top_cn if cn_probs[top_cn] >= args.threshold else -1

            low_cn_probs = cn_probs[0], cn_probs[1], cn_probs[2]

            valid_ratio = float(cov_row['valid_ratio'])
            valid_len = int(cov_row['valid_len'])
            if valid_ratio < args.min_valid_ratio \
               or valid_len < args.min_valid_len:
                cn = -1


            gts = list(norm.iter_genotypes())[0]
            gt_keys, gt_vals = get_allelotype_records(gts)
            # get at records
        else:
            top_cn = -1
            cn = -1

        out_row = list(cov_row[col] for col in cov_header)
        out_row.append(cn)
        out_row.append(top_cn)
        out_row.append(gt_keys[0])
        out_row.append('{0:.4f}'.format(gt_vals[0]))
        out_row.append(':'.join(gt_keys))
        out_row.append(':'.join(map('{0:.4f}'.format, gt_vals)))
        out_row.extend(tuple(map('{0:.4f}'.format, low_cn_probs)))
        out_row.append(':'.join(map('{0:.4f}'.format, cn_probs)))
        print (*out_row, sep='\t')


def get_allelotype_records(gts, psum_max=.99):
    """
    gts: {'0': 1.0}
         {'1/0': 0.8, '2/0': 0.1}
    """
    gt_keys = []
    gt_vals = []
    psum = .0
    for key, p in sorted(gts.items(), key=lambda x: x[1], reverse=True):
        gt_keys.append(key)
        gt_vals.append(p)
        psum += p
        if psum > psum_max:
            break
    return gt_keys, gt_vals
