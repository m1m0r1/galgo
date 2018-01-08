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
from ..utils import iter_tabs, Counter, with_header, attrdict


_InputRow = namedtuple('_InputRow', 'chrom start end region_id samples covs valid_ratio valid_len'.split(' '))

def gen_covfit_inputs(args):
    if args.gender_info:
        sample_genders = load_gender_info(args.gender_info)
    else:
        sample_genders = defaultdict(lambda : 1)
    zoom = 2. / args.diploid_value
    max_copy = getattr(args, 'max_copy', None)

    if args.format == 'vertical':
        def iter_inputs(it):
            for region_id, rows in groupby(it, lambda x: x[idx_region_id]):
                rows = list(rows)
                chrom = rows[0][0]
                start = int(rows[0][1]),
                end = int(rows[0][2]),
                samples = np.array([row[idx_sample_id] for row in rows])
                covs = np.array([row[idx_cov] for row in rows], dtype=float)
                valid_ratio = np.array([row[idx_valid_ratio] for row in rows], dtype=float)
                valid_len = np.array([row[idx_valid_len] for row in rows], dtype=int)
                yield _InputRow(chrom=chrom, start=start, end=end, region_id=region_id, samples=samples, covs=covs, valid_ratio=valid_ratio, valid_len=valid_len)

    elif args.format == 'horizontal':
        def iter_inputs(it):
            for row in it:
                chrom = row[0]
                start = int(row[1])
                end = int(row[2])
                region_id = row[idx_region_id]
                samples = np.array(row[idx_sample_id].split(','))
                covs = np.array(row[idx_cov].split(','), dtype=float)
                valid_ratio = np.array(row[idx_valid_ratio].split(','), dtype=float)
                valid_len = np.array(row[idx_valid_len].split(','), dtype=int)
                yield _InputRow(chrom=chrom, start=start, end=end, region_id=region_id, samples=samples, covs=covs, valid_ratio=valid_ratio, valid_len=valid_len)

    else:
        raise NotImplementedError(args.format)

    with open(args.coverage_bed) as fp:
        it = iter_tabs(fp)
        header = next(it)
        idx_region_id = header.index(args.col_region_id) if args.col_region_id else 3
        idx_sample_id = header.index(args.col_sample_id)
        idx_cov = header.index(args.col_cov)
        idx_valid_ratio = header.index('valid_ratio')
        idx_valid_len = header.index('valid_len')

        for row in iter_inputs(it):
            covs = row.covs * zoom # zoom

            # filter condition
            sel = (~np.isnan(covs) & (covs <= args.max_copy))
            if getattr(args, 'min_valid_ratio', None):
                sel &= row.valid_ratio >= args.min_valid_ratio
            if getattr(args, 'min_valid_len', None):
                sel &= row.valid_len >= args.min_valid_len

            mod_covs = covs[sel]
            mod_samples = row.samples[sel]

            genders = [sample_genders[s] for s in mod_samples]
            ploidies = get_ploidies(row.chrom, genders)
            ploidy_counts = Counter(ploidies)
            data = attrdict(
                    row = row,
                    mod_covs = mod_covs,
                    mod_samples = mod_samples,
                    args = args,
                    ploidies = ploidies,
                    nsamples = ','.join([str(ploidy_counts[cn]) for cn in (0, 1, 2)]),
            )
            yield data


@command.add_sub
@argument('coverage_bed')
@argument('--col-region-id')
@argument('--col-sample-id', default='sampleid')
@argument('--col-cov', default='mean_cov')
@argument('-F', '--format', default='vertical', choices=['vertical', 'horizontal'])
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
#TODO maybe min_valid_ratio and min_valid_len should be add to the options
def covfit_gmm(args):
    """

    1. chrom
    2. start
    3. end
    4. region id
    """
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
            task_inputs = gen_covfit_inputs(args)

            if args.process == 1:
                results = itertools.imap(task, task_inputs)
            else:
                results = pool.imap(task, task_inputs, 1)
                #results = imap_wrap(pool, task, task_inputs, 100)

            row_header = ['chrom', 'start', 'end', 'region_id']
            data_header = ['nsamples']
            param_header = ['status', 'ia', 'ib', 'ic', 'a', 'b', 'c', 'top_cn', 'top_theta', 'theta', 'll', 'step']
            out_header = row_header + data_header + param_header
            print (*out_header, sep='\t')
            for data, params in results:
                row = [getattr(data.row, c) for c in row_header] \
                    + [getattr(data, c) for c in data_header] \
                    + [params[c] for c in param_header]
                print (*row, sep='\t')
            pool.close()
        except Exception as e:
            pool.terminate()
            raise e
        finally:
            pool.join()

    def task(data):
        chrom = data.row.chrom
        start = data.row.start
        end = data.row.end
        name = data.row.region_id
        covs = data.mod_covs
        la = data.args.la
        lb = data.args.lb
        lc = data.args.lc
        init_as = data.args.a
        init_bs = data.args.b
        init_abs = [(a, b) for a in init_as for b in init_bs]
        init_c = data.args.c
        delta = data.args.delta
        max_diff = data.args.max_diff
        step_limit = data.args.step_limit
        max_copy = data.args.max_copy
        ploidies = data.ploidies

        #ploidy_counts = Counter(ploidies)  # 0, 1, 2
        #nsamples_str =  ','.join(str(ploidy_counts[i]) for i in (0, 1, 2))

        n_callable = len(covs)
        logging.info('[%s] n_callable: %s', name, n_callable)
        logging.debug('%s', data)
        #init_as = (0.,)
        results = []
        for init_a, init_b in init_abs:
            init_params = {
                    'status': 0,
                    'ia': init_a,
                    'ib': init_b,
                    'ic': init_c,
                    'a': init_a,
                    'b': init_b,
                    'c': init_c,
                    'theta': '',
                    'step': 0,
                    'top_cn': 0,
                    'top_theta': 0,
                    'll': float('-inf'),
            }
            norm = None
            try:
                assert n_callable > 0, 'The number of callable samples ({0}) is less than the minimum value'.format(n_callable)
                norm = gmm.Normalizer(covs, ploidies=ploidies, max_copy=max_copy, a=init_a, b=init_b, c=init_c, la=la, lb=lb, lc=lc, max_diff=max_diff, step_limit=step_limit, name=name)
                #init_a, init_b, _, _ = get_best_fit(mod_depths)
                result = norm.iteration()
                top_cn = np.argmax(result.theta)
                params = dict(init_params,
                        status=1, a=result.a, b=result.b, c=result.c,
                        theta = ','.join(map(str, result.theta)),
                        step = result.step,
                        top_cn = top_cn,
                        top_theta = result.theta[top_cn],
                        ll = result.ll,
                )
            except AssertionError as e:
                logging.warning('Assert error: %s', e)
                logging.warning('Cannot optimize model for %s', data)
                step = 0 if norm is None else norm.step
                params = dict(init_params, step=step)
            results.append(params)
        best_params = max(results, key=lambda v: (v['status'], v['ll']))

        logging.info('Emitting %s:%s-%s, %s, %s', chrom, start, end, name, best_params)
        return (data, best_params)

    run()


@command.add_sub
@argument('coverage_bed')
@argument('-p', '--fit-params', required=True)
@argument('-F', '--format', default='vertical', choices=['vertical', 'horizontal'])
@argument('--min-valid-ratio', type=float, default=.5)
@argument('--min-valid-len', type=int, default=100)
@argument('-t', '--threshold', type=float, default=.8)
@argument('--col-sample-id', default='sampleid')
@argument('--col-region-id')
@argument('--col-cov', default='mean_cov', help='column name of normalized coverage')
@argument('--delta', type=float, default=0.3, help='variance for zero copy')    # TODO incorporate into the result
@argument('-d', '--diploid-value', type=float, default=1.0)
@argument('-m', '--max-copy', type=int, default=10)
@argument('--gender-info', help='col1: sample_id, col2: {1:male, 2:female}')
def covfit_gmm_genotype(args):
    """
    gt        # filtered by specified condition
    map_gt    #
    likelihoods
    """
    fit_it = with_header(iter_tabs(open(args.fit_params)), use_header=True)
    fit_header = next(fit_it)
    #col_region_id = args.col_region_id if args.col_region_id else 'name'
    cov_header = 'chrom start end region_id'.split(' ')
    header = cov_header + ['samples', 'copy_number', 'gts', 'gt_probs']
    print (*header, sep='\t')
    covfit_inputs = gen_covfit_inputs(args)

    NA = float('nan')
    for data, params in zip(covfit_inputs, fit_it):
        chrom = data.row.chrom
        start = data.row.start
        end = data.row.end
        name = data.row.region_id
        assert (chrom, start, end, name) == (params['chrom'], int(params['start']), int(params['end']), params['region_id']), 'chrom, start, and end should be aligned for the two input files'
        call_covs = data.mod_covs
        call_samples = data.mod_samples
        ploidies = data.ploidies
        chrom = params['chrom']
        start = int(params['start'])
        end = int(params['end'])

        try:
            a = float(params['a'])
            b = float(params['b'])
            c = float(params['c'])
            theta = np.array(map(float, params['theta'].split(',')))
        except ValueError:
            theta = NA

        sample_gt_keys = {}
        sample_gt_probs = {}
        sample_top_cns = {}
        if sum(ploidies) > 0 and not np.any(np.isnan(theta)):
            norm = gmm.Normalizer(call_covs, ploidies=ploidies, max_copy=args.max_copy, a=a, b=b, c=c, theta=theta, name=name, delta=args.delta)
            cn_probs = np.array(list(norm.iter_cns()))
            top_cns = np.argmax(cn_probs, axis=1)
            top_cn_probs = np.array([cn_probs[i, c] for i, c in enumerate(top_cns)])
            logging.info(cn_probs)
            logging.info(top_cns)
            logging.info(top_cn_probs)
            called_sel = top_cn_probs >= args.threshold
            sample_top_cns.update(zip(call_samples[called_sel], top_cns[called_sel]))
            for sample, gts in zip(call_samples, norm.iter_genotypes()):
                gt_keys, gt_vals = get_allelotype_records(gts)
                sample_gt_keys[sample] = gt_keys
                sample_gt_probs[sample] = gt_vals

        top_cns = [str(sample_top_cns.get(s, '.')) for s in data.row.samples]
        gt_keys = [';'.join(map(str, sample_gt_keys.get(s, ()))) or '.' for s in data.row.samples]
        gt_probs = [';'.join(map(str, sample_gt_probs.get(s, ()))) or '.' for s in data.row.samples]

        out_row = list(getattr(data.row, col) for col in cov_header)
        out_row.append(','.join(data.row.samples))
        out_row.append(','.join(top_cns))
        out_row.append(','.join(gt_keys))
        out_row.append(','.join(gt_probs))
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
    header = list(cov_header) + ['copy_number', 'top_cn', 'top_gt', 'top_gt_prob', 'samples', 'gts', 'cn_probs']
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
        gt_keys = ['.']
        gt_vals = [1.]
        if ploidy == 0:
            top_cn = 0
            cn = 0
            gt_keys = ['0']
            gt_vals = [1.]
        elif np.any(np.isnan(theta)):
            max_copy = len(theta) - 1
            norm = gmm.Normalizer(np.array([mod_depth]), ploidies=[ploidy], max_copy=max_copy, a=a, b=b, c=c, theta=theta, name=name, delta=args.delta)
            cn_probs = np.array(list(norm.iter_cns()))
            top_cns = np.argmax(cn_probs, axis=1)
            cn = top_cn if cn_probs[top_cn] >= args.threshold else -1

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
