from __future__ import absolute_import
from __future__ import print_function
import sys
import logging
from argtools import command, argument
#import vcf as VCF
import cyvcf2 as VCF
import pandas as pd
from io import open
from ..plot import PdfGrid
from ..utils import safediv
import numpy as np
from scipy import stats
from sklearn.linear_model import LinearRegression
from sklearn.feature_selection import f_regression


def read_sample_file(fname):
    """
    ID_1 ID_2 missing cov_1 cov_2 cov_3 cov_4 pheno1 bin1
    0 0 0 D D C C P B
    1 1 0.007 1 2 0.0019 -0.008 1.233 1
    2 2 0.009 1 2 0.0022 -0.001 6.234 0
    3 3 0.005 1 2 0.0025 0.0028 6.121 1
    4 4 0.007 2 1 0.0017 -0.011 3.234 1
    5 5 0.004 3 2 -0.012 0.0236 2.786 0
    """
    default_types = {'sex': 'D', 'age': 'C'}
    with open(fname) as fp:
        header = next(fp).rstrip().split(' ')
        types = next(fp).rstrip().split(' ')

        types = [default_types.get(col, t) for col, t in zip(header, types)]  # override with default types
        dtypes = {'ID_1': 'str', 'ID_2': 'str', 'missing': 'float'}
        for col, t in zip(header, types):
            if t == '0':
                continue
            elif t in ('C', 'P'):
                dtypes[col] = 'float'
            elif t in ('D', 'B'):
                dtypes[col] = 'category'

        tab = pd.read_table(fp, sep=' ', names=header, dtype=dtypes)
    return tab


def get_vcf_target_rec(vcf, id):
    reader = VCF.Reader(vcf)
    for rec in reader:
        if rec.ID == id:
            return reader, rec

def isnan(obj):
    return isinstance(obj, float) and np.isnan(obj)


@command.add_sub
@argument('sample')
@argument('-f', '--vcf', help='VCF or gzipped VCF')
@argument('-i', '--id')
@argument('-o', '--output', default='test.pdf')
def pheno_plot(args):
    tab = read_sample_file(args.sample)

    min_qual = 90
    col_pheno = tab.columns[-1]  # default target pheno

    import matplotlib as mp
    mp.use('Agg')
    import matplotlib.pyplot as plt
    import seaborn as sns

    def with_gt_dosage(tab, vcf, id):
        vcf_reader, rec = get_vcf_target_rec(vcf, id)
        if not rec:
            logging.error('Cannot find {} in vcf: {}'.format(id, vcf))
            raise Exception()

        samples = vcf_reader.samples
        if rec.INFO.get('SVTYPE') == 'CNV':
            dosage = rec.format('CN')[:, 0]
            quality = rec.format('CNQ')[:, 0]  # 0-100
            logging.info('Variant is CNV')
        else:
            raise NotImplementedError

        dosage_tab = pd.DataFrame({'id': samples, 'dosage': dosage, 'dosage_qual': quality})
        tab = tab.merge(dosage_tab, left_on='ID_1', right_on='id', how='left')
        return tab

    if args.vcf and args.id:
        tab = with_gt_dosage(tab, vcf=args.vcf, id=args.id)
    logging.info(tab)

    # TODO category phenotype
    def pheno_plot(tab):
        # SNP: 0, 1, 2
        # CN: 0, 1, 2, ...
        tab[col_pheno].hist()

    def pheno_sex_plot(tab):
        sns.violinplot(x=tab['sex'], y=tab[col_pheno])

    def pheno_sex_age_plot(tab):
        plt.scatter(x=tab['age'], y=tab[col_pheno], c=tab['sex'])

    def pheno_sex_dosage_plot(tab):
        # SNP: 0, 1, 2
        # CN: 0, 1, 2, ...
        #plt.scatter(x=tab['dosage'], y=tab[col_pheno], c=tab['sex'])
        sns.swarmplot(x=tab['dosage'], y=tab[col_pheno], hue=tab['sex'])

    logging.info('Plot to %s', args.output)
    plt.figure(figsize=(11.7, 8.3))
    with PdfGrid(args.output, row=2, col=2) as pg:
        plt.gcf().subplots_adjust()
        next(pg)
        pheno_plot(tab)
        next(pg)
        pheno_sex_plot(tab)
        next(pg)
        pheno_sex_age_plot(tab)
        next(pg)
        pheno_sex_dosage_plot(tab)



class VCFAssocCNV(object):
    _header = ['chrom', 'start', 'end', 'name', 'nsample', 'ncalled', 'call_rate', 'alt_sample', 'alt_rate', 'nassoc']

    def __init__(self, vcf, sample, min_qual=90, sample_id='ID_1', pheno_id=None, covariates=None, sex='both'):
        self.min_qual = min_qual
        sample_tab = read_sample_file(sample)
        if pheno_id is None:
            pheno_id = sample_tab.columns[-1]   # take last column as assoc phenotype
        self.pheno_id = pheno_id
        self.covariates = covariates or []

        self._reader = reader = VCF.Reader(vcf)
        vcf_samples = reader.samples

        # TODO restrict sex
        if sex == 'both':
            pass
        elif sex == 'male':
            sample_tab = sample_tab[sample_tab.sex == '1']
        elif sex == 'female':
            sample_tab = sample_tab[sample_tab.sex == '2']
        else:
            logging.warning('sex is one of 0/1/2')
            raise NotImplementedError

        assoc_tab = pd.DataFrame({'vcf_sample': vcf_samples})
        assoc_tab = assoc_tab.merge(sample_tab, left_on='vcf_sample', right_on=sample_id, how='left')
        self._assoc_tab = assoc_tab
        self.sample_ok = assoc_tab[sample_id].map(lambda x: not isnan(x))   # ID_1 is exist in sample file
        self.pheno_ok = assoc_tab[pheno_id].map(lambda x: not isnan(x))   # phenotype is not an NA
        if covariates:
            self.covs_ok = np.prod([assoc_tab[cov_id].map(lambda x: not isnan(x)) for cov_id in self.covariates], dtype=bool, axis=0)  # covariates are not an NA
        else:
            self.covs_ok = np.ones(len(assoc_tab), dtype=bool)
        #print (self.covs_ok)

    def iter_assocs(self):
        for rec in self._reader:
            arec = self._assoc(rec)
            if arec:
                yield arec

    def _get_data(self, rec):
        if rec.INFO.get('SVTYPE') != 'CNV':
            return

        data = {}
        data['chrom'] = rec.CHROM
        data['start'] = rec.POS - 1
        data['end'] = rec.INFO['END']
        data['name'] = rec.ID
        data['nsample'] = self.sample_ok.sum()

        dosage = pd.Series(rec.format('CN')[:, 0])
        quality = pd.Series(rec.format('CNQ')[:, 0])  # 0-100
        is_called = dosage.map(lambda x: not isnan(x))
        has_good_call = is_called & (quality >= self.min_qual)

        total_call = np.sum(self.sample_ok & is_called)
        valid_call = has_good_call & self.sample_ok
        valid_call_nsample = valid_call.sum()
        data['ncalled'] = valid_call_nsample

        alt_sample = (dosage[valid_call] != 2).sum()
        data['alt_sample'] = alt_sample
        data['alt_rate'] = safediv(1. * alt_sample, valid_call_nsample, error=0.)
        data['call_rate'] = safediv(1. * valid_call_nsample, data['nsample'], error=0.)

        is_assoc_target = has_good_call & self.sample_ok & self.pheno_ok & self.covs_ok
        data['a_dosage'] = a_dosage = dosage[is_assoc_target]
        data['a_covs'] = self._assoc_tab[self.covariates][is_assoc_target] if self.covariates else None
        data['a_pheno'] = a_pheno = self._assoc_tab[self.pheno_id][is_assoc_target]
        data['nassoc'] = len(a_dosage)

        return data


class VCFAssocCNVPearson(VCFAssocCNV):
    @property
    def header(self):
        return self._header + ['cor', 'pvalue']

    def _assoc(self, rec):
        data = self._get_data(rec)
        if data is None:
            return

        if data['nassoc']:
            #print (a_dosage, a_pheno)
            cor, pvalue = stats.pearsonr(data['a_dosage'], data['a_pheno'])
        else:
            cor, pvalue = float('nan'), float('nan')

        data['cor'] = cor
        data['pvalue'] = pvalue
        return data


def hierarchical_regression(X, y):
    import statsmodels.api as sm
    X = X.astype(float)
    X = sm.add_constant(X)
    results = []
    pvalues = []
    for icol in range(X.shape[1]):
        X1 = X[:, :icol + 1]
        res = sm.OLS(y, X1).fit()   # TODO normalize?
        if results:
            lr_test, p_value, df_diff = res.compare_lr_test(results[-1])
        else:
            p_value = res.pvalues[0]  # pvalue for intercept
        pvalues.append(p_value)
        results.append(res)

    return {'results': results, 'pvalues': pvalues}


class VCFAssocCNVLinear(VCFAssocCNV):
    @property
    def header(self):
        exts = []
        exts.append('b0')
        for cov_id in self.covariates:
            exts.append('cov_b_{0}'.format(cov_id))
            exts.append('cov_p_{0}'.format(cov_id))
        exts.append('beta')
        exts.append('pvalue')
        return self._header + exts

    def _assoc(self, rec):
        data = self._get_data(rec)
        if data is None:
            return

        y = data['a_pheno']
        if len(data['a_covs']):
            X = np.c_[data['a_covs'], data['a_dosage'].reshape((-1, 1))]
            #X = np.c_[data['a_dosage'].reshape((-1, 1)), data['a_covs']]
        else:
            X = np.c_[data['a_dosage']].reshape((-1, 1))

        #logging.info(X)
        lin_reg= LinearRegression(normalize=True)
        if len(X):
            res = hierarchical_regression(X, y.values)
            b0 = res['results'][-1].params[0]
            theta = res['results'][-1].params[1:]
            pvals = res['pvalues'][1:]  # ignore intercept
            # 
            # lin_reg.fit(X, y)
            # b0 = lin_reg.intercept_
            # theta = lin_reg.coef_
            # F, pvals = f_regression(X, y, center=True)
        else:
            b0 = float('nan')
            theta = [float('nan')] * (len(self.covariates) + 1)
            pvals = [float('nan')] * (len(self.covariates) + 1)

        data['b0'] = b0
        for i, cov_id in enumerate(self.covariates): #theta[1:]):
            data['cov_b_{0}'.format(cov_id)] = theta[i]
            data['cov_p_{0}'.format(cov_id)] = pvals[i]
        data['beta'] = theta[-1]
        data['pvalue'] = pvals[-1]
        return data


@command.add_sub
@argument('sample')
@argument('vcf', help='VCF or gzipped VCF')
@argument('-m', '--method', choices=['pearson', 'spearman', 'linear', 't', 'welch', 'logistic'])
@argument('-s', '--sex', choices=['both', 'male', 'female'], default='both')
@argument('-p', '--pheno-id')
@argument('-c', '--covs')
def vcf_assoc_cnv(args):
    """
    Assume autosomals

    Output TSV:
        chrom
        start
        end
        name
        (regino_type, strand, region_id, region_name) <= need gene annotation file

        # basic sample info
        call_rate : change with min_cnq
        alt_sample :
        alt_rate :
        cns :
        nsample

        # pearson (cor)
        cor
        pvalue

        # linear (standardized linear)
        b0
        cov_b_{cov_name1}
        cov_p_{cov_name1}
        cov_b_{cov_name2}
        cov_p_{cov_name2}
        beta                # genotype
        pvalue

        # logistic regression  # TODO
    """
    vcffile = '/dev/stdin' if args.vcf == '-' else args.vcf
    min_cnq = 90

    cls = {
          'pearson': VCFAssocCNVPearson,
          'linear': VCFAssocCNVLinear,
    }[args.method]

    covs = [] if not args.covs else args.covs.split(',')
    vcf_assoc = cls(args.vcf, args.sample, sex=args.sex, pheno_id=args.pheno_id, covariates=covs)

    header = vcf_assoc.header
    print (*header, sep='\t')
    for rec in vcf_assoc.iter_assocs():
        print (*[rec[col] for col in header], sep='\t')
