from __future__ import print_function
import sys
from argtools import command, argument
from pysam import Samfile
import numpy as np
import logging
from operator import attrgetter
from itertools import islice, groupby
from collections import defaultdict
from builtins import range, filter, map
from ..utils import unzip, Counter
from ..utils import format as fmt
from .. import sh
from .sam_loader import SamDepthAnalyzer, SamPosInfo, SamModelBuilder
from .model import InferModel
from tqdm import tqdm


import networkx as nx
class VariantGraph(object):
    def __init__(self, refname):
        self.refname = refname
        self._graph = nx.MultiGraph()
        self._aln_vars = defaultdict(list)   # {aln_rec: [variant]}

    @classmethod
    def from_depth_analyzer(cls, depth_analyzer):
        refs = depth_analyzer.iter_ref_records()
        for refname, it in refs:
            self = cls(refname)
            logging.info('Storing variants')
            it = tqdm(it)
            #it = islice(it, 20000)
            for pos_info in it:
                if pos_info.nubases > 1 and list(sorted(pos_info.unsegs, reverse=True))[1] > 1:   #  at least two second major alleles were needed
                    self.add_variant(pos_info)
            logging.info('Storing variants is done')
            yield self

    def add_variant(self, pos_info):
        self._graph.add_node(pos_info.pos, **{'pos_info': pos_info})
        for info in pos_info.get_read_info():
            pcol = info.pop('pcol')
            data = dict(info, pos_id=pos_info.pos)
            aln = pcol.alignment
            self._aln_vars[aln].append(data)

            # edge
            if len(self._aln_vars[aln]) > 1:
                p1, p2 = self._aln_vars[aln][-2:]
                p1_id = p1['pos_id']
                p2_id = p2['pos_id']
                data = {
                    'qname': aln.qname,
                    'bases': [p1['bases'], p2['bases']],
                    'quals': [p1['quals'], p2['quals']],
                }
                self._graph.add_edge(p1_id, p2_id, **data)

    def get_graph(self):
        return self._graph

    def get_subgraphs(self):
        return list(nx.connected_component_subgraphs(self._graph))


@command.add_sub
@argument('bam')
def mc_pos_info(args):
    cols = 'contig pos nsegs nubases ubases unsegs uquals_mean'.split(' ')
    cols = 'contig pos nsegs nubases ubases unsegs wsegs uwsegs uquals_mean'.split(' ')

    getvals = {col: attrgetter(col) for col in cols}
    getvals['nsegs'] = lambda rec: rec.nsegs
    getvals['ubases'] = lambda rec: ','.join(rec.ubases)
    getvals['unsegs'] = lambda rec: ','.join(map(str, rec.unsegs))
    getvals['uwsegs'] = lambda rec: ','.join(map('{:.02f}'.format, rec.uwsegs))
    getvals['uquals_mean'] = lambda rec: ','.join(map('{:.02f}'.format, rec.uquals_mean))
    getvals['contig'] = lambda rec: tid_refnames[rec.tid]

    print (*cols, sep='\t')  # header
    with Samfile(args.bam) as sam:
        tid_refnames = {sam.get_tid(name): name for name in sam.references}
        da = SamDepthAnalyzer(sam)
        for rec in da.iter_records():
            vals = [getvals[col](rec) for col in cols]
            print (*vals, sep='\t')


import gzip
from dill import dill
@command.add_sub
@argument('bam')
@argument('-o', '--outdir', default='mc_vg_output')
def mc_variant_graph(args):
    sh.call(['mkdir', '-p', args.outdir])
    with Samfile(args.bam) as sam:
        da = SamDepthAnalyzer(sam)
        for vg in VariantGraph.from_depth_analyzer(da):
            outname = fmt('{args.outdir}/{vg.refname}.vg')
            logging.info('Saving graph to %s', outname)
            with open(outname, 'w+') as fo:
                dill.dump(vg, fo)
            subs = vg.get_subgraphs()
            graph = vg.get_graph()
            nnodes = len(graph.nodes())
            logging.info('Number of nodes for %s : %s', vg.refname, nnodes)
            logging.info('Number of subgraphs for %s : %s', vg.refname, len(subs))


def show_model(model):
    for ref in model.refs:
        logging.info('Reference: %s', ref.name)
        logging.info('Number of variants: %s', len(ref.variants))
        logging.info('Number of candidates: %s', len(ref.candidates))
        logging.info('Number of paired candidates: %s', len([c for c in ref.candidates if c.is_paired]))

    ref_counts = Counter()
    for frag in tqdm(model.fragments):
        key = frozenset([c.ref.name for c in frag.candidates])
        ref_counts[key] += 1

    # for frag in tqdm(model.fragments):
    #     key = [(c.ref.name, tuple(sorted(c.bases.items()))) for c in frag.candidates]
    #     diff = len(key) - len(set(key))
    #     if len(set(key)) > 1:
    #         logging.info('Number of redundant candidates is %s of %s (%s)', diff, len(key), set(key))

    for key in sorted(ref_counts, key=lambda key: (len(key), tuple(sorted(key)))):
        key1 = tuple(sorted(key))
        count = ref_counts[key]
        logging.info('Number of fragments for %s is %s', ','.join(key1), count)


def show_initial_genotype_info(infer_model):
    v_poss = {}  # {(rname, vpos): phase_start_pos}
    for ref in infer_model.refs:
        rname = ref.name
        blocks = infer_model.get_het_blocks(ref.name)
        for bl in blocks:
            for v in bl['variants']:
                v_poss[rname, v.pos] = bl['pos']

    for v in infer_model.variants:
        rname = v.ref.name
        pos = v.pos
        #depths = infer_model.get_full_depths(v)
        tot_depth = infer_model.get_total_depth(v)
        sinfo = infer_model.get_site_info(v)
        gen_scores = infer_model.get_init_gen_scores(v) # [((b1, b2, ..), prob, score)]

        v_pos = v_poss.get((rname, v.pos))

        print (rname, pos,
               v_pos,
               sinfo['ncands'],
               sinfo['site_score'],
               sinfo['base_counts'],
               sinfo['mean_depths'],
               gen_scores,
               sep='\t')


def show_variant_info(infer_model):
    for v in infer_model.variants:
        ref = v.ref.name
        pos = v.pos
        sinfo = infer_model.get_site_info(v)
        gen_probs = infer_model.get_variant_probs(v)   # [((b1, b2, ..), prob)]

        gens = []
        hap_probs = []
        for gp in gen_probs:
            probs = gp['prob']
            b = ','.join(gp['bases'])
            gens.append(b)
            #hap_probs.append(','.join(map('{0:.3f}'.format, probs)))
            hap_probs.append(','.join(map('{0:.3f}'.format, probs)))

        #print (gens, probs)
        print (ref, pos, sinfo,
               '|'.join(gens),
               '|'.join(hap_probs),
                sep='\t')

@command.add_sub
@argument('bam')
@argument('-o', '--outname', default='mc_vg_output.vg')
@argument('-r', '--regions', nargs='*')
@argument('--hap-depth', type=float, default=20.)
def mc_vg_create(args):
    import sys
    sys.setrecursionlimit(10000)
    with Samfile(args.bam) as sam:
        smb = SamModelBuilder(sam, regions=args.regions)
        show_model(smb.model)
    im = InferModel(smb.model, args.hap_depth)
    show_initial_genotype_info(im)
    #im.fix_
    show_variant_info(im)   # show variant, alts, depth, probs.

        # with gzip.open(args.outname, 'wb+') as fo:
        #     pickle.dump(smb.model, fo)


#TODO dill.load is segfault
import glob
@argument('vc_info')
@argument('-r', '--refnames', required=True, nargs='+')  # should change option in future ?
@argument('-f', '--fastas', required=True, nargs='+')
@argument('--subtype-map', help='tab delimited list of name to subtype map')
@argument('-t', '--threshold', type=int, default=5)
@argument('-s', '--score', default='edit', choices=['edit', 'match', 'hmm'])
@argument('--no-phase-probs', dest='use_phase_probs', action='store_false', default=True, help='only affected for hmm score')
def mc_ref_match(args):
    """
    Result info include
    - ploidy
    - variant type
    - local phase
    """
    vc_info = VarCallInfo(args.vc_info)
    ref_db = RefDB(args.refnames, args.fastas)
    hm = hapmatch.HapMatch(vc_info, ref_db)
    hapmatch.USE_PHASE_PROBS = args.use_phase_probs
    logging.info('HapMatch prepared')
    name_subtypes = {}
    if args.subtype_map:
        tab = pd.read_table(args.subtype_map)
        name_subtypes = dict(zip(tab['name'], tab['subtype']))
    print ('refname', 'ref_cn', 'gid', 'grank', 'gmatch', 'gedit', 'gscore', 'nblocks', 'hap_id', 'name', 'subtype', 'match', 'edit', 'edits', 'block_ids', sep='\t')
    for refname in hm.refnames:
        recs = hm.get_bests(refname, threshold=args.threshold, score=args.score)
        if not recs:
            logging.warning('No matched contigs for %s', refname)
        for gen_id, rec in enumerate(recs, 1):
            for hap_id, hrec in enumerate(rec['haps'], 1):
                name = hrec['name']
                print (
                    refname,
                    rec['ref_cn'],
                    gen_id,
                    rec['rank'],
                    rec['match'],
                    rec['edit'],
                    rec['score'],
                    rec['nblocks'],
                    hap_id,
                    name,
                    name_subtypes.get(name, name),
                    hrec['match'],
                    hrec['edit'],
                    ','.join('{0}:{1}>{2}'.format(*ed) for ed in sorted(hrec['edit_details'])),
                    ','.join(map(str, rec['block_ids'])),
                sep='\t')


@command.add_sub
@argument('input_dir')
def mc_vg_read(args):
    fnames = sorted(glob.glob(fmt('{args.input_dir}/*.vg')))
    for fname in fnames:
        logging.info('Loading %s', fname)
        #sh.call('sleep 1')
        with open(fname) as fi:
            vg = dill.load(fi)
            subs = vg.get_subgraphs()
            graph = vg.get_graph()
            nnodes = len(graph.nodes())
            logging.info('Number of nodes for %s : %s', vg.refname, nnodes)
            logging.info('Number of subgraphs for %s : %s', vg.refname, len(subs))


@command.add_sub
@argument('vc_info')
@argument('--bam', required=True)  # should change option in future ?
@argument('-r', '--regions', nargs='*')
@argument('-m', '--min-second-bases', type=int, default=2)
def mc_read_match(args):
    """
    Result info include
    - ploidy
    - variant type
    - local phase
    """
    vc_info = VarCallInfo(args.vc_info)
    with Samfile(args.bam) as sam:
        smb = SamModelBuilder(sam, regions=args.regions, min_second_bases=args.min_second_bases)

    print ('refname', 'block_id', 'pos', 'hap_bases', 'error', 'match', 'no_error',
            'ncands', 'errors', 'matches', 'no_errors', 'pair_errors', sep='\t')

    for ref in smb.model.refs:
        rname = ref.name
        pos_vars = {}
        for v in ref.variants:
            pos_vars[v.pos] = v # read selection is ignored
        #print (ref)

        vc1 = vc_info.get_var_call(rname)
        if vc1 is None:
            logging.info('No variant calls for %s', rname)
            continue
        blocks = list(vc1.iter_blocks())
        for rec in blocks: #vc1.iter_blocks():
            block_id = rec['block_id']

            phased_pos_bases = []   # [{pos: base}]
            for bases in rec['phased']:
                phased_pos_bases.append(dict(zip(rec['pos'], bases)))

            for i, pos in enumerate(rec['pos']):
                phased_bases = [hap[i] for hap in rec['phased']]
                if pos not in pos_vars:
                    continue
                v = pos_vars[pos]
                errors = [0] * vc1.copy_number
                matches = [0] * vc1.copy_number
                no_errors = [0] * vc1.copy_number
                ncands = len(v.candidates)

                pair_sites = defaultdict(lambda : [0] * vc1.copy_number)
                pair_errors = defaultdict(lambda : [0] * vc1.copy_number)

                for cand in v.candidates:
                    cand_bases = [(pos1, cand.bases[pos1]) for pos1 in cand.pos]
                    es = []
                    ms = []
                    no_es = []
                    # pair_bases = {}   # pos1:0|4:0|2  # pos1 : phased | error
                    # pair_matches = {}
                    pair_bases = {}
                    pair_matches = {}

                    # loop for each haplotype
                    for hap_id, pos_bases in enumerate(phased_pos_bases):
                        e = len(list(filter(lambda pos1, base1: pos_bases.get(pos1) and pos_bases[pos1] != base1, cand_bases)))
                        m = len(list(filter(lambda pos1, base1: pos_bases.get(pos1) and pos_bases[pos1] == base1, cand_bases)))

                        if cand.bases[pos] == pos_bases[pos]:
                            for pos1, base1 in cand_bases:
                                if pos_bases.get(pos1):
                                    pair_sites[pos1][hap_id] += 1
                                    pair_errors[pos1][hap_id] += (base1 != pos_bases[pos1])

                        es.append(e)
                        ms.append(m)
                        no_es.append(e == 0 and m >= 2)  # no errors and phase informative

                    #hap_id = np.argmin(es)  # select minimum error haplotype
                    hap_id = np.argmax(ms)  # select maximum match haplotype
                    errors[hap_id] += es[hap_id]
                    matches[hap_id] += ms[hap_id]
                    no_errors[hap_id] += no_es[hap_id]

                error = sum(errors)
                match = sum(matches)
                no_error = sum(no_errors)
                print (rname, block_id, pos,
                    '|'.join(phased_bases),
                    error, match, no_error, ncands,
                    '|'.join(map(str, errors)),
                    '|'.join(map(str, matches)),
                    '|'.join(map(str, no_errors)),
                    ' '.join('{0}:{1}:{2}'.format(
                        pos,
                        '|'.join(map(str, pair_sites[pos])),
                        '|'.join(map(str, pair_errors[pos])),
                    ) for pos in sorted(pair_sites.keys()) if sum(pair_errors[pos]) > 0),  # only pairs with errors
                    sep='\t',
                )


class VarCallTrioCheck(object):
    def __init__(self, parents, children):
        """
        parents: [VarCallInfo]
        children: [VarCallInfo]
        """
        assert len(parents) == 2
        assert len(children) >= 1
        self._parents = parents
        self._children = children
        self.variants = []
        self._ref_pos_set = defaultdict(set)   # {ref: set(pos)}

        logging.info('Merging multiple variants')
        for vc_info in chain(parents, children):
            for rname in vc_info.refnames:
                vc = vc_info.get_var_call(rname)
                self._ref_pos_set[rname].update(vc.iter_poss())

        for rname in self._ref_pos_set.keys():
            p_var_calls = [p.get_var_call(rname) for p in parents]
            c_var_calls = [c.get_var_call(rname) for c in children]
            for pos in sorted(self._ref_pos_set[rname]):
                p_pos_infos = [call and call.get_pos_info(pos) for call in p_var_calls]
                c_pos_infos = [call and call.get_pos_info(pos) for call in c_var_calls]

                rec = {
                    'refname': rname,
                    'pos': pos,
                    'parent_blocks': [info and info['block_id'] for info in p_pos_infos],
                    'parent_bases': [info and info['bases'] for info in p_pos_infos],
                    'parent_has_calls': [info and info['has_call'] for info in p_pos_infos],
                    'child_blocks': [info and info['block_id'] for info in c_pos_infos],
                    'child_bases': [info and info['bases'] for info in c_pos_infos],
                    'child_has_calls': [info and info['has_call'] for info in c_pos_infos],
                    'child_matches': [self._check_pos_info_match(info, p_pos_infos[0], p_pos_infos[1]) for info in c_pos_infos],
                    'child_contained': [self._check_pos_info_contained(info, p_pos_infos[0], p_pos_infos[1]) for info in c_pos_infos],
                }
                self.variants.append(rec)

    def _check_pos_info_match(self, child_pos_info, parent_pos_info1, parent_pos_info2):
        """
        child_pos_info: PosInfo
        parent_pos_info1: PosInfo
        parent_pos_info2: PosInfo

        PosInfo = {'block_id', 'bases', 'has_call'}
        """
        if not (child_pos_info and parent_pos_info1 and parent_pos_info2):
            return
        cbases = child_pos_info['bases'].split('|')
        pbases1 = parent_pos_info1['bases'].split('|')
        pbases2 = parent_pos_info2['bases'].split('|')

        if cbases == ['']:  # copy number: 0
            return True
        if len(cbases) == 1:   # copy number: 1
            return cbases[0] in pbases1 or cbases[0] in pbases2
        if len(cbases) == 2:   # copy number: 2
            return (
                    (cbases[0] in pbases1 and cbases[1] in pbases2)
                 or (cbases[0] in pbases2 and cbases[1] in pbases1)
            )

    def _check_pos_info_contained(self, child_pos_info, parent_pos_info1, parent_pos_info2):
        """
        child_pos_info: PosInfo
        parent_pos_info1: PosInfo
        parent_pos_info2: PosInfo

        PosInfo = {'block_id', 'bases', 'has_call'}
        """
        if not (child_pos_info and parent_pos_info1 and parent_pos_info2):
            return
        cbases = Counter(child_pos_info['bases'].split('|'))
        pbases1 = Counter(parent_pos_info1['bases'].split('|'))
        pbases2 = Counter(parent_pos_info2['bases'].split('|'))
        cbases.pop('', 0)  # remove empty case
        return all([c <= pbases1[b] + pbases2[b] for b, c in cbases.items()])


@command.add_sub
@argument('-p', '--parents', required=True, nargs=2, help='vc_info')
@argument('-c', '--children', required=True, nargs='+', help='vc_info')
@argument('-N', '--parent-names', nargs='+')
@argument('-C', '--children-names', nargs='+')
def mc_trio_check(args):
    """
    Emit VCF ?
    """
    pos_infos = PosInfo.parse(args.pos_info)
    vg = VariantGraph.parse(args.graph)
    infer = Inferer(vg, pos_infos)
    while not infer.is_converged:
        infer.update()
        infer.save_info()
    infer.emit_vcf()
    pname = ','.join(args.parent_names) if args.parent_names else '.'
    cname = ','.join(args.children_names) if args.children_names else '.'
    pinfos = [VarCallInfo(parent) for parent in args.parents]
    cinfos = [VarCallInfo(child) for child in args.children]

    check = VarCallTrioCheck(parents=pinfos, children=cinfos)
    header = 'refname pos parents children pblocks pcalls pbases cblocks ccalls cbases matches contained'.split(' ')

    print (*header, sep='\t')
    matches = [Counter() for _ in args.children]
    het_matches = [Counter() for _ in args.children]
    for v in check.variants:
        print (
            v['refname'],
            v['pos'],
            pname,
            cname,
            ','.join(str(-1 if x is None else x) for x in v['parent_blocks']),
            ','.join(str(int(False if x is None else x)) for x in v['parent_has_calls']),
            ','.join(map(str, v['parent_bases'])),
            ','.join(str(-1 if x is None else x) for x in v['child_blocks']),
            ','.join(str(int(False if x is None else x)) for x in v['child_has_calls']),
            ','.join(map(str, v['child_bases'])),
            ','.join(str(-1 if x is None else int(x)) for x in v['child_matches']),
            ','.join(str(-1 if x is None else int(x)) for x in v['child_contained']),
            #sum(v.oks),
            #','.join(map(str, v['oks'])),
            #','.join(map(str, v['ibs'])),
            sep='\t')

        # stat
        for i, m in enumerate(v['child_matches']):
            matches[i][m] += 1
        for i, (m, cb) in enumerate(zip(v['child_matches'], v['child_bases'])):
            if cb and len(set(cb.split('|'))) >= 2:  # over 2 alleles
                het_matches[i][m] += 1

    for i in range(len(args.children)):
        mc = matches[i]
        hmc = het_matches[i]
        logging.info('Performance for child %s ::: match: %s/%s (%.3f) uncalls: %s het_match: %s/%s (%.3f)',
                i + 1,
                mc[1],
                mc[0] + mc[1],  # calls
                1. * mc[1] / (mc[0] + mc[1]) if (mc[0] + mc[1]) else 0,
                mc[None],
                hmc[1],
                hmc[0] + hmc[1],  # het calls
                1. * hmc[1] / (hmc[0] + hmc[1]) if (hmc[0] + hmc[1]) else 0,
        )


class HapCallTrioChecker(object):
    """
        - sample
        - ref_cn
        - gid
        - grank
        - subtype
        - edit
        - edits
    """
    def __init__(self, hap_calls):
        self._tab = pd.read_table(hap_calls)

    def iter_trios(self, sample, father=None, mother=None, children=None):
        """
        """
        tab1 = self._tab[self._tab['grank'] == 1].fillna('')
        #tab1.sample sample
        fsubs = []  #  [(gid, subtype, edits)]
        msubs = []  #  [(gid, subtype, edits)]

        if father:
            tab2 = tab1[tab1['sample'] == father]
            fsubs = [(row.gid, row.subtype, row.edits) for row in tab2.itertuples()]

        if mother:
            tab2 = tab1[tab1['sample'] == mother]
            msubs = [(row.gid, row.subtype, row.edits) for row in tab2.itertuples()]

        # TODO children

        #logging.info((tab1['sample'].unique(), sample))
        for rec in tab1[tab1['sample'] == sample].sort_values('gid').itertuples():
            fcounts = Counter(gid for gid, sub, edits in fsubs if sub == rec.subtype and edits == rec.edits)
            mcounts = Counter(gid for gid, sub, edits in msubs if sub == rec.subtype and edits == rec.edits)
            fgids = [gid for gid, cn in sorted(fcounts.items())]
            fcns = [cn for gid, cn in sorted(fcounts.items())]
            mgids = [gid for gid, cn in sorted(mcounts.items())]
            mcns = [cn for gid, cn in sorted(mcounts.items())]

            yield {
                    'sample': rec.sample,
                    'father': father or '.',
                    'mother': mother or '.',
                    'ref_cn': rec.ref_cn,
                    'subtype': rec.subtype,
                    'gid': rec.gid,
                    'grank': rec.grank,
                    'match': rec.match,
                    'edit': rec.edit,
                    'edits': rec.edits,
                    'father_gids': ','.join(str(t) for t in fgids) or '.',
                    'father_cns': ','.join(str(t) for t in fcns) or '.',
                    'mother_gids': ','.join(str(t) for t in mgids) or '.',
                    'mother_cns': ','.join(str(t) for t in mcns) or '.',
            }


@command.add_sub
@argument('hap_calls')
@argument('-f', '--fam-file', required=True, help='Family info')
def mc_trio_hap_check(args):
    """
    Fam file:
        1. Family ID
        2. Sample ID (cannot be 0)
        3. Father ID (0 if NA)
        4. Mother ID (0 if NA)
        5. Sex code (1=male, 2=female, 0=unknown)
        (6+ not used)

    Required fileds in hap_calls:
        # print ('refname', 'ref_cn', 'gid', 'grank', 'gmatch', 'gedit', 'gscore', 'nblocks', 'hap_id', 'name', 'subtype', 'match', 'edit', 'edits', 'block_ids', sep='\t')
        - sample
        - ref_cn
        - grank
        - subtype
        - edit
        - edits
    """
    checker = HapCallTrioChecker(args.hap_calls)
    names = ['family', 'sample', 'father', 'mother', 'sex']
    fam_tab = pd.read_table(args.fam_file, header=0, names=names, sep='\t')   # TODO 0 => None?

    header = 'sample father mother gid ref_cn subtype edit match edits father_gids father_cns mother_gids mother_cns'.split(' ')  # list only matched
    print (*header, sep='\t')
    for rec in fam_tab.itertuples():
        sample = rec.sample
        father = rec.father
        mother = rec.mother
        sample = str(sample)
        assert str(sample) != '0', 'Invalid fam file ("0" is not allowed in FAM file)'
        father = None if str(father) == '0' else father
        mother = None if str(mother) == '0' else mother

        if father is None and mother is None:
            logging.info('No data of parents were found for %s', sample)
            continue

        for rec in checker.iter_trios(sample=sample, father=father, mother=mother):
            print (*(rec[name] for name in header), sep='\t')



@command.add_sub
@argument('bam')
@argument('-r', '--regions', nargs='*')
@argument('--hap-depth', type=float, default=20.)
@argument('-g', '--gref', required=True)
@argument('--max-steps', type=int, help='max steps after every node has been updated')
@argument('-m', '--min-second-bases', type=int, default=2)
def mc_path_call2(args):
    with open(args.gref) as fp:
        fasta = Fasta(fp)
        contigs = {contig.name: contig.seq.upper() for contig in fasta.contigs}

    with Samfile(args.bam) as sam:
        smb = SamModelBuilder2(sam, regions=args.regions, min_second_bases=args.min_second_bases, contigs=contigs,
                ploidies={}, # TODO
                hap_depth=args.hap_depth)
        summary = smb.summary
        show_model(smb.model, verbosity=args.verbose)
    from infer4 import InferModel
    im = InferModel(smb.model, hap_depth=args.hap_depth)
    #im.init_best_het()
    im.init_all_variants()

    def show_status(step, status):
        v = status['variant']
        gen_probs = im.get_variant_probs(v)   # [{'paths': [[bases]], 'bases': [joined bases], 'prob'; prob}]
        gen_idx = np.argmax([rec['prob'] for rec in gen_probs])
        best_bases = [''.join(path) for path in gen_probs[gen_idx]['paths']]
        best_prob = gen_probs[gen_idx]['prob']

        logging.info('Step %s updated %s:%s (upd:%s, rem:%s) to %s (%.3f)',
                step, v.ref.name, v.pos, status['queued_count'], len(status['remained']), '|'.join(best_bases), best_prob)
        if args.verbose >= 2:
            show_variant_info2(im, summary, contigs, stdout=sys.stderr)   # show variant, alts, depth, probs.

    def run_loop(step=0, read_propagate=False):
        while 1:
            status = im.update(read_propagate=read_propagate)  # more info
            if status['done']:
                break
            step += 1
            show_status(step, status)
            if len(status['remained']) == 0 and args.max_steps and step >= args.max_steps:
                break
        return step

    #step = run_loop()
    if hasattr(im, 'enqueue_all_het_variants'):
        logging.info('Starting second round inference with read propagation')
        #im.enqueue_all_het_variants()
        #im.enqueue_all_variants()
        #step = run_loop(step=step, read_propagate=True)

        #logging.info('Starting third round inference with read propagation')
        #im.enqueue_all_variants()
        #run_loop(step=step, read_propagate=True)

    show_variant_info2(im, summary, contigs, show_all_variant=True)   # show variant, alts, depth, probs.


def show_variant_info2(infer_model, sam_summary, contigs, show_all_variant=False, stdout=sys.stdout):
    print (
        'refname',
        'pos',
        'end',
        'block_id',
        'last_step',
        'min_link_pos',
        'max_link_pos',
        'link_pos',
        'bases',   # candidates
        'counts',  # 
        'ref',
        'zygosity',     # number of different bases
        'hap_bases',     # called
        'depths',
        'unmapped',
        'probs',
        sep='\t', file=stdout)

    ref_pos_variants = {(v.ref.name, v.pos): v for v in infer_model.variants}
    for rname in sam_summary.refnames:
        seq = contigs[rname]
        ploidy = infer_model.get_ploidy(rname)
        last_end = 0
        for pos, ref_base in enumerate(seq):
            if pos < last_end:
                continue
            if (rname, pos) in ref_pos_variants:
                v = ref_pos_variants[rname, pos]
                last_end = end = v.end
                ref_base = seq[pos:end]
                sinfo = infer_model.get_site_info(v)
                gen_probs = infer_model.get_variant_probs(v)   # [((b1, b2, ..), prob)]
                gen_idx = np.argmax([rec['prob'] for rec in gen_probs])
                #logging.info(gen_probs)
                #logging.info(gen_idx)
                #best_bases = [''.join(path) for path in gen_probs[gen_idx]['paths']]
                best_bases = gen_probs[gen_idx]['bases']
                best_prob = gen_probs[gen_idx]['prob']
                hap_prob = '{0:.3f}'.format(best_prob)

                phase_block = infer_model.get_phase_block(v)
                linked_pos = sinfo['linked_pos']
                last_step = sinfo['update_id']
                cand_bases = sinfo['bases']
                base_counts =  sinfo['base_counts']
                mean_depths = sinfo['mean_depths']
                um_depths = sinfo.get('depths_unmapped')
            elif ref_base == 'N' and show_all_variant:
                end = pos + 1
                phase_block = None
                linked_pos = None
                last_step = 0
                info = sam_summary.get_info(rname, pos)
                cand_bases = info['cand_bases']  if info else [] # TODO
                base_counts = info['base_counts'] if info else Counter()
                # need ploidy
                # assume homozygous of most common bases
                best_bases = [base_counts.most_common(1)[0][0] if base_counts else 'N'] * ploidy
                mean_depths = None
                um_depths = None
                hap_prob = None
            else:
                continue
            zygosity = len(set(best_bases))

            print (rname, pos, end, #sinfo,
                   phase_block,
                   last_step,
                   linked_pos and min(linked_pos),
                   linked_pos and max(linked_pos),
                   linked_pos and len(linked_pos),
                   ','.join(cand_bases),
                   ','.join(str(base_counts[b]) for b in cand_bases),
                   ref_base,
                   zygosity,
                   '|'.join(best_bases),
                   mean_depths and '|'.join('{0:.1f}'.format(d) for d in mean_depths),
                   um_depths and '|'.join('{0:.1f}'.format(d) for d in um_depths),
                   hap_prob,
                   sep='\t',
                   file=stdout)


# Using prior, without lambda, 
@command.add_sub
@argument('bam')
@argument('-r', '--regions', nargs='*')
@argument('--hap-depth', type=float, default=20.)
@argument('-g', '--gref', required=True)
@argument('--max-steps', type=int, help='max steps after every node has been updated')
@argument('--show-allocs', action='store_true')
@argument('--show-scores', action='store_true')
@argument('-m', '--min-second-bases', type=int, default=2)
def mc_path_call5(args):
    with open(args.gref) as fp:
        fasta = Fasta(fp)
        contigs = {contig.name: contig.seq.upper() for contig in fasta.contigs}

    with Samfile(args.bam) as sam:
        smb = SamModelBuilder2(sam, regions=args.regions, min_second_bases=args.min_second_bases, contigs=contigs,
                ploidies=None, # TODO
                hap_depth=args.hap_depth)
        hap_depths = {ref.name: args.hap_depth for ref in smb.model.refs}
        ploidies = {ref.name: 2 for ref in smb.model.refs}
        summary = smb.summary
        show_model(smb.model, verbosity=args.verbose)
    from infer5 import InferModel
    im = InferModel(smb.model, hap_depths=hap_depths, ploidies=ploidies)
    #im.init_best_het()
    im.init()

    def show_status(step, status):  # TODO
        v = status['variant']
        gen_probs = im.get_variant_probs(v)   # [{'paths': [[bases]], 'bases': [joined bases], 'prob'; prob}]
        gen_idx = np.argmax([rec['prob'] for rec in gen_probs])
        best_bases = [''.join(path) for path in gen_probs[gen_idx]['paths']]
        best_prob = gen_probs[gen_idx]['prob']

        #if len(set(best_bases)) > 1:  # HET
        if True:
            si = im.get_site_info(v)
            logging.info('Step %s updated %s:%s (upd:%s, rem:%s) to %s (%.3f), %s (%.3f), dp: %s, um: %s from %s',
                    step, v.ref.name, v.pos, status['queued_count'], len(status['remained']), '|'.join(best_bases), best_prob, 
                    '/'.join(''.join(bs) for bs in si['summary']['best_gt_bases']),
                    si['summary']['best_gt_prob'],
                    '|'.join('{0:.1f}'.format(x) for x in si['summary']['depths_mapped']),
                    '|'.join('{0:.1f}'.format(x) for x in si['summary']['depths_unmapped']),
                    si['summary']['base_counts'])
            if args.show_scores:
                scores = defaultdict(list)
                for c in v.candidates:
                    bs = c.bases[v.pos]
                    scores[bs].append(c.score)
                for bs, scs in scores.items():
                    logging.info('%s num: %s, min: %s, max: %s, avg: %s', bs, len(scs), min(scs), max(scs), np.mean(scs))
            if args.show_allocs:
                for sup in si['summary']['supports']:
                    c = sup['cand']
                    logging.info('u: %.2f, z: %s, pos: %s, base: %s, bs: %s',
                            sup['QU'],
                            '|'.join('{0:.2f}'.format(q) for q in sup['QZ']),
                            ''.join('.' if b is None else b for b in c.bases[v.pos]),
                            ','.join(str(x) for x in c.pos),
                            ','.join(''.join('.' if b is None else b for b in c.bases[pos]) for pos in c.pos))

        if args.verbose >= 2:
            show_variant_info2(im, summary, contigs, stdout=sys.stderr)   # show variant, alts, depth, probs.

    def run_loop(step=0):
        while 1:
            status = im.update()  # more info
            if status['done']:
                break
            elif status['variant'] is None:
                im.init()
                continue
            step += 1
            show_status(step, status)
            if len(status['remained']) == 0 and args.max_steps and step >= args.max_steps:
                break
        return step

    vcs = im.run_phasing()
    show_variant_info5(im, vcs, summary, contigs, show_all_variant=True)
    #step = run_loop()
    #show_variant_info2(im, summary, contigs, show_all_variant=True)   # show variant, alts, depth, probs.

def show_variant_info5(infer_model, vcs, sam_summary, contigs, show_all_variant=False, stdout=sys.stdout, prob_threshold=.9):
    print (
        'refname',
        'pos',
        'end',
        'block_id',
        'last_step',
        'min_link_pos',
        'max_link_pos',
        'link_pos',
        'bases',   # candidates
        'counts',  # 
        'ref',
        'zygosity',     # number of different bases
        'hap_bases',     # called
        'depths',
        'unmapped',
        'probs',
        sep='\t', file=stdout)

    ref_pos_variants = {(v.ref.name, v.pos): v for v in infer_model.variants}
    vc_variants = {(v.ref.name, v.pos): v for vc in vcs for v in vc.get_variants()}
    var_phase_blocks = {v: vc.get_variants()[0].pos for vc in vcs for v in vc.get_variants()}
    for rname in sam_summary.refnames:
        seq = contigs[rname]
        ploidy = infer_model.get_ploidy(rname)
        last_end = 0
        for pos, ref_base in enumerate(seq):
            if pos < last_end:
                continue
            if (rname, pos) in ref_pos_variants:
                v = ref_pos_variants[rname, pos]
                last_end = end = v.end
                ref_base = seq[pos:end]
                sinfo = infer_model.get_site_info(v)
                gen_probs = infer_model.get_variant_probs(v)   # [((b1, b2, ..), prob)]
                gen_idx = np.argmax([rec['prob'] for rec in gen_probs])  # TODO fix it
                best_bases = gen_probs[gen_idx]['bases']
                best_prob = gen_probs[gen_idx]['prob']
                #logging.info((v.pos, best_prob))
                if best_prob is None or best_prob < prob_threshold:
                    best_bases = ['.'.join('N' for b in bs.split('.')) for bs in best_bases]
                hap_prob = '{0:.3f}'.format(best_prob)

                phase_block = var_phase_blocks.get(v)
                #phase_block = infer_model.get_phase_block(v)
                linked_pos = sinfo['linked_pos']
                last_step = sinfo['update_id']
                cand_bases = sinfo['bases']    # TODO
                base_counts =  sinfo['base_counts']  # TODO
                mean_depths = sinfo['mean_depths']  # TODO
                um_depths = sinfo.get('depths_unmapped')  # TODO
            elif ref_base == 'N' and show_all_variant:
                end = pos + 1
                phase_block = None
                linked_pos = None
                last_step = 0
                info = sam_summary.get_info(rname, pos)
                cand_bases = info['cand_bases']  if info else [] # TODO
                base_counts = info['base_counts'] if info else Counter()
                # need ploidy
                # assume homozygous of most common bases
                best_bases = [base_counts.most_common(1)[0][0] if base_counts else 'N'] * ploidy
                mean_depths = None
                um_depths = None
                hap_prob = None
            else:
                continue
            zygosity = len(set(best_bases))

            print (rname, pos, end, #sinfo,
                   phase_block,
                   last_step,
                   linked_pos and min(linked_pos),
                   linked_pos and max(linked_pos),
                   linked_pos and len(linked_pos),
                   ','.join(cand_bases),
                   ','.join(str(base_counts[b]) for b in cand_bases),
                   ref_base,
                   zygosity,
                   '|'.join(best_bases),
                   mean_depths and '|'.join('{0:.1f}'.format(d) for d in mean_depths),
                   um_depths and '|'.join('{0:.1f}'.format(d) for d in um_depths),
                   hap_prob,
                   sep='\t',
                   file=stdout)


@command.add_sub
@argument('bam')
@argument('-r', '--regions', nargs='*')
@argument.exclusive(
    (argument('-c', '--copy-number', type=int, default=2)
    .add_argument('--hap-depth', type=float, default=20.)),
    argument('-t', '--table', help='contig, hap_depth, and copy_number are required')
)
@argument('-g', '--gref', required=True)
@argument('--no-phase', action='store_true')
@argument('-m', '--min-second-bases', type=int, default=2)
def mc_path_call6(args):
    with open(args.gref) as fp:
        fasta = Fasta(fp)
        contigs = {contig.name: contig.seq.upper() for contig in fasta.contigs}

    with Samfile(args.bam) as sam:
        smb = SamModelBuilder2(sam, regions=args.regions, min_second_bases=args.min_second_bases, contigs=contigs)

    if not args.table:
        hap_depths = {ref.name: args.hap_depth for ref in smb.model.refs}
        ploidies = {ref.name: args.copy_number for ref in smb.model.refs}
    else:
        tab = pd.read_table(args.table)
        hap_depths = dict(zip(tab.contig, tab.hap_depth))
        ploidies = dict(zip(tab.contig, tab.copy_number))
    show_model(smb.model, verbosity=args.verbose)
    from .infer6 import InferModel
    im = InferModel(smb, hap_depths=hap_depths, ploidies=ploidies)
    #im.init_best_het()
    #im.init()
    #im.run_through_variants()
    if args.no_phase:
        im.run_genotyping()
    else:
        im.run_haplotyping()
        start_var = None
        #start_var = smb.model.refs[0].get_variant(1203)
        #im.run_phase_variants(start_var=start_var)
    im.show_variant_info(show_all_variant=True)
