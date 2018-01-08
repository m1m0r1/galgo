from __future__ import print_function
from argtools import command, argument
from builtins import zip
import logging
import pandas as pd
import os
import gzip
import tempfile
import glob
from collections import namedtuple
from ..utils import file_split, iter_tabs, collect_while
from ..utils import format as fmt
from .. import sh


@command.add_sub
@argument('input')
@argument('--sep', default='\t')
def transpose(args):
    with open(args.input) as fp:
        for c in zip(*(l.strip().split(args.sep) for l in fp.readlines() if l.strip())):
             print (*c, sep=args.sep)


@command.add_sub   # similar to 'csvcut -t -n table'
@argument('input')
def header(args):
    with open(args.input) as fp:
        header = next(iter_tabs(fp), None)
        if header:
            for name in header:
                print (name)


@command.add_sub   # similar to 'csvcut -t -n table'
@argument('value')
@argument('-i', '--input', default='/dev/stdin')
@argument('-n', '--name', help='column name to add')
@argument('-c', '--column', type=int, help='Where to add column. default is last')
# argument('--no-template')
# argument('--eval')
# argument('--no-header')
def add_column(args):
    """
    add_column 'column_{col_name or col_number}' -n name
    add_column 'column_{__index__:05d}' -n id
    """
    use_header = bool(args.name)

    if args.column is None:  # last
        add_col = lambda row, value: row.append(value)
    else:
        add_col = lambda row, value: row.insert(args.column, value)

    with open(args.input) as fp:
        it = iter_tabs(fp)
        if use_header:
            header = next(it)
            add_col(header, args.name)
            print (*header, sep='\t')

        kwds = {}
        for __index__, row in enumerate(it, 1):
            if use_header:
                kwds = dict(zip(header, row))
            value = args.value.format(__index__=__index__, **kwds)
            add_col(row, value)
            print (*row, sep='\t')


@command.add_sub
@argument('file1')
@argument('file2')
def set_isec(args):
    """
    * tsv file
    * key is leftmost column

    # TODO impl. for many files
    """

    key_idx = 0
    with open(args.file1) as fp:
        keys1 = [row[key_idx] for row in iter_tabs(fp)]
    with open(args.file2) as fp:
        keys2 = [row[key_idx] for row in iter_tabs(fp)]

    key_set1 = set(keys1)
    key_set2 = set(keys2)
    key_isec = key_set1 & key_set2

    for key in keys1:
        if key in key_isec:
            print (key)


@command.add_sub
@argument('input', nargs='?', default='/dev/stdin')
@argument('-c', '--columns')
@argument('-C', '--not-columns')
def cut(args):
    with open(args.input) as fp:
        it = iter_tabs(fp)
        header = next(it, None)
        if header is None:
            return
        col_idxs = dict((col, i) for i, col in enumerate(header))
        if args.columns:
            columns = args.columns.split(',')
            idxs = [col_idxs[col] for col in columns]
        else:
            idxs = list(range(len(header)))

        if args.not_columns:
            not_idxs = set([col_idxs[col] for col in args.not_columns.split(',')])
            idxs = [i for i in idxs if i not in not_idxs]

        print (*(header[i] for i in idxs), sep='\t')
        for row in it:
            print (*(row[i] for i in idxs), sep='\t')


@command.add_sub
@argument('inputs', nargs='*', default='/dev/stdin')
@argument('--header', action='store_true')
@argument('--lenient', action='store_true', help='if file does not exist, just warn')
def cat(args):
    """ concat inputs with header
    """
    def iter_files():
        for fn in args.inputs:
            if not os.path.exists(fn):
                if args.lenient:
                    logging.warning('%s does not exist and skipped', fn)
                    continue
                else:
                    raise Exception('File {0} does not found'.format(fn))
            yield fn

    for i, fn in enumerate(iter_files()):
        if i == 0:
            sh.call(['cat', fn])
        else:
            if args.header:
                sh.call(['sed', '1d', fn])
            else:
                sh.call(['cat', fn])


@command.add_sub
@argument('input', nargs='?', default='/dev/stdin')
@argument('-m', '--max-char', default=20, type=int)
@argument('--sep', default='\t')
def truncate(args):
    """ concat inputs with header
    """
    sep = args.sep
    max_char = args.max_char
    with open(args.input) as fp:
        it = iter_tabs(fp, sep=sep)
        for row in it:
            row = [s[:max_char] for s in row]  # truncate
            print (*row, sep=sep)

# TODO .gz
# - running gzip after split may be sufficient
# - using piped gzip inside of file_split is better
@command.add_sub
@argument('input')
@argument('-n', '--nlines', default=1000, type=int)
@argument('-p', '--prefix', default='')
@argument('-s', '--suffix', default='')
@argument('-d', '--dir', default=None)
@argument('--header', default='0', help='header char "#" or number of header lines')
@argument('-c', '--compress', choices=['gzip', 'bgzip', 'bzip2'], help='specify compressor')
def split(args):
    logging.info('The number of lines per file: %s', args.nlines)
    input = '/dev/stdin' if args.input == '-' else args.input

    with open(input) as fp:
        if args.header.isdigit():
            nheader = int(args.header)
            assert nheader >= 0
            logging.info('Header line: %s', nheader)
            header = '\n'.join([next(fp).rstrip() for _ in range(nheader)])
        elif args.header == '#':
            head, fp = collect_while(lambda x: x.startswith(args.header), fp)
            header = '\n'.join([line.rstrip() for line in head])
        else:
            logging.warning('Invalid header choice: %s', args.header)
            raise NotImplementedError

        for name in file_split(fp, prefix=args.prefix, suffix=args.suffix, dir=args.dir, nlines=args.nlines, header=header):
            logging.info('Wrote into %s', name)
            if args.compress:   # adhoc
                sh.call([args.compress, name])
                csuffix = {'gzip': '.gz', 'bgzip': '.gz', 'bzip2': '.bz2'}[args.compress]
                name += csuffix
                logging.info('Compressed the file to %s', name)
            print (name)


class _group_fns:
    @staticmethod
    def first(vals):
        return vals[0]

    @staticmethod
    def concat(vals, sep=','):
        return sep.join(vals)

    @staticmethod
    def length(vals):
        return len(vals)

    @staticmethod
    def sum(vals):
        return sum(map(float, vals))

_group_fn_list = [_attr for _attr in dir(_group_fns)
        if (not _attr.startswith('__') and not _attr.endswith('__') and callable(getattr(_group_fns, _attr)))]

@command.add_sub
@argument('-c', '--columns', nargs='+', help='column specs')
@argument('--no-header', '-H', dest='use_header', action='store_false')
@argument('--default-fn', '-f', default='first', choices=_group_fn_list, help='default group function')
@argument.exclusive(
argument('-i', '--inputs', nargs='+', help='input files'),
argument('-l', '--input-list', help='list of input file names'),
required=True
)
def summarize_row(args):
    """
    -c {col1} {col2} ...
    Availble specifiers :
        col1:col2               ---  concat(col1)
        col1::first
        col1:col1-sum:sum
        col1::concat            ---  concat(col1) as col1
        col1::length            ---  length(col1) as col1
        0                       ---  default_fn(filelist) (1st column) as filelist
        +n:sample               ---  (n+1)-th column of filelist
    """
    if not args.use_header:
        raise NotImplementedError('Not implemented yet')

    input_tab = None
    if args.inputs:
        logging.info('Reading %s files', len(args.inputs))
        fnames = list(args.inputs)
    elif args.input_list:
        logging.info('Read from %s', args.input_list)
        input_tab = pd.read_table(args.input_list, header=None)
        fnames = list(input_tab[0])
    else:
        logging.error('Require input file list')
        return 1

    rows_list = []  # [rows]
    headers = []  # [header]
    # TODO read using threads ?
    for i, fname in enumerate(fnames):
        logging.info('Loading file %s', fname)
        if fname.endswith('.gz'):
            _open = gzip.open
        else:
            _open = open
        with _open(fname) as fp:
            it = iter_tabs(fp)
            header = next(it)
            headers.append(header)
            if i == 0:
                for row in it:
                    rows_list.append([row])
            else:
                for j, row in enumerate(it):
                    rows_list[j].append(row)

    col_specs = []  # (col_name, getter, fn_name)
    if args.columns:
        columns = args.columns
    else:
        columns = headers[0]  # use header of first sample

    for col_spec in columns:
        col_specs.append(parse_col_spec(col_spec, headers=headers, fnames=fnames, input_tab=input_tab, default_fn=args.default_fn))

    print (*[sp.name for sp in col_specs], sep='\t')
    for rows in rows_list:
        out_row = []
        for sp in col_specs:
            vals = [sp.getter(row, fidx) for fidx, row in enumerate(rows)]
            out_val = sp.fn(vals)
            out_row.append(out_val)
        print (*out_row, sep='\t')


Colspec = namedtuple('ColSpec', 'name getter fn')
def parse_col_spec(col_spec, headers, fnames, input_tab, default_fn='first'):
    tokens = col_spec.split(':')
    if len(tokens) == 3:
        col, col_name, fn_name = tokens
        if not col_name:
            col_name = col
    elif len(tokens) == 2:
        col, col_name = tokens
        fn_name = None
    elif len(tokens) == 1:
        col, = tokens
        col_name = col
        fn_name = None
    else:
        raise NotImplementedError('Unreadable col_spec: {0}'.format(col_spec))

    if col == '0':
        getter = lambda row, fidx: fnames[fidx]
    elif col[0] == '+' and col[1:].isdigit():
        col_idx = int(col[1:])
        getter = lambda row, fidx: input_tab[col_idx][fidx]
    else:
        col_idxs = [header.index(col) for header in headers]
        getter = lambda row, fidx: row[col_idxs[fidx]]

    if fn_name is None:
        fn_name = default_fn
    fn = getattr(_group_fns, fn_name)

    return Colspec(name=col_name, getter=getter, fn=fn)


@command.add_sub
@argument('-c', '--columns', nargs='+', help='column specs')
@argument('--no-header', '-H', dest='use_header', action='store_false')
@argument('--default-fn', '-f', default='first', choices=_group_fn_list, help='default group operator')
@argument('-l', '--input-list', help='list of input file names', required=True)
@argument('-s', '--chunk-size', default=5000, help='default operator')
@argument('--dir', help='working directory')
@argument('--memory-reduce', type=float, default=8.)
@argument('--read-cmd', default='cat {0}', help='read cmd for input list')
#@argument('--keep-dir', action='store_true')
def summarize_row_chunked(args):
    """
    """
    if not args.use_header:
        raise NotImplementedError('Not implemented yet')

    from ..job_queue import UGEQueue
    queue = UGEQueue()

    if args.dir:
        work_dir = args.dir
        sh.call('mkdir -p {dir}'.format(dir=work_dir))
    else:
        work_dir = tempfile.mkdtemp()
    logging.info('Working directory: %s', work_dir)

    logging.info('Read from %s', args.input_list)
    input_tab = pd.read_table(args.input_list, header=None)
    fnames = list(input_tab[0])
    input_len = len(fnames)
    input_fn = '{dir}/input.fn'.format(dir=work_dir)
    input_tab[[0]].to_csv(input_fn, header=False, index=False, sep='\t')  # save only file names

    import sys
    galgo = ' '.join([sys.executable, sys.argv[0]])

    # create dirs
    sh.call('mkdir -p {dir}/inputs'.format(dir=work_dir))
    sh.call('mkdir -p {dir}/splits'.format(dir=work_dir))

    # split
    # less is used for read gz file if possible
    read_cmd = args.read_cmd.format('$fname')
    cmd = '''
    mkdir -p {dir}/inputs/$SGE_TASK_ID
    fname="$(sed -n ${{SGE_TASK_ID}}p {input_fn})"
    {galgo} split --dir {dir}/inputs/$SGE_TASK_ID --header 1 -n {nlines} --suffix .txt <({read_cmd})
    '''.format(galgo=galgo, dir=work_dir, nlines=args.chunk_size, input_fn=input_fn, read_cmd=read_cmd)
    queue.array_call(cmd, name='summarize_row_split', memory=8, ntask=input_len)

    first_splits = glob.glob('{dir}/inputs/1/*.txt'.format(dir=work_dir))
    split_len = len(first_splits)
    logging.info('Split to %s files.', split_len)

    # create input_list for each split job
    for split_idx in range(1, split_len + 1):
        input_files = [
            '{dir}/inputs/{input_idx}/{split_idx}.txt'.format(dir=work_dir, input_idx=input_idx, split_idx=split_idx)
            for input_idx in range(1, input_len + 1)
        ]
        input_list1 = '{dir}/splits/{split_idx}.list'.format(dir=work_dir, split_idx=split_idx)
        input_tab1 = input_tab.copy()
        input_tab1[0] = input_files  # overwrite with split files
        logging.info('Saving list %s', input_list1)
        input_tab1.to_csv(input_list1, header=False, index=False, sep='\t')

    # setup summarize_row options
    summary_opts = []
    if not args.use_header:
        summary_opts.append('-H')
    if args.default_fn:
        summary_opts.append('-f {0}'.format(args.default_fn))
    if args.columns:
        summary_opts.append('-c {0}'.format(' '.join(args.columns)))  #TODO quote

    cmd = '''
    {galgo} summarize_row -v -l {dir}/splits/$SGE_TASK_ID.list {opts} | gzip -c > {dir}/splits/$SGE_TASK_ID.txt.gz
    '''.format(galgo=galgo, dir=work_dir, opts=' '.join(summary_opts))
    queue.array_call(cmd, name='summarize_row_reduce', memory=args.memory_reduce, ntask=split_len)

    # concat
    logging.info('concat all')
    for split_idx in range(1, split_len + 1):
        if split_idx == 1:
            sh.call('zcat {dir}/splits/{split_idx}.txt.gz'.format(dir=work_dir, split_idx=split_idx))
        else:
            sh.call('zcat {dir}/splits/{split_idx}.txt.gz | sed 1d'.format(dir=work_dir, split_idx=split_idx))


@command.add_sub
@argument('input_xsls', help='input.xsls')
@argument('-l', '--list', action='store_true', help='list sheetnames')
@argument('-n', '--sheet-name', default='0')
@argument('-s', '--sep', default='\t')
@argument('-o', '--output', default='/dev/stdout')
def excel(args):
    xls = pd.ExcelFile(args.input_xsls)
    if args.list:
        for name in xls.sheet_names:
            print (name)
        return 0

    n = args.sheet_name
    if args.sheet_name.isdigit():
        n = int(args.sheet_name)
    parsed = xls.parse(sheetname=n)
    parsed.to_csv(args.output, header=True, index=False, sep=args.sep)


@command.add_sub
@argument('-o', '--output', default='output.xls', help='output.xsls')
@argument('-i', '--inputs', required=True, nargs='+', default='tsv files')
@argument('-n', '--names', nargs='+')
def excelify(args):
    logging.info('Writing to %s', args.output)
    writer = pd.ExcelWriter(args.output, engine='xlsxwriter')
    names = args.names
    for i, ifile in enumerate(args.inputs):
        tab = pd.read_table(ifile, dtype='str')
        if names:
            sheet_name = names[i]
        else:
            sheet_name = 'Sheet{0}'.format(i + 1)
        logging.info('Writing %s', sheet_name)
        tab.to_excel(writer, sheet_name, index=False)
    writer.save()


@command.add_sub
@argument('table')
@argument('-c', '--column')
@argument('-o', '--output', default='plot_hist.pdf')
@argument('-b', '--bins', default=50, type=int)
@argument('--xlim', nargs=2, type=float)
@argument('--ylim', nargs=2, type=float)
def plot_hist(args):
    logging.info('Save plot to %s', args.output)
    tab = pd.read_table(args.table)
    print (tab.describe())
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.figure()
    ax = plt.gca()
    if args.xlim:
        ax.set_xlim(args.xlim[0], args.xlim[1])
    if args.ylim:
        ax.set_ylim(args.ylim[0], args.ylim[1])
    tab[args.column].hist(ax=ax, bins=args.bins)
    plt.savefig(args.output)
    plt.close()


import argparse
@command.add_sub
@argument('-s', '--sort', choices=['time', 'calls', 'cumulative', 'stdname'], default='time')
@argument('-o', '--output', default=None)
@argument('args', nargs=argparse.REMAINDER)
def profile(args):
    """ Utility for parse cProfile result
    """

    def run(output):
        import os
        import sys
        output_opt = fmt('-o {output}')
        argv = ' '.join(map(sh.quote, args.args))
        cmd = fmt('{sys.executable} -m cProfile -s {args.sort} {output_opt} {sys.argv[0]} {argv}')
        sh.call(cmd)
        import pstats
        st = pstats.Stats(output, stream=sys.stderr)
        st.sort_stats(args.sort).print_stats()  # stdname calls time cumulative

    if args.output:
        run(args.output)
    else:
        with tempfile.NamedTemporaryFile() as fp:
            run(fp.name)


@command.add_sub
@argument('pstat_file')
@argument('-s', '--sort', choices=['time', 'calls', 'cumulative', 'stdname'], default='time')
def show_pstats(args):
    """ Utility for parse cProfile result
    """
    import pstats
    st = pstats.Stats(args.pstat_file)
    st.sort_stats(args.sort).print_stats()  # stdname calls time cumulative


@command.add_sub
@argument('-m', '--memory', default=4, type=float, help='Requesting memory in gigabyte')
#@argument('--smem-boost', default=1.5, type=float, help='difference of smem')
@argument('-q', '--queue', default=None)
@argument('-N', '--name', default=None)
@argument('-t', '--task', default=None)
#@argument('--hold_jid')
@argument('--hold', default=None)
@argument('-s', '--slot', default=1, type=int)
@argument('--dry-run', action='store_true')
@argument('-e', '--engine', choices=['uge', 'local'], default='local')
@argument('args', nargs='+')
def run(args):
    from ..job_queue import UGEQueue, LocalQueue
    if args.engine == 'local':
        queue = LocalQueue()
    elif args.engine == 'uge':
        queue = UGEQueue()
    else:
        raise NotImplementedError(args.engine)

    if len(args.args) == 1:
        cmd = args.args[0]
    else:
        cmd = args.args

    def parse_task(task):
        # 'n[-m][:s]'
        sp = task.split(':')
        if len(sp) >= 2:
            rg, step = sp[0], int(sp[1])
        else:
            rg, step = sp[0], None
        sp = rg.split('-')
        if len(sp) >= 2:
            start, ntask = int(sp[0]), int(sp[1])
        else:
            start = ntask = int(sp[0])
        return start, ntask, step

    if args.task:
        start, ntask, step = parse_task(args.task)
        queue.array_call(cmd, queue=args.queue, memory=args.memory, name=args.name, hold=args.hold, 
                start=start, ntask=ntask, step=step, slot=args.slot, dry_run=args.dry_run)
    else:
        queue.call(cmd, queue=args.queue, memory=args.memory, name=args.name, hold=args.hold, slot=args.slot, dry_run=args.dry_run)

