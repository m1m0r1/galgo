from __future__ import print_function
import os
from six import string_types
import logging
import multiprocess
import signal
import datetime
import time
from . import sh

class BaseQueue(object):
    def call(self, args, **kwds):
        raise NotImplementedError

    def array_call(self, args, **kwds):
        raise NotImplementedError


# deprecated
class ShellQueue(BaseQueue):
    def __init__(self, max_ntask=None, dry_run=False, **kwds):
        self.job_dir = None
        self.max_ntask = max_ntask
        self.dry_run = dry_run

    def _strip_ignores(self, kwds):
        ignores = ['name', 'hold', 'memory', 'queue', 'shell', 'script_name']
        for name in ignores:
            if kwds.pop(name, None):
                logging.info('Ignore option: {0}'.format(name))
        return kwds

    def call(self, args, slot=None, dry_run=None, **kwds):
        kwds = self._strip_ignores(kwds)
        if dry_run is None:
            dry_run = self.dry_run
        return sh.call(args, dry_run=dry_run, **kwds)

    def array_call(self, args, ntask, slot=None, max_ntask=None, dry_run=None, **kwds):
        kwds = self._strip_ignores(kwds)
        if dry_run is None:
            dry_run = self.dry_run
        assert isinstance(args, string_types)
        def init_pool():
            signal.signal(signal.SIGINT, signal.SIG_IGN)

        def task(n):
            task_id = n + 1
            env = dict(os.environ, TASK_ID=task_id)
            return self.call(args, dry_run=dry_run, env=env)

        try:
            pool = Pool(self.max_ntask, init_pool)
            for result in pool.imap_unordered(task, range(ntask)):
                pass
            pool.close()
        except Exception as e:
            logging.error(e)
            pool.terminate()
            raise e
        finally:
            pool.join()

LocalQueue = ShellQueue

# TODO load queue.ini

class UGEQueue(BaseQueue):
    """
    .job/{}
    """
    def __init__(self, max_ntask=None, job_dir='.job', dry_run=False):
        self.job_dir = job_dir
        self.dry_run = dry_run
        self.queue = ''
        self.max_ntask = max_ntask
        self.shell = '/bin/bash'
        self.script_name = 'job.sh'
        self.s_vmem_boost = 1.5
        self.s_vmem_boost = 1.  # for safety
        self.slot = 1
        self.memory = 4.

    def _create_job_script(self, script_body, name=None, job_dir=None, script_name=None, shell=None, slot=None, memory=None, queue=None, hold=None, is_array=False, dry_run=None):
        if job_dir is None:
            job_dir = self.job_dir
        if name is None:
            name = 'job'
        if shell is None:
            shell = self.shell
        if script_name is None:
            script_name = self.script_name
        if dry_run is None:
            dry_run = self.dry_run
        queue = self.queue if queue is None else queue
        memory = self.memory if memory is None else memory
        slot = self.slot if slot is None else slot

        qopts = []
        # shell
        qopts.append(('-S', self.shell))

        # resource
        assert slot >= 1
        if slot > 1:
            qopts.append(('-pe', 'def_slot', str(slot)))
            memory = 1. * memory / slot
        qopts.append(('-l', '{queue},mem_req={mem_req}G,s_vmem={s_vmem}G'.format(
            queue=queue, mem_req=memory, s_vmem=memory * self.s_vmem_boost)))

        # job/log directory
        dt = datetime.datetime.now().strftime('%Y.%m.%d_%H.%M.%S.%f')
        data_dir = os.path.join(job_dir, '{name}_{dt}'.format(name=name, dt=dt))
        time.sleep(.00001)
        logging.info('Save job files to %s', data_dir)
        if not dry_run:
            os.makedirs(data_dir)

        if is_array:
            log_suffix = '.$TASK_ID.$JOB_ID.$HOSTNAME'
        else:
            log_suffix = '.$JOB_ID.$HOSTNAME'

        stdout_path = os.path.join(data_dir, 'stdout{0}'.format(log_suffix))
        stderr_path = os.path.join(data_dir, 'stderr{0}'.format(log_suffix))

        # other environment
        qopts.append(('-cwd',))
        qopts.append(('-V',))
        qopts.append(('-N', name))
        qopts.append(('-o', stdout_path))
        qopts.append(('-e', stderr_path))

        if hold:
            qopts.append(('-hold_jid', hold))

        lines = []
        lines.append('#!' + shell)
        for opt in qopts:
            opt = ' '.join(opt)
            lines.append('#$' + opt)

        lines.append(script_body)
        content = '\n'.join(lines)
        script_path = os.path.join(data_dir, script_name)
        logging.info('Script path: %s', script_path)
        logging.info('Stdout path: %s', stdout_path)
        logging.info('Stderr path: %s', stderr_path)
        logging.info('Script content:\n%s', content)
        if not dry_run:
            with open(script_path, 'w+') as writer:
                print (content, file=writer)
        return {'data_dir': data_dir, 'script_path': script_path}

    def _qsub_sync(self, script, ext_opts=None, dry_run=False):
        opts = ['-terse', '-sync', 'y']
        cmd = ['qsub'] + opts + list(ext_opts or []) + [script]
        sh.call(cmd, stdout='/dev/stderr', dry_run=dry_run)
        #from cStringIO import StringIO
        #out = StringIO()
        #sh.call(cmd, stdout=out)
        #if not dry_run:
        #    out.seek(0)
        #    job_id = out.read().rstrip()
        #    logging.info('%s (job_id: %s)', cmd, job_id)

    def _get_script_body(self, args):
        if isinstance(args, string_types):
            script_body = args
        else:
            script_body = ' '.join(map(sh.quote, args))
        return script_body

    def call(self, args, name=None, job_dir=None, slot=None, memory=None, hold=None, **kwds):
        script_body = self._get_script_body(args)
        dry_run = kwds.get('dry_run')
        info = self._create_job_script(script_body, name=name, hold=hold, job_dir=job_dir, slot=slot, memory=memory, dry_run=dry_run)
        script = info['script_path']
        return self._qsub_sync(script, dry_run=dry_run)

    def array_call(self, args, ntask, start=1, step=None, max_ntask=None, name=None, job_dir=None, slot=None, memory=None, hold=None, **kwds):
        script_body = self._get_script_body(args)

        assert ntask > 0, 'task number should be > 0'
        assert start > 0, 'start should be > 0'
        dry_run = kwds.get('dry_run')
        info = self._create_job_script(script_body, name=name, hold=hold, job_dir=job_dir, slot=slot, memory=memory, dry_run=dry_run, is_array=True)
        script = info['script_path']

        opts = []
        if step > 0:
            opts.extend(['-t', '{0}-{1}:{2}'.format(start, ntask, step)])
        else:
            opts.extend(['-t', '{0}-{1}'.format(start, ntask)])

        if max_ntask is not None:
            opts.extend(['-tc', max_ntask])

        return self._qsub_sync(script, dry_run=dry_run, ext_opts=opts)
