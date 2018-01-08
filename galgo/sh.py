import subprocess
import logging
import sys
import shlex
from six import string_types
try:
    from shlex import quote
except ImportError:
    from pipes import quote

CalledProcessError = subprocess.CalledProcessError
dry_run = False
default_shell = '/bin/bash'
default_shell_mode = True

# sync api
def call(args, **kwds):
    stdout = kwds.pop('stdout', sys.stdout)
    _dry_run = kwds.get('dry_run', dry_run)
    if not _dry_run and isinstance(stdout, string_types):
        with open(stdout, 'wb+', buffering=0) as stdout_fp:
            logging.debug('Opend for stdout: %s', stdout)
            ret = call(args, stdout=stdout_fp, **kwds)
        logging.debug('Closed for stdout: %s', stdout)
        return ret

    stderr = kwds.pop('stderr', sys.stderr)
    #stderr = kwds.pop('stderr', subprocess.PIPE)
    if not _dry_run and isinstance(stderr, string_types):
        with open(stderr, 'wb+', buffering=0) as stderr_fp:
            logging.debug('Opened for stderr: %s', stderr)
            ret = call(args, stderr=stderr_fp, **kwds)
        logging.debug('Closed for stderr: %s', stderr)
        return ret

    kwds.setdefault('executable', default_shell)
    kwds.setdefault('shell', default_shell_mode)
    if isinstance(args, string_types):
        cmd = args
    else:
        args = map(str, args)
        cmd = ' '.join(map(quote, args))
        #args = shlex.split(args[0])
    logging.info('$ %s', cmd)
    #logging.debug('$ %s', ' '.join(map(quote, args)))
    if not _dry_run:
        kwds.pop('dry_run', None)

        # adhoc
        # stdout_is_buffer = not hasattr(stdout, 'fileno')
        # stderr_is_buffer = not hasattr(stderr, 'fileno')
        # if stdout_is_buffer:
        #     stdout_buffer = stdout
        #     stdout = subprocess.PIPE
        # if stderr_is_buffer:
        #     stderr_buffer = stderr
        #     stderr = subprocess.PIPE
        p = subprocess.Popen(cmd, stdout=stdout, stderr=stderr, **kwds)
        logging.info('pid: %s starts', p.pid)
        # if stdout_is_buffer:
        #     stdout_buffer.write(p.stdout.read())
        # if stderr_is_buffer:
        #     stderr_buffer.write(p.stderr.read())
        ret_code = p.wait()
        logging.info('pid: %s ends with returncode: %s', p.pid, ret_code)
        if ret_code != 0:
            raise CalledProcessError(ret_code, cmd)
    #return p  # TODO change to return code ?
