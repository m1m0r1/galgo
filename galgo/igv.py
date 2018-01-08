from __future__ import print_function, absolute_import
import socket
import os
import sys
import signal
import logging
from .utils import get_empty_port


class IGVServer(object):
    def __init__(self, igv_cmd='igv.sh', port=None, outdir=None):
        import subprocess
        from threading import Thread
        import time

        if port is None:
            port = get_empty_port(default=60151)

        self.port = port
        def readit(ffrom, fto, self):
            for line in iter(ffrom.readline, ''):
                if "Listening on port" in line:
                    self._wait = False
                if line.rstrip():
                    fto.write(line)
            ffrom.close()

        if self.port:
            igv_cmd = ' '.join((igv_cmd, '--port', str(self.port)))
        logging.info('igv_cmd: %s', igv_cmd)
        self._proc = p = subprocess.Popen(igv_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        self._wait = True
        self._tout = tout = Thread(target=readit, args=(p.stdout, sys.stderr, self))
        self._terr = terr = Thread(target=readit, args=(p.stderr, sys.stderr, self))
        #tout.daemon = terr.daemon = True
        tout.start()
        terr.start()
        while p.poll() is None and self._wait:
            time.sleep(10)
            logging.info("waiting...")
        logging.info('port: %s', self.port)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()

    def close(self):
        self._proc.send_signal(signal.SIGINT)
        while self._tout.is_alive():
            self._tout.join(1)
        while self._terr.is_alive():
            self._terr.join(1)

    def get_client(self):
        return IGVClient(port=self.port)


class IGVClient(object):
    r"""
    Simple wrapper to the IGV (http://www.broadinstitute.org/software/igv/home)
    socket interface (http://www.broadinstitute.org/software/igv/PortCommands)
    Successful commands return 'OK'
    """
    _socket = None
    _outdir = None

    def __init__(self, host=None, port=None, outdir=None):
        self.host = host or '127.0.0.1'
        self.port = port or 60151
        self.commands = []
        self._connect()
        if outdir:
            self._set_outdir(outdir)

    def _connect(self):
        if self._socket:
            self._socket.close()
        self._socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        self._socket.connect((self.host, self.port))

    def close(self):
        self._socket.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()

    def _set_outdir(self, outdir):
        outdir = os.path.abspath(outdir)
        if not os.path.exists(outdir):
            os.makedirs(outdir)
            logging.info('Created %s', outdir)

        self.send('snapshotDirectory %s' % outdir)
        self._outdir = outdir

    def send(self, cmd):
        logging.info('cmd: %s', cmd)
        self.commands.append(cmd)
        self._socket.send(cmd + '\n')
        return self._socket.recv(4096).rstrip('\n')

    def save(self, path=None):
        """ .png, .jpg, or .svg
        """
        if self._outdir is None:
            self._set_outdir(os.curdir)
        if path is not None:
            # igv assumes the path is just a single filename, but
            # we can set the snapshot dir. then just use the filename.
            dirname = os.path.dirname(path)
            dirname = os.path.abspath(dirname)
            org_dir = self._outdir
            if dirname != org_dir:
                self._set_outdir(dirname)
            status = self.send('snapshot ' + os.path.basename(path))
            if org_dir is not None and org_dir != dirname:
                self._set_outdir(org_dir)
            return status
        else:
            return self.send('snapshot')

    snapshot = save

    def new(self):
        return self.send('new')

    def genome(self, name):
        return self.send('genome {0}'.format(name))

    def load(self, url):
        return self.send('load {0}'.format(url))

    def goto(self, *loci):
        if not loci:
            loci = ['all']
        return self.send('goto {0}'.format(' '.join(map(str, loci))))

    def set_height(self, height=1000):
        return self.send('maxPanelHeight {0}'.format(height))

    def sort(self, option='base', locus=''):
        """
        options is one of: base, position, strand, quality, sample, and
        readGroup.
        base, position, strand, quality, sample, readGroup, AMPLIFICATION, DELETION, EXPRESSION, SCORE, MUTATION_COUNT
        """
        assert option in ("base", "position", "strand", "quality", "sample",
                         "readGroup")
        return self.send('sort {0} {1}'.format(option, locus))

    def expand(self, track=''):
        self.send('expand {0}'.format(track))

    def collapse(self, track=''):
        self.send('collapse {0}'.format(track))

    def squish(self, track=''):
        self.send('squish {0}'.format(track))

    def view_as_pairs(self, track=''):
        self.send('viewaspairs {0}'.format(track))

    def mark_region(self, chrom, start, end):
        self.send('region {0} {1} {2}'.format(chrom, start, end))
