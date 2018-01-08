import matplotlib

def init(use='Agg'):
    matplotlib.use(use)


class PdfGrid(object):
    def __init__(self, output, row=4, col=4):
        self.output = output
        self.row = row
        self.col = col
        self._saved = False

    def __enter__(self):
        self.open()
        return iter(self)

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()

    def open(self):
        from matplotlib.backends.backend_pdf import PdfPages
        self._pp = PdfPages(self.output)

    def close(self):
        import matplotlib.pyplot as plt
        if not self._saved:
            self._pp.savefig(plt.gcf())
        self._pp.close()

    def __iter__(self):
        """
        """
        import matplotlib.pyplot as plt
        page = 1
        while 1:
            plt.figure()
            for i in range(self.row * self.col):
                plt.subplot(self.row, self.col, i+1)  # change current axis
                self._saved = False
                yield page, i+1  # offset
            self._pp.savefig(plt.gcf())
            self._saved = True
            page += 1
