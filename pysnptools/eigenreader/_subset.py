from pysnptools.eigenreader import EigenReader
from pysnptools.pstreader._subset import _PstSubset

class _EigenSubset(_PstSubset,EigenReader):
    def __init__(self, *args, **kwargs):
        super(_EigenSubset, self).__init__(*args, **kwargs)
