"""
cmk update Tools for reading and manipulating SNP distribution data. For each individual and SNP location, it gives three probabilities: P(AA),P(AB),P(BB),
where A and B are the alleles. The probabilities should sum to 1.0. Missing data is represented by three numpy.NaN's.
"""

from pysnptools.eigenreader.eigenreader import EigenReader
from pysnptools.eigenreader.eigendata import EigenData
#from pysnptools.eigenreader.distnpz import DistNpz
#from pysnptools.eigenreader.disthdf5 import DistHdf5
#from pysnptools.eigenreader.distmemmap import DistMemMap
#from pysnptools.eigenreader.bgen import Bgen
#from pysnptools.eigenreader.distgen import DistGen
#from pysnptools.eigenreader._distmergesids import _DistMergeSIDs

