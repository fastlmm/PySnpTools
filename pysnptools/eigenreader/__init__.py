"""
cmk update Tools for reading and manipulating SNP distribution data. For each individual and SNP location, it gives three probabilities: P(AA),P(AB),P(BB),
where A and B are the alleles. The probabilities should sum to 1.0. Missing data is represented by three numpy.NaN's.
"""

from pysnptools.eigenreader.eigenreader import EigenReader
from pysnptools.eigenreader.eigendata import EigenData
from pysnptools.eigenreader.eigennpz import EigenNpz
from pysnptools.eigenreader.eigenmemmap import EigenMemMap

#!!!cmk make sure all doc tests get run