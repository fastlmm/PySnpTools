"""
Tools for reading and manipulating SNP distribution data. For each individual and SNP location, it gives three probabilities: P(AA),P(AB),P(BB),
where A and B are the alleles. The probabilities should sum to 1.0. Missing data is represented by three numpy.NaN's.
"""

from pysnptools.distreader.distreader import DistReader  # noqa: F401, E402
from pysnptools.distreader.distdata import DistData  # noqa: F401, E402
from pysnptools.distreader.distnpz import DistNpz  # noqa: F401, E402
from pysnptools.distreader.disthdf5 import DistHdf5  # noqa: F401, E402
from pysnptools.distreader.distmemmap import DistMemMap  # noqa: F401, E402
from pysnptools.distreader.bgen import Bgen  # noqa: F401, E402
from pysnptools.distreader.distgen import DistGen  # noqa: F401, E402
from pysnptools.distreader._distmergesids import _DistMergeSIDs  # noqa: F401, E402

