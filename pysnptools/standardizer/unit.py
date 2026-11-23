import numpy as np
import logging
from pysnptools.standardizer import Standardizer
from pysnptools.standardizer.unittrained import UnitTrained
import warnings
import pysnptools.util as pstutil

# TODO Remove the native constant-column caveat in this docstring once PySnpTools requires a bed-reader release whose CPU standardizer zeros constant fractional SNPs.
class Unit(Standardizer):
    """A :class:`.Standardizer` to unit standardize SNP data. Observed nonconstant values are centered and scaled to standard deviation one.
    Observed values use population standard deviation (ddof=0); NaNs then become
    zero, so final variance can be below one. Constant columns become zero;
    wholly missing training columns raise ValueError. NumPy/CuPy arrays stay on
    their device, regardless of ARRAY_MODULE. See :doc:`/gpu` for complete examples.

    Native-path caveat: with bed-reader 1.3.0b2, the default NumPy path
    (``force_python_only=False``) does not recognize a constant column of
    fractional dosages, such as all 0.1. It leaves small nonzero values (up to
    about 3e-3 in float32) and trains a tiny finite standard deviation instead
    of infinity, so applying the trained standardizer to new values of that SNP
    gives nonzero, often very large, results. Integer genotypes, CuPy arrays, and
    ``force_python_only=True`` are unaffected.

    See :class:`.Standardizer` for more information about standardization.

    >>> from pysnptools.standardizer import Unit
    >>> from pysnptools.snpreader import Bed
    >>> from pysnptools.util import example_file # Download and return local file name
    >>> bedfile = example_file("tests/datasets/all_chr.maf0.001.N300.*","*.bed")
    >>> snpdata1 = Bed(bedfile,count_A1=False).read().standardize(Unit())
    >>> print('{0:.6f}'.format(snpdata1.val[0,0]))
    0.229416
    """
    def __init__(self):
        super(Unit, self).__init__()

    def __repr__(self):
        return "{0}()".format(self.__class__.__name__)

    def standardize(self, snps, block_size=None, return_trained=False, force_python_only=False, num_threads=None):
        '''
        Standardize on the array's current device, in place.

        See :meth:`.Standardizer.standardize` for shared options and :doc:`/gpu`
        for complete examples. Returns the input, or (input, UnitTrained) when
        return_trained=True. The trained stats are NumPy arrays. Direct
        block_size arguments warn and are ignored; force_python_only and
        num_threads do not change CuPy execution.
        '''

        if block_size is not None:
            warnings.warn("block_size is deprecated (and not needed, since standardization is in-place", DeprecationWarning)

        if hasattr(snps,"val"):
            val = snps.val
        else:
            warnings.warn("standardizing an ndarray instead of a SnpData is deprecated", DeprecationWarning)
            val = snps

        stats = self._standardize_unit_and_beta(val, is_beta=False, a=np.nan, b=np.nan, apply_in_place=True,
                                                use_stats=False,stats=None,num_threads=num_threads, force_python_only=force_python_only)

        if return_trained:
            assert hasattr(snps,"val"), "return_trained=True requires that snps be a SnpData"
            return snps, UnitTrained(snps.sid, pstutil.asnumpy(stats))
        else:
            return snps

    def _merge_trained(self, trained_list):
        sid = np.concatenate([trained.sid for trained in trained_list])
        stats = np.concatenate([trained.stats for trained in trained_list])
        return UnitTrained(sid, stats)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    import doctest
    doctest.testmod(optionflags=doctest.ELLIPSIS)
