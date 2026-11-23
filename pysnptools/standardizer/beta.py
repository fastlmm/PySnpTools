import pysnptools.util as pstutil
import numpy as np
import logging
import warnings
from pysnptools.standardizer import Standardizer


# TODO Remove the native constant-column caveat in this docstring once PySnpTools requires a bed-reader release whose CPU standardizer zeros constant fractional SNPs.
class Beta(Standardizer):
    """
    A :class:`.Standardizer` to beta standardize SNP data.

    Subtract the observed mean and multiply by the beta density at the minor
    allele frequency; this does not divide by standard deviation. Native and
    CuPy paths zero constant columns, impute NaNs with zero, and reject wholly
    missing training SNPs and means outside [0, 2]. The CuPy path keeps the
    matrix on its device but computes per-SNP beta weights with SciPy on CPU.
    See :class:`.Standardizer` and :doc:`/gpu` for options and complete examples.

    Native-path caveat: with bed-reader 1.3.0b2, the default NumPy path
    (``force_python_only=False``) does not recognize a constant column of
    fractional dosages, such as all 0.1. It leaves small nonzero values (up to
    about 3e-3 in float32) and trains a tiny finite standard deviation instead
    of infinity, so applying the trained standardizer to new values of that SNP
    gives nonzero results. Integer genotypes, CuPy arrays, and
    ``force_python_only=True`` are unaffected.

    **Constructor:**
        :Parameters: * **a** (*float*) -- The *a* parameter of the beta distribution
                     * **b** (*float*) -- The *b* parameter of the beta distribution

    >>> from pysnptools.standardizer import Beta
    >>> from pysnptools.snpreader import Bed
    >>> from pysnptools.util import example_file # Download and return local file name
    >>> bedfile = example_file("tests/datasets/all_chr.maf0.001.N300.*","*.bed")
    >>> snpdata1 = Bed(bedfile,count_A1=False).read().standardize(Beta(1,25))
    >>> print('{0:.6f}'.format(snpdata1.val[0,0]))
    0.680802
    """

    def __init__(self, a, b):
        super(Beta, self).__init__()
        self.a = a
        self.b = b

    def __repr__(self):
        return "{0}(a={1},b={2})".format(self.__class__.__name__, self.a, self.b)

    def standardize(
        self,
        snpdata,
        block_size=None,
        return_trained=False,
        force_python_only=False,
        num_threads=None,
    ):  # !!!later why is the 2nd argument called 'snpdata' here, but 'snps' in unit.py?
        """Apply Beta in place on the current array device.

        See :meth:`.Standardizer.standardize` for options and :doc:`/gpu` for
        examples. Returns the input, or (input, BetaTrained) with NumPy stats when
        return_trained=True. Direct block_size warns and is ignored. CuPy ignores
        force_python_only and num_threads. Constant training columns remain zero,
        including singular endpoint beta densities, on native and Python CPU paths
        (but see the native-path caveat for constant fractional dosages in :class:`.Beta`)."""
        if block_size is not None:
            warnings.warn(
                "block_size is deprecated (and not needed, since standardization is in-place",
                DeprecationWarning,
            )

        if hasattr(snpdata, "val"):
            val = snpdata.val
        else:
            warnings.warn(
                "standardizing an nparray instead of a SnpData is deprecated",
                DeprecationWarning,
            )
            val = snpdata

        stats = self._standardize_unit_and_beta(
            val,
            is_beta=True,
            a=self.a,
            b=self.b,
            apply_in_place=True,
            use_stats=False,
            stats=None,
            num_threads=num_threads,
            force_python_only=force_python_only,
        )
        if return_trained:
            from pysnptools.standardizer import BetaTrained

            assert hasattr(
                snpdata, "val"
            ), "return_trained=True must be used with SnpData"
            return snpdata, BetaTrained(self.a, self.b, snpdata.sid, pstutil.asnumpy(stats))
        else:
            return snpdata

    def _merge_trained(self, trained_list):
        from pysnptools.standardizer import BetaTrained

        sid = np.concatenate([trained.sid for trained in trained_list])
        stats = np.concatenate([trained.stats for trained in trained_list])
        a_set = set([trained.a for trained in trained_list])
        b_set = set([trained.b for trained in trained_list])
        # Raise explicitly, not with assert, so these checks survive python -O.
        if len(a_set) > 1:
            raise AssertionError("Expect all BetaTrained's to have the same 'a'")
        if len(b_set) > 1:
            raise AssertionError("Expect all BetaTrained's to have the same 'b'")
        a = list(a_set)[0] if a_set else None
        b = list(b_set)[0] if b_set else None
        return BetaTrained(a, b, sid, stats)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    import doctest

    doctest.testmod(optionflags=doctest.ELLIPSIS)
