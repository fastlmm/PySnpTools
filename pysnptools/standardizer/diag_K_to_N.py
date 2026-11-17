import pysnptools.util as pstutil
import numpy as np
from pysnptools.standardizer import Standardizer
import logging
import warnings


class DiagKtoN(Standardizer):
    """
    Both a :class:`.Standardizer` and a :class:`.KernelStandardizer`.

    When applied to a :class:`.SnpData`, it multiplies the SNP values by the square root of a factor such that a kernel
    constructed from the SNP values will have a diagonal that sums to iid_count. This class thus
    standardizes the kernel before it is even constructed.

    When applied to a :class:`.KernelData`, it multiplies the kernel values by a factor such that a kernel
    will have a diagonal that sums to iid_count.

    Scaling follows the values' existing NumPy or CuPy device. This operation
    rescales the whole matrix; it does not center each SNP. Use it on data with
    a finite, positive kernel trace. Zero, negative, or nonfinite normalization
    denominators, including empty inputs, raise ValueError before values change.
    See :doc:`/gpu` for a complete GPU example
    and :class:`.Standardizer` for more about standardization.

    **Constructor:**
        :Parameters: * **deprecated_iid_count** (*int*) (optional) -- Deprecated.


    Example of DiagKtoN to :class:`.SnpData`:

    >>> import numpy as np
    >>> from pysnptools.standardizer import DiagKtoN, Unit, Identity
    >>> from pysnptools.snpreader import Bed
    >>> from pysnptools.util import example_file # Download and return local file name
    >>> bedfile = example_file("tests/datasets/all_chr.maf0.001.N300.*","*.bed")
    >>> snpdata1 = Bed(bedfile,count_A1=False).read().standardize(Unit()).standardize(DiagKtoN())
    >>> kernel1 = snpdata1.read_kernel(Identity())
    >>> print('{0:.6f}'.format(np.diag(kernel1.val).sum()))
    300.000000

    Example of DiagKtoN to :class:`.KernelData`:

    >>> import numpy as np
    >>> from pysnptools.standardizer import DiagKtoN, Unit, Identity
    >>> from pysnptools.snpreader import Bed
    >>> from pysnptools.util import example_file # Download and return local file name
    >>> bedfile = example_file("tests/datasets/all_chr.maf0.001.N300.*","*.bed")
    >>> snpdata1 = Bed(bedfile,count_A1=False).read().standardize(Unit())
    >>> kernel1 = snpdata1.read_kernel(DiagKtoN(),block_size=None)
    >>> print('{0:.6f}'.format(np.diag(kernel1.val).sum()))
    300.000000

    """

    """diag(K)=N standardization of the data"""

    def __init__(self, deprecated_iid_count=None):
        super(DiagKtoN, self).__init__()
        if deprecated_iid_count is not None:
            warnings.warn(
                "'iid_count' is deprecated (and not needed, since can get iid_count from SNPs val's first dimension",
                DeprecationWarning,
            )

    def _standardize_kernel(
        self,
        kerneldata,
        return_trained=False,
        force_python_only=False,
        num_threads=None,
    ):
        xp = pstutil.get_array_module(kerneldata.val)
        trace = (np.diag(kerneldata.val).sum() if xp is np
                 else float(xp.diag(kerneldata.val).sum(dtype=np.float64)))
        factor = _normalization_factor(kerneldata.iid_count, trace)

        if abs(factor - 1.0) > 1e-15:
            kerneldata._val *= factor

        if return_trained:
            return kerneldata, DiagKtoNTrained(factor)
        else:
            return kerneldata

    def standardize(
        self,
        input,
        block_size=None,
        return_trained=False,
        force_python_only=False,
        num_threads=None,
    ):
        from pysnptools.kernelreader import KernelReader

        if block_size is not None:
            warnings.warn(
                "block_size is deprecated (and not needed, since standardization is in-place",
                DeprecationWarning,
            )
        if isinstance(input, KernelReader) and hasattr(input, "val"):
            return self._standardize_kernel(
                input,
                return_trained=return_trained,
                force_python_only=force_python_only,
                num_threads=num_threads,
            )
        else:
            return self._standardize_snps(
                input,
                return_trained=return_trained,
                force_python_only=force_python_only,
                num_threads=num_threads,
            )

    def _standardize_snps(
        self, snps, return_trained=False, force_python_only=False, num_threads=None
    ):
        if hasattr(snps, "val"):
            val = snps.val
        else:
            warnings.warn(
                "standardizing an nparray instead of a SnpData is deprecated",
                DeprecationWarning,
            )
            val = snps

        vec = val.reshape(-1, order="A")
        # make sure no copy was made
        assert not vec.flags["OWNDATA"]
        # Accumulate in float64 without a matrix-sized float64 temporary.
        xp = pstutil.get_array_module(val)
        if xp is np:
            squared_sum = vec.dot(vec)
        else:
            from pysnptools.standardizer._cupy import _sum_squares_kernel
            squared_sum = float(_sum_squares_kernel()(vec))
        factor = _normalization_factor(val.shape[0], squared_sum)

        if abs(factor - 1.0) > 1e-15:
            val *= np.sqrt(factor)

        if return_trained:
            return snps, DiagKtoNTrained(factor)
        else:
            return snps

    def __repr__(self):
        return "{0}()".format(self.__class__.__name__)


# !!!move to its own file
class DiagKtoNTrained(Standardizer):
    """Both a :class:`.Standardizer` and a :class:`.KernelStandardizer`.

    When applied to a :class:`.SnpData`, DiagKtoN standardizes one set of SNP data based on another set of SNP data.

    When applied to a :class:`.KernelData`, DiagKtoN standardizes one kernel data based on another kernel data.

    Reuses the fitted scalar on the values' existing NumPy or CuPy device,
    without refitting it from the target data. The target kernel's trace need
    not equal its sample count. See :class:`.DiagKtoN` for a fitting example.

    **Constructor:**
        :Parameters: * **factor** (*float*) -- Fitted kernel scale. SNP values are multiplied by its square root.

    """

    def __init__(self, factor):
        super(DiagKtoNTrained, self).__init__()
        self.factor = factor

    def standardize(
        self,
        input,
        block_size=None,
        return_trained=False,
        force_python_only=False,
        num_threads=None,
    ):
        if block_size is not None:
            warnings.warn(
                "block_size is deprecated (and not needed, since standardization is in-place",
                DeprecationWarning,
            )

        from pysnptools.kernelreader import KernelReader

        if isinstance(input, KernelReader) and hasattr(input, "val"):
            return self._standardize_kernel(
                input,
                return_trained=return_trained,
                force_python_only=force_python_only,
                num_threads=num_threads,
            )
        else:
            return self._standardize_snps(
                input,
                return_trained=return_trained,
                force_python_only=force_python_only,
                num_threads=num_threads,
            )  # LATER test coverage

    @property
    def is_constant(self):
        return abs(self.factor - 1.0) < 1e-15

    def _standardize_snps(
        self, snps, return_trained=False, force_python_only=False, num_threads=None
    ):
        if hasattr(snps, "val"):
            val = snps.val  # LATER test coverage
        else:
            warnings.warn(
                "standardizing an nparray instead of a SnpData is deprecated",
                DeprecationWarning,
            )  # LATER test coverage
            val = snps

        if not self.is_constant:
            val *= np.sqrt(self.factor)

        if return_trained:
            return snps, self
        else:
            return snps

    def _standardize_kernel(
        self,
        kerneldata,
        return_trained=False,
        force_python_only=False,
        num_threads=None,
    ):
        if not self.is_constant:
            kerneldata._val *= self.factor

        if return_trained:
            return kerneldata, self
        else:
            return kerneldata

    def __repr__(self):
        return "{0}({1})".format(self.__class__.__name__, self.factor)


def _normalization_factor(sample_count, denominator):
    # Validate the existing reduction result: no extra matrix scan or transfer.
    if sample_count == 0 or not np.isfinite(denominator) or denominator <= 0:
        raise ValueError("DiagKtoN requires a finite, positive normalization denominator and at least one sample")
    # Preserve the CPU reduction's scalar dtype and its established rounding.
    factor = float(sample_count) / denominator
    if not np.isfinite(factor) or factor <= 0:
        raise ValueError("DiagKtoN normalization factor must be finite and positive")
    return factor


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    import doctest

    doctest.testmod(optionflags=doctest.ELLIPSIS)
