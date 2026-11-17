import numpy as np
from numbers import Integral
from itertools import *  # noqa: F403
import logging
from pysnptools.standardizer import Identity as SS_Identity
from pysnptools.kernelreader import KernelReader
from pysnptools.kernelreader import KernelData
from pysnptools.kernelstandardizer import DiagKtoN
import pysnptools.util as pstutil


class SnpKernel(KernelReader):
    """
    A :class:`.KernelReader` that creates a kernel from a :class:`.SnpReader`. No SNP data will be read until
    the :meth:`SnpKernel.read` method is called. A positive block_size limits SNP
    blocks during kernel construction; the dense kernel and temporary products
    still need memory. Reading a subset can construct the full kernel first.
    See :doc:`/gpu` for explicit read backends, dtype/order options, and examples.

    See :class:`.KernelReader` for general examples of using KernelReaders.

    **Constructor:**
        :Parameters: * **snpreader** (:class:`SnpReader`) -- The SNP data
                     * **standardizer** (:class:`Standardizer`) -- How the SNP data should be standardized
                     * **block_size** (optional, positive int) -- The number of SNPs to read at a time. Invalid values raise ValueError.

        If **block_size** is not given, then all SNP data will be read at once.

        :Example:

        >>> from pysnptools.snpreader import Bed
        >>> from pysnptools.standardizer import Unit
        >>> from pysnptools.util import example_file # Download and return local file name
        >>> bed_file = example_file('pysnptools/examples/toydata.5chrom.*','*.bed')
        >>> snp_on_disk = Bed(bed_file,count_A1=False)     # A Bed file is specified, but nothing is read from disk
        >>> kernel_on_disk = SnpKernel(snp_on_disk, Unit(),block_size=500)  # A kernel is specified, but nothing is read from disk
        >>> print(kernel_on_disk) #Print the specification
        SnpKernel(Bed(...pysnptools/examples/toydata.5chrom.bed',count_A1=False),standardizer=Unit(),block_size=500)
        >>> print(kernel_on_disk.iid_count)                                  # iid information is read from disk, but not SNP data
        500
        >>> kerneldata = kernel_on_disk.read().standardize()                # SNPs are read and Unit standardized, 500 at a time, to create a kernel, which is then standardized
        >>> print('{0:.6f}'.format(kerneldata.val[0,0]))
        0.992307
    """

    def __init__(self, snpreader, standardizer=None, block_size=None):
        super(SnpKernel, self).__init__()

        assert standardizer is not None, "'standardizer' must be provided"
        if block_size is not None and (
            isinstance(block_size, (bool, np.bool_))
            or not isinstance(block_size, Integral)
            or block_size <= 0
        ):
            raise ValueError("block_size must be a positive integer or None")

        self.snpreader = snpreader
        self.standardizer = standardizer
        self.block_size = block_size

    @property
    def row(self):
        return self.snpreader.iid

    @property
    def col(self):
        if hasattr(self, "_col"):
            return self._col
        else:
            return self.snpreader.iid

    def __repr__(self):
        return self._internal_repr(self.standardizer)

    def _internal_repr(self, standardizer):  # !!! merge this with __repr__
        s = "SnpKernel({0},standardizer={1}".format(self.snpreader, standardizer)
        if self.block_size is not None:
            s += ",block_size={0}".format(self.block_size)
        s += ")"
        return s

    def copyinputs(self, copier):
        # Doesn't need run_once
        copier.input(self.snpreader)
        copier.input(self.standardizer)

    def _read(
        self,
        row_index_or_none,
        col_index_or_none,
        order,
        dtype,
        force_python_only,
        view_ok,
        num_threads,
    ):
        return self._read_array(
            row_index_or_none, col_index_or_none, order, dtype,
            force_python_only, view_ok, num_threads, pstutil.array_module())

    def _read_array(self, row_index, col_index, order, dtype,
                    force_python_only, view_ok, num_threads, xp):
        if (self.standardizer.is_constant and row_index is not None
                and col_index is not None and np.array_equal(row_index, col_index)):
            return self.snpreader[row_index, :]._read_kernel(
                self.standardizer, block_size=self.block_size, order=order, dtype=dtype,
                force_python_only=force_python_only, view_ok=view_ok,
                num_threads=num_threads, xp=xp)
        # LATER: Building just the requested rectangle could save memory.
        whole = self.snpreader._read_kernel(
            self.standardizer, block_size=self.block_size, order=order, dtype=dtype,
            force_python_only=force_python_only, view_ok=view_ok,
            num_threads=num_threads, xp=xp)
        val, shares_memory = self._apply_sparray_or_slice_to_val(
            whole, row_index, col_index, order, dtype, force_python_only, num_threads)
        return val

    def __getitem__(self, iid_indexer_and_snp_indexer):
        if isinstance(iid_indexer_and_snp_indexer, tuple):
            row_index_or_none, col_index_or_none = iid_indexer_and_snp_indexer
        else:
            row_index_or_none = iid_indexer_and_snp_indexer
            col_index_or_none = row_index_or_none

        # Special case: If square and constant, can push the subsetting into the SnpReader
        if (
            self.standardizer.is_constant
            and row_index_or_none is not None
            and col_index_or_none is not None
            and np.array_equal(row_index_or_none, col_index_or_none)
        ):
            return SnpKernel(
                self.snpreader[row_index_or_none, :],
                self.standardizer,
                block_size=self.block_size,
            )
        else:
            return KernelReader.__getitem__(self, iid_indexer_and_snp_indexer)

    def _read_with_standardizing(
        self,
        to_kerneldata,
        kernel_standardizer=DiagKtoN(),
        return_trained=False,
        num_threads=None,
        dtype=None,
        xp=None,
    ):
        """
        Reads a SnpKernel with two cases
              If returning KernelData,
                 just calls snpreader._read_kernel, package, kernel_standardize
              If returning simple SnpKernel that needs no more standardization
                  read the reference and learn both standardization (but can't this cause multiple reads?)
        Note that snp_standardizer should be None or the standardizer instead the SnpKernel should have the placeholder value Standardizer()

        Will choose which array module to use based on the ARRAY_MODULE environment variable (e.g. 'numpy' (default) or 'cupy')

        """
        logging.info("Starting '_read_with_standardizing'")
        xp = pstutil.array_module(xp)
        dtype = pstutil.array_dtype(dtype)

        if to_kerneldata:
            val, snp_trained = self.snpreader._read_kernel(
                self.standardizer,
                block_size=self.block_size,
                dtype=dtype, xp=xp,
                return_trained=True,
                num_threads=num_threads,
            )
            kernel = KernelData(iid=self.snpreader.iid, val=val, name=str(self), xp=xp)
            kernel, kernel_trained = kernel.standardize(
                kernel_standardizer, return_trained=True, num_threads=num_threads
            )
        else:
            snpdata, snp_trained = self.snpreader._read_with_standardizing(
                self.standardizer, return_trained=True, num_threads=num_threads,
                dtype=dtype, xp=xp
            )
            snpdata, kernel_trained = snpdata.standardize(
                kernel_standardizer, return_trained=True, num_threads=num_threads
            )
            kernel = SnpKernel(snpdata, SS_Identity())

        logging.info("Ending '_read_with_standardizing'")
        if return_trained:
            return kernel, snp_trained, kernel_trained
        else:
            return kernel

    @property
    def sid(self):
        """The :attr:`SnpReader.sid` property of the SNP data."""
        return self.snpreader.sid

    @property
    def sid_count(self):
        """The :attr:`SnpReader.sid_count` property of the SNP data."""
        return self.snpreader.sid_count

    @property
    def pos(self):
        """The :attr:`SnpReader.pos` property of the SNP data."""
        return self.snpreader.pos

    # LATER this needs a test
    def read_snps(
        self,
        order="F",
        dtype=None,
        force_python_only=False,
        view_ok=False,
        num_threads=None,
        xp=None,
    ):
        """Reads the SNP values, applies the standardizer, and returns a :class:`.SnpData`.

        Values are standardized on the selected device. See :doc:`/gpu` for
        GPU prerequisites, training semantics, and complete workflow examples.

        :param order: {'F' (default), 'C', 'A'}, optional -- Specify the order of the ndarray. If order is 'F' (default),
            then the array will be in F-contiguous order (iid-index varies the fastest).
            If order is 'C', then the returned array will be in C-contiguous order (sid-index varies the fastest).
            If order is 'A', then the :attr:`SnpData.val`
            ndarray may be in any order (either C-, Fortran-contiguous).
        :type order: string or None

        :param dtype: {None, numpy.float64, numpy.float32}, optional -- The data-type for the :attr:`SnpData.val` ndarray.
             None uses ARRAY_DTYPE, then float64; explicit dtype overrides it.
        :type dtype: data-type

        :param force_python_only: Request a Python CPU path where supported.
            Native GPU BED reads reject True; this flag does not select a CPU backend.
        :type force_python_only: bool

        :param view_ok: Accepted for compatibility but ignored. This method reads
            independent values before standardizing them in place.
        :type view_ok: bool

        :param num_threads: optional -- The CPU thread count, ignored by GPU decoding and CuPy math. Defaults to all available
            processors. Can also be set with these environment variables (listed in priority order):
            'PST_NUM_THREADS', 'NUM_THREADS', 'MKL_NUM_THREADS'.
        :type num_threads: None or int

        :param xp: NumPy or CuPy module or name. An explicit choice overrides
            ARRAY_MODULE; None uses the environment, then NumPy.
        :type xp: module, str, or None

        :rtype: :class:`.SnpData`

        >>> from pysnptools.kernelreader import SnpKernel
        >>> from pysnptools.snpreader import SnpData
        >>> from pysnptools.standardizer import Unit
        >>> source = SnpData(iid=[["f", "0"], ["f", "1"]], sid=["s"], val=[[0.], [2.]], xp="numpy")
        >>> standardized = SnpKernel(source, Unit()).read_snps(xp="numpy")
        >>> standardized.val.tolist()
        [[-1.0], [1.0]]
        """
        dtype = pstutil.array_dtype(dtype)
        return self.snpreader._read_with_standardizing(
            self.standardizer, order=order, dtype=dtype, xp=xp,
            force_python_only=force_python_only, num_threads=num_threads)



if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    import doctest

    doctest.testmod(optionflags=doctest.ELLIPSIS)
    # There is also a unit test case in 'pysnptools\test.py' that calls this doc test
