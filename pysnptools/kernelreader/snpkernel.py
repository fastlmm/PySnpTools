import numpy as np
from itertools import *
import pandas as pd
import logging
from pysnptools.standardizer import Unit, Standardizer, UnitTrained, BetaTrained
from pysnptools.standardizer import Identity as SS_Identity
from pysnptools.kernelreader import KernelReader
from pysnptools.kernelreader import KernelData
from pysnptools.kernelstandardizer import DiagKtoN
import pysnptools.util as pstutil

class SnpKernel(KernelReader):
    '''
    A :class:`.KernelReader` that creates a kernel from a :class:`.SnpReader`. No SNP data will be read until
    the :meth:`SnpKernel.read` method is called. Use block_size to avoid ever reading all the SNP data into memory.
    at once.

    See :class:`.KernelReader` for general examples of using KernelReaders.

    **Constructor:**
        :Parameters: * **snpreader** (:class:`SnpReader`) -- The SNP data
                     * **standardizer** (:class:`Standardizer`) -- How the SNP data should be standardized
                     * **block_size** (optional, int) -- The number of SNPs to read at a time.

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
    '''
    def __init__(self, snpreader, standardizer=None, block_size=None):
        super(SnpKernel, self).__init__()

        assert standardizer is not None, "'standardizer' must be provided"

        self.snpreader = snpreader
        self.standardizer = standardizer
        self.block_size = block_size

    @property
    def row(self):
        return self.snpreader.iid

    @property
    def col(self):
        if hasattr(self,'_col'):
            return self._col
        else:
            return self.snpreader.iid

    def __repr__(self):
        return self._internal_repr(self.standardizer)

    def _internal_repr(self,standardizer): #!!! merge this with __repr__
        s = "SnpKernel({0},standardizer={1}".format(self.snpreader,standardizer)
        if self.block_size is not None:
            s += ",block_size={0}".format(self.block_size)
        s += ")"
        return s

    def copyinputs(self, copier):
        #Doesn't need run_once
        copier.input(self.snpreader)
        copier.input(self.standardizer)

    def _read(self, row_index_or_none, col_index_or_none, order, dtype, force_python_only, view_ok, num_threads):
        dtype = np.dtype(dtype)
        #Special case: If square and constant, can push the subsetting into the SnpReader
        if (self.standardizer.is_constant and row_index_or_none is not None and col_index_or_none is not None and np.array_equal(row_index_or_none,col_index_or_none)):
            return self.snpreader[row_index_or_none,:]._read_kernel(self.standardizer,self.block_size,order, dtype, force_python_only, view_ok, num_threads)#!!!LATER can this ever be reached?
        else:
            #LATER: If it was often that case that we wanted to standardize on all the data, but then only return a slice of the result,
            #       that could be done with less memory by working in blocks but not tabulating for all the iids.
            whole = self.snpreader._read_kernel(self.standardizer,self.block_size,order, dtype, force_python_only, view_ok, num_threads)
            val, shares_memory = self._apply_sparray_or_slice_to_val(whole, row_index_or_none, col_index_or_none, order, dtype, force_python_only, num_threads)
            return val

    def __getitem__(self, iid_indexer_and_snp_indexer):
        if isinstance(iid_indexer_and_snp_indexer,tuple):
            row_index_or_none, col_index_or_none = iid_indexer_and_snp_indexer
        else:
            row_index_or_none = iid_indexer_and_snp_indexer
            col_index_or_none = row_index_or_none

        #Special case: If square and constant, can push the subsetting into the SnpReader
        if (self.standardizer.is_constant and row_index_or_none is not None and col_index_or_none is not None and np.array_equal(row_index_or_none,col_index_or_none)):
            return SnpKernel(self.snpreader[row_index_or_none,:], self.standardizer, block_size=self.block_size)
        else:
            return KernelReader.__getitem__(self,iid_indexer_and_snp_indexer)


    def _read_with_standardizing(self, to_kerneldata, kernel_standardizer=DiagKtoN(), return_trained=False, num_threads=None):
        '''
        Reads a SnpKernel with two cases
              If returning KernelData,
                 just calls snpreader._read_kernel, package, kernel_standardize
              If returning simple SnpKernel that needs no more standardization
                  read the reference and learn both standardization (but can't this cause multiple reads?)
        Note that snp_standardizer should be None or the standardizer instead the SnpKernel should have the placeholder value Standardizer()

        Will choose which array module to use based on the ARRAY_MODULE environment variable (e.g. 'numpy' (default) or 'cupy')

        '''
        logging.info("Starting '_read_with_standardizing'")
        xp = pstutil.array_module()

        if to_kerneldata:
            val, snp_trained = self.snpreader._read_kernel(self.standardizer,block_size=self.block_size,return_trained=True,num_threads=num_threads)
            kernel = KernelData(iid=self.snpreader.iid, val=val, name=str(self), xp=xp)
            kernel, kernel_trained = kernel.standardize(kernel_standardizer,return_trained=True,num_threads=num_threads)
        else:
            snpdata, snp_trained = self.snpreader.read().standardize(self.standardizer, return_trained=True, num_threads=num_threads)
            snpdata, kernel_trained = snpdata.standardize(kernel_standardizer, return_trained=True, num_threads=num_threads)
            kernel = SnpKernel(snpdata, SS_Identity())

        logging.info("Ending '_read_with_standardizing'")
        if return_trained:
            return kernel, snp_trained, kernel_trained
        else:
            return kernel

    @property
    def sid(self):
        '''The :attr:`SnpReader.sid` property of the SNP data.
        '''
        return self.snpreader.sid

    @property
    def sid_count(self):
        '''The :attr:`SnpReader.sid_count` property of the SNP data.
        '''
        return self.snpreader.sid_count

    @property
    def pos(self):
        '''The :attr:`SnpReader.pos` property of the SNP data.
        '''
        return self.snpreader.pos

    #LATER this needs a test
    def read_snps(self, order='F', dtype=np.float64, force_python_only=False, view_ok=False, num_threads=None):
        """Reads the SNP values, applies the standardizer, and returns a :class:`.SnpData`.

        :param order: {'F' (default), 'C', 'A'}, optional -- Specify the order of the ndarray. If order is 'F' (default),
            then the array will be in F-contiguous order (iid-index varies the fastest).
            If order is 'C', then the returned array will be in C-contiguous order (sid-index varies the fastest).
            If order is 'A', then the :attr:`SnpData.val`
            ndarray may be in any order (either C-, Fortran-contiguous).
        :type order: string or None

        :param dtype: {numpy.float64 (default), numpy.float32}, optional -- The data-type for the :attr:`SnpData.val` ndarray.
        :type dtype: data-type

        :param force_python_only: optional -- If False (default), may use outside library code. If True, requests that the read
            be done without outside library code.
        :type force_python_only: bool

        :param view_ok: optional -- If False (default), allocates new memory for the :attr:`SnpData.val`'s ndarray. If True,
            if practical and reading from a :class:`SnpData`, will return a new 
            :class:`SnpData` with a ndarray shares memory with the original :class:`SnpData`.
            Typically, you'll also wish to use "order='A'" to increase the chance that sharing will be possible.
            Use these parameters with care because any change to either ndarray (for example, via :meth:`.SnpData.standardize`) will effect
            the others. Also keep in mind that :meth:`read` relies on ndarray's mechanisms to decide whether to actually
            share memory and so it may ignore your suggestion and allocate a new ndarray anyway.
        :type view_ok: bool

        :param num_threads: optional -- The number of threads with which to read data. Defaults to all available
            processors. Can also be set with these environment variables (listed in priority order):
            'PST_NUM_THREADS', 'NUM_THREADS', 'MKL_NUM_THREADS'.
        :type num_threads: None or int

        :rtype: :class:`.SnpData`

        """
        dtype = np.dtype(dtype)
        return self.snpreader.read(order=order, dtype=dtype, force_python_only=force_python_only, view_ok=view_ok, num_threads=num_threads).standardize(self.standardizer, num_threads=num_threads)

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    import doctest
    doctest.testmod(optionflags=doctest.ELLIPSIS)
    # There is also a unit test case in 'pysnptools\test.py' that calls this doc test
