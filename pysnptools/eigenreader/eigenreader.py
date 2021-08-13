import numpy as np
#import os
#import os.path
#from itertools import *
#import pandas as pd
import logging
#import time
#import pysnptools.util as pstutil
from pysnptools.pstreader import PstReader
#from pysnptools.snpreader import SnpData
#import warnings
#import pysnptools.standardizer as stdizer
#from pysnptools.snpreader._eigen2snp import _eigen2Snp
#cmk

#!!why do the examples use ../tests/datasets instead of "examples"?
class EigenReader(PstReader):
    """A EigenReader is one of three things:

    * A class such as :class:`.EigenNpz` for you to specify a file with data. For example,

        >>> from pysnptools.eigenreader import EigenNpz
        >>> from pysnptools.util import example_file # Download and return local file name
        >>> npz_file = example_file("pysnptools/examples/toydata.eigen.npz")
        >>> eigen_on_disk = EigenNpz(npz_file)
        >>> print(eigen_on_disk) # prints the name of the file reader
        EigenNpz('...pysnptools/examples/toydata.eigen.npz')
        >>> eigen_on_disk.eid_count # prints the number of eigenvalues and eigenvectors (but doesn't necessarily read any eigenvectors)
        500

    * A :class:`.EigenData` class that holds column eigenvectors in memory, typically after reading it from disk:
        *The eigenvalues, which are typically much smaller, are available even without reading.)

        >>> from pysnptools.util import example_file # Download and return local file name
        >>> npz_file = example_file("pysnptools/examples/toydata.eigen.npz")
        >>> eigen_on_disk = EigenNpz(npz_file)
        >>> eigendata1 = eigen_on_disk.read() #reads the column eigenvectors
        >>> print(type(eigendata1.vectors)) # The vectors property is an 2-D ndarray of column eigenvectors
        <class 'numpy.ndarray'>
        >>> print(eigendata1) # prints the name in-memory eigenvector reader.
        EigenData(EigenNpz('...pysnptools/examples/toydata.eigen.npz'))
        >>> eigendata1.iid_count #prints the number of iids (number of individuals) in this in-memory data
        500

    * A subset of any EigenReader, specified with "[ *iid_index* , *eid_index* ]", to read only some column eigenvectors or
      parts of column vectors. It can also be used to re-order the rows (individuals) or columns (vectors).

        >>> from pysnptools.util import example_file # Download and return local file name
        >>> npz_file = example_file("pysnptools/examples/toydata.eigen.npz")
        >>> eigen_on_disk = EigenNpz(npz_file)
        >>> subset_on_disk = eigen_on_disk[[3,4],::2] # specification for a subset of the data on disk. No eigenvectors are read yet.
        >>> print(subset_on_disk.eid_count) # prints the number of eids in this subset (but still doesn't read any column eigenvectors)
        250
        >>> print(subset_on_disk) #prints a specification of 'subset_on_disk'
        EigenNpz('...pysnptools/examples/toydata.eigen.npz')[[3,4],::2]
        >>> eigendata_subset = subset_on_disk.read() # efficiently reads the specified subset of values from the disk
        >>> print(eigendata_subset) # prints the specification of the in-memory eigenvector information
        EigenData(EigenNpz('...pysnptools/examples/toydata.eigen.npz')[[3,4],::2])
        >>> print((int(eigendata_subset.vectors.shape[0]), int(eigendata_subset.vectors.shape[1]))) # The dimensions of the ndarray of column eigenvectors
        (2, 250)

    The EigenReaders Classes

        ========================= =================== ====================== ================== ======================
        *Class*                   *Format*            *Random Access*        *Suffixes*         *Write* method?
        :class:`.EigenData`               in-memory          Yes                    *n/a*              *n/a*
        :class:`.EigenNpz`                binary             No                     .eigen.npz        Yes
        ========================= =================== ====================== ================== ======================
    
  
    Methods & Properties:

        Every EigenReader, such as :class:`.EigenNpz` and :class:`.EigenData`, has these properties: :attr:`iid`, :attr:`iid_count`, :attr:`eid`, :attr:`eid_count`,
        :attr:`values` and these methods: :meth:`read`, :meth:`iid_to_index`, :meth:`eid_to_index`. See below for details.

        :class:`.EigenData` is a EigenReader so it supports the above properties and methods. In addition, it supports property :attr:`EigenData.values`.

        See below for details.

        Many of the classes, such as :class:`.EigenNpz`, also provide a static :meth:`EigenNpz.write` method for writing :class:`.EigenReader` or :class:`.EigenData` to disk.

        >>> from pysnptools.eigenreader import EigenNpz, EigenData
        >>> from pysnptools.kernelreader import KernelNpz
        >>> import pysnptools.util as pstutil
        >>> from pysnptools.util import example_file #  Download and return local file name
        >>> kernel_npz_file = example_file('pysnptools/examples/toydata.kernel.npz')
        >>> kernel_data = KernelNpz(kernel_npz_file).read()
        >>> values, vectors = np.linalg.eigh(kernel_data.read().val)
        >>> eigendata = EigenData(values=values,vectors=vectors,iid=kernel_data.iid)
        >>> pstutil.create_directory_if_necessary("tempdir/toydata.eigen.npz")
        >>> EigenNpz.write("tempdir/toydata.eigen.npz",eigendata)          #  Write data in EigenNpz format
        EigenNpz('tempdir/toydata.eigen.npz')


    iids and eids:

        Individuals are identified with an iid, which is an ndarray of two strings: a family ID and a individual ID.
        Pairs of eigenvalues and eigenvectors are identified with a (perhaps arbitrary) eid string in an ndarray. For example:

        >>> from pysnptools.eigenreader import EigenNpz
        >>> from pysnptools.util import example_file # Download and return local file name
        >>> npz_file = example_file("pysnptools/examples/toydata.eigen.npz")
        >>> eigen_on_disk = EigenNpz(npz_file)
        >>> print(eigen_on_disk.iid[:3]) # print the first three iids
        [['POP1' '0']
         ['POP1' '12']
         ['POP1' '44']]
        >>> print(eigen_on_disk.eid[:4]) # print the first four eids
        ['eid0', 'eid1', 'eid2', 'eid4']
        >>> print(eigen_on_disk.iid_to_index([['POP1','44'],['POP1','12']])) #Find the indexes for two iids.
        [2 1]

    Selecting and Reordering Individuals and Eigenvectors/Values

        Selecting and reordering individuals is like  :class:`.SnpReader`.
        Selecting and reordering eigenvectors and eigenvalues is similar to :class:`.SnpReader` except instead of selecting and reordering SNPs,
        we select and reorder eigenvectors and eigenvalues.
        
    When Data is Read:

        *same as with* :class:`.SnpReader`

    When Data is Re-Read and Copied:

        *same as with* :class:`.SnpReader`

    Avoiding Unwanted ndarray Allocations

        *same as with* :class:`.SnpReader`

    The :meth:`read` Method
  
        By default the :meth:`read` returns a 2-D ndarray of numpy.float64 laid out in memory in F-contiguous order
        (iid-index varies the fastest). You may, instead,
        ask for numpy.float32 or for C-contiguous order or any order. See :meth:`read` for details.


    Details of Methods & Properties:
    """

    def __init__(self, *args, **kwargs):
        super(EigenReader, self).__init__(*args, **kwargs)

    @property
    def iid(self):
        """An ndarray of the iids. Each iid is an ndarray of two strings (a family ID and a individual ID) that identifies an individual.

        :rtype: ndarray of strings with shape [:attr:`.iid_count`,2]

        This property (to the degree practical) reads only iid and eid data from the disk, not column eigenvectors. Moreover, the iid and eid data is read from file only once.

        :Example:

        >>> from pysnptools.eigenreader import EigenNpz
        >>> from pysnptools.util import example_file # Download and return local file name
        >>> npz_file = example_file("pysnptools/examples/toydata.eigen.npz")
        >>> eigen_on_disk = EigenNpz(npz_file)
        >>> print(eigen_on_disk.iid[:3]) # print the first three iids
        [['0' 'iid_0']
         ['0' 'iid_1']
         ['0' 'iid_2']]
        """
        return self.row

    @property
    def iid_count(self):
        """number of iids

        :rtype: integer

        This property (to the degree practical) reads only iid and eid data from the disk, not column eigenvectors. Moreover, the iid and eid data is read from file only once.
        """
        return self.row_count

    @property
    def eid(self):
        """An ndarray of the eids. Each eid is a string that identifies a eigenvector/value pair.

        :rtype: ndarray (length :attr:`.eid_count`) of strings

        This property (to the degree practical) reads only iid, eid, and eigenvalue data from the disk, not eigenvectors. Moreover, the iid, eid, and eigenvalue data is read from file only once.

        :Example:

        >>> from pysnptools.eigenreader import EigenNpz
        >>> from pysnptools.util import example_file # Download and return local file name
        >>> npz_file = example_file("pysnptools/examples/toydata.eigen.npz")
        >>> eigen_on_disk = EigenNpz(npz_file)
        >>> print(eigen_on_disk.eid[:4]) # print the first four eids
        ['eid_0' 'eid_1' 'eid_2' 'eid_3']
        """
        return self.col

    @property
    def eid_count(self):
        """number of eids which is the number of eigenvector and eigenvalue pairs.

        :rtype: integer

        This property (to the degree practical) reads only iid, eid, eigenvalue data from the disk, not eigenvectors. Moreover, the iid, eid, and eigenvalue data is read from file only once.

        """
        return self.col_count

    @property
    def values(self):
        """An ndarray of the position information for each eid. Each element is an ndarray of three numpy.numbers (chromosome, genetic eigenance, basepair eigenance).

        :rtype: ndarray of float64 with shape [:attr:`.eid_count`, 3]

        This property (to the degree practical) reads only iid and eid data from the disk, not column eigenvectors. Moreover, the iid and eid data is read from file only once.

        :Example:

        >>> from pysnptools.eigenreader import EigenNpz
        >>> from pysnptools.util import example_file # Download and return local file name
        >>> npz_file = example_file("pysnptools/examples/toydata.eigen.npz")
        >>> eigen_on_disk = EigenNpz(npz_file)
        >>> print(eigen_on_disk.pos[:4,].astype('int')) # print position information for the first four eids: #The '...' is for possible space char
        [[        1         0   9270273]
         [        1         0  39900273]
         [        1         0  70530273]
         [        1         0 101160273]]
        """
        return self.col_property

    @property
    def row_property(self):
        """Defined as a zero-width array for compatibility with :class:`PstReader`, but not used.
        """
        if not hasattr(self,'_row_property'):
            self._row_property = np.empty((self.row_count,0))
        return self._row_property


    def _read(self, iid_index_or_none, eid_index_or_none, order, dtype, force_python_only, view_ok, num_threads):
        raise NotImplementedError
    
    #!!check that views always return contiguous memory by default
    def read(self, order='F', dtype=np.float64, force_python_only=False, view_ok=False, num_threads=None):
        """Reads the SNP values and returns a :class:`.EigenData` (with :attr:`EigenData.vectors` property containing a new 2D ndarray of the column eigenvectors).

        :param order: {'F' (default), 'C', 'A'}, optional -- Specify the order of the ndarray. If order is 'F' (default),
            then the array will be in F-contiguous order (iid-index varies the fastest).
            If order is 'C', then the returned array will be in C-contiguous order (eid-index varies the fastest).
            If order is 'A', then the :attr:`EigenData.vectors`
            ndarray may be in any order (either C-, Fortran-contiguous).
        :type order: string or None

        :param dtype: {numpy.float64 (default), numpy.float32}, optional -- The data-type for the :attr:`EigenData.vectors` ndarray.
        :type dtype: data-type

        :param force_python_only: optional -- If False (default), may use outeide library code. If True, requests that the read
            be done without outeide library code.
        :type force_python_only: bool

        :param view_ok: optional -- If False (default), allocates new memory for the :attr:`EigenData.vectors`'s ndarray. If True,
            if practical and reading from a :class:`EigenData`, will return a new 
            :class:`EigenData` with an ndarray shares memory with the original :class:`EigenData`.
            Typically, you'll also wish to use "order='A'" to increase the chance that sharing will be possible.
            Use these parameters with care because any change to either ndarraywill effect
            the others. Also keep in mind that :meth:`read` relies on ndarray's mechanisms to decide whether to actually
            share memory and so it may ignore your suggestion and allocate a new ndarray anyway.
        :type view_ok: bool

        :param num_threads: optional -- The number of threads with which to read data. Defaults to all available
            processors. Can also be set with these environment variables (listed in priority order):
            'PST_NUM_THREADS', 'NUM_THREADS', 'MKL_NUM_THREADS'.
        :type num_threads: None or int

        :rtype: :class:`.EigenData`

        Calling the method again causes the column eigenvectors to be re-read and creates a new in-memory :class:`.EigenData` with a new ndarray of SNP values.

        If you request the values for only a subset of the eids or iids, (to the degree practical) only that subset will be read from disk.

        :Example:

        >>> from pysnptools.eigenreader import EigenNpz
        >>> from pysnptools.util import example_file # Download and return local file name
        >>> npz_file = example_file("pysnptools/examples/toydata.eigen.npz")
        >>> eigen_on_disk = EigenNpz(npz_file) # Specify SNP data on disk
        >>> eigendata1 = eigen_on_disk.read() # Read all the SNP data returning a EigenData instance
        >>> print(type(eigendata1.vectors).__name__) # The EigenData instance contains an ndarray of the data.
        ndarray
        >>> subset_eigendata = eigen_on_disk[:,::2].read() # From the disk, read SNP values for every other eid
        >>> print(subset_eigendata.vectors[0,0]) # Print the first SNP value in the subset
        [0.466804   0.38812848 0.14506752]
        >>> subsub_eigendata = subset_eigendata[:10,:].read(order='A',view_ok=True) # Create an in-memory subset of the subset with SNP values for the first ten iids. Share memory if practical.
        >>> import numpy as np
        >>> # print np.may_share_memory(subset_eigendata.vectors, subsub_eigendata.vectors) # Do the two ndarray's share memory? They could. Currently they won't.       
        """
        dtype = np.dtype(dtype)
        vectors = self._read(None, None, order, dtype, force_python_only, view_ok, num_threads)
        from pysnptools.eigenreader import EigenData
        ret = EigenData(self.iid,self.eid,vectors,values=self.values,name=str(self))
        return ret

    def iid_to_index(self, list):
        """Takes a list of iids and returns a list of index numbers

        :param list: list of iids
        :type order: list of list of strings

        :rtype: ndarray of int
        
        This method (to the degree practical) reads only iid and eid data from the disk, not SNP value data. Moreover, the iid and eid data is read from file only once.

        :Example:

        >>> from pysnptools.eigenreader import EigenNpz
        >>> from pysnptools.util import example_file # Download and return local file name
        >>> npz_file = example_file("pysnptools/examples/toydata.eigen.npz")
        >>> eigen_on_disk = EigenNpz(npz_file) # Specify SNP data on disk
        >>> print(eigen_on_disk.iid_to_index([['0','iid_2'],['0','iid_1']])) #Find the indexes for two iids.
        [2 1]
        """
        return self.row_to_index(list)

    def eid_to_index(self, list):
        """cmkTakes a list of eids and returns a list of index numbers

        :param list: list of eids
        :type list: list of strings

        :rtype: ndarray of int
        
        This method (to the degree practical) reads only iid and eid data from the disk, not SNP value data. Moreover, the iid and eid data is read from file only once.

        :Example:

        >>> from pysnptools.eigenreader import EigenNpz
        >>> from pysnptools.util import example_file # Download and return local file name
        >>> npz_file = example_file("pysnptools/examples/toydata.eigen.npz")
        >>> eigen_on_disk = EigenNpz(npz_file) # Specify SNP data on disk
        >>> print(eigen_on_disk.eid_to_index(['eid_2','eid_9'])) #Find the indexes for two eids.
        [2 9]
        """
        return self.col_to_index(list)

    @property
    def val_shape(self):
        '''
        Tells the shape of value for a given individual and SNP. For EigenReaders always returns 3.
        '''
        return 2


    def __getitem__(self, iid_indexer_and_snp_indexer):
        from pysnptools.eigenreader._subset import _EigenSubset
        iid_indexer, snp_indexer = iid_indexer_and_snp_indexer
        return _EigenSubset(self, iid_indexer, snp_indexer)

    @staticmethod
    def _as_eigendata(eigenreader, force_python_only, order, dtype, num_threads):
        '''
        Like 'read' except won't read if already a EigenData
        '''
        from pysnptools.eigenreader import EigenData #must import here to avoid cycles
        dtype = np.dtype(dtype)

        if hasattr(eigenreader,'vectors') and eigenreader.vectors.dtype==dtype and (order=="A" or (order=="C" and eigenreader.vectors.flags["C_CONTIGUOUS"]) or (order=="F" and eigenreader.vectors.flags["F_CONTIGUOUS"])):
            return eigenreader
        else:
            return eigenreader.read(order=order,dtype=dtype,view_ok=True, num_threads=num_threads)
    
    def copyinputs(self, copier):
        raise NotImplementedError

    def _assert_iid_eid_values(self,check_vectors):
        if check_vectors:
            assert len(self._val.shape)==2, "vectors should have 2 dimensions"
            assert self._val.shape == (len(self._row),len(self._col)), "vectors shape should match that of iid_count x eid_count"
        assert self._row.dtype.type is np.str_ and len(self._row.shape)==2 and self._row.shape[1]==2, "iid should be dtype str, have two dimensions, and the second dimension should be size 2"
        assert self._col.dtype.type is np.str_ and len(self._col.shape)==1, "eid should be of dtype of str and one dimensional"


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    if False: # cmk
        from pysnptools.eigenreader import cmk
        eigen_on_disk = EigenNpz('../examples/toydata.eigen.npz')
        print(eigen_on_disk.pos[:4,].astype('int')) # print position information for the first three eids: #The '...' is for possible space char

    import doctest
    doctest.testmod(optionflags=doctest.ELLIPSIS|doctest.NORMALIZE_WHITESPACE)
    # There is also a unit test case in 'pysnptools\test.py' that calls this doc t
    print("done")