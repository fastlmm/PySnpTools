import numpy as np
import os
import os.path
from itertools import *
import pandas as pd
import logging
import time
import pysnptools.util as pstutil
from pysnptools.pstreader import PstReader
from pysnptools.snpreader import SnpData
import warnings
import pysnptools.standardizer as stdizer
from pysnptools.snpreader._dist2snp import _Dist2Snp

#!!why do the examples use ../tests/datasets instead of "examples"?
class DistReader(PstReader):
    """A DistReader is one of three things:

    * A class such as :class:`.Bgen` for you to specify a file with data. For example,

        >>> from pysnptools.distreader import Bgen
        >>> from pysnptools.util import example_file # Download and return local file name
        >>> bgen_file = example_file("pysnptools/examples/2500x100.bgen")
        >>> dist_on_disk = Bgen(bgen_file)
        >>> print(dist_on_disk) # prints the name of the file reader
        Bgen('...pysnptools/examples/2500x100.bgen')
        >>> dist_on_disk.sid_count # prints the number of SNPS (but doesn't read any SNP distribution values)
        100

    * A :class:`.DistData` class that holds SNP distribution data in memory, typically after reading it from disk:

        >>> from pysnptools.util import example_file # Download and return local file name
        >>> bgen_file = example_file("pysnptools/examples/2500x100.bgen")
        >>> dist_on_disk = Bgen(bgen_file)
        >>> distdata1 = dist_on_disk.read() #reads the SNP distribution values
        >>> print(type(distdata1.val)) # The val property is an 3-D ndarray of SNP distribution values
        <class 'numpy.ndarray'>
        >>> print(distdata1) # prints the name in-memory SNP distribution reader.
        DistData(Bgen('...pysnptools/examples/2500x100.bgen'))
        >>> distdata1.iid_count #prints the number of iids (number of individuals) in this in-memory data
        2500

    * A subset of any DistReader, specified with "[ *iid_index* , *sid_index* ]", to read only some SNP distribution values. It can
      also be used to re-order the values.

        >>> from pysnptools.util import example_file # Download and return local file name
        >>> bgen_file = example_file("pysnptools/examples/2500x100.bgen")
        >>> dist_on_disk = Bgen(bgen_file)
        >>> subset_on_disk = dist_on_disk[[3,4],::2] # specification for a subset of the data on disk. No SNP distriubtion values are read yet.
        >>> print(subset_on_disk.sid_count) # prints the number of sids in this subset (but still doesn't read any SNP distribution values)
        50
        >>> print(subset_on_disk) #prints a specification of 'subset_on_disk'
        Bgen('...pysnptools/examples/2500x100.bgen')[[3,4],::2]
        >>> distdata_subset = subset_on_disk.read() # efficiently reads the specified subset of values from the disk
        >>> print(distdata_subset) # prints the specification of the in-memory SNP distribution information
        DistData(Bgen('...pysnptools/examples/2500x100.bgen')[[3,4],::2])
        >>> print((int(distdata_subset.val.shape[0]), int(distdata_subset.val.shape[1]))) # The dimensions of the ndarray of SNP distriubtion values
        (2, 50)

    The DistReaders Classes

        ========================= =================== ====================== ================== ======================
        *Class*                   *Format*            *Random Access*        *Suffixes*         *Write* method?
        :class:`.DistData`        in-memory floats    Yes                    *n/a*              *n/a*              
        :class:`.Bgen`            binary, floats      Yes (by sid)           \*.bgen            Yes              
        :class:`.DistNpz`         binary, floats      No                     .dist.npz          Yes
        :class:`.DistHdf5`        binary, floats      Yes (by sid or iid)    .dist.hdf5         Yes
        :class:`.DistMemMap`      mem-mapped floats   Yes                    .dist.memmap       Yes              
        ========================= =================== ====================== ================== ======================
    
  
    Methods & Properties:

        Every DistReader, such as :class:`.Bgen` and :class:`.DistData`, has these properties: :attr:`iid`, :attr:`iid_count`, :attr:`sid`, :attr:`sid_count`,
        :attr:`pos` and these methods: :meth:`read`, :meth:`iid_to_index`, :meth:`sid_to_index`, :meth:`as_snp`. See below for details.

        :class:`.DistData` is a DistReader so it supports the above properties and methods. In addition, it supports property :attr:`DistData.val`.

        See below for details.

        Many of the classes, such as :class:`.Bgen`, also provide a static :meth:`Bgen.write` method for writing :class:`.DistReader` or :class:`.DistData` to disk.

        >>> from pysnptools.distreader import DistHdf5, Bgen
        >>> import pysnptools.util as pstutil
        >>> hdf5_file = example_file("pysnptools/examples/toydata.snpmajor.dist.hdf5")
        >>> distreader = DistHdf5(hdf5_file)[:,:10] # A reader for the first 10 SNPs in Hdf5 format
        >>> pstutil.create_directory_if_necessary("tempdir/toydata10.bgen")
        >>> Bgen.write("tempdir/toydata10.bgen",distreader)        # Write data in BGEN format
        Bgen('tempdir/toydata10.bgen')


    iids and sids:

        *same as with* :class:`.SnpReader`

    Selecting and Reordering Individuals and SNPs

        *same as with* :class:`.SnpReader`
        
    When Data is Read:

        *same as with* :class:`.SnpReader`

    When Data is Re-Read and Copied:

        *same as with* :class:`.SnpReader`

    Avoiding Unwanted ndarray Allocations

        *same as with* :class:`.SnpReader`

    The :meth:`read` Method
  
        By default the :meth:`read` returns a 3-D ndarray of numpy.float64 laid out in memory in F-contiguous order
        (iid-index varies the fastest). You may, instead,
        ask for numpy.float32 or for C-contiguous order or any order. See :meth:`read` for details.


    Details of Methods & Properties:
    """

    def __init__(self, *args, **kwargs):
        super(DistReader, self).__init__(*args, **kwargs)

    @property
    def iid(self):
        """A ndarray of the iids. Each iid is a ndarray of two strings (a family ID and a individual ID) that identifies an individual.

        :rtype: ndarray of strings with shape [:attr:`.iid_count`,2]

        This property (to the degree practical) reads only iid and sid data from the disk, not SNP distribution data. Moreover, the iid and sid data is read from file only once.

        :Example:

        >>> from pysnptools.distreader import Bgen
        >>> from pysnptools.util import example_file # Download and return local file name
        >>> bgen_file = example_file("pysnptools/examples/2500x100.bgen")
        >>> dist_on_disk = Bgen(bgen_file)
        >>> print(dist_on_disk.iid[:3]) # print the first three iids
        [['0' 'iid_0']
         ['0' 'iid_1']
         ['0' 'iid_2']]
        """
        return self.row

    @property
    def iid_count(self):
        """number of iids

        :rtype: integer

        This property (to the degree practical) reads only iid and sid data from the disk, not SNP distribution data. Moreover, the iid and sid data is read from file only once.
        """
        return self.row_count

    @property
    def sid(self):
        """A ndarray of the sids. Each sid is a string that identifies a SNP.

        :rtype: ndarray (length :attr:`.sid_count`) of strings

        This property (to the degree practical) reads only iid and sid data from the disk, not SNP distribution data. Moreover, the iid and sid data is read from file only once.

        :Example:

        >>> from pysnptools.distreader import Bgen
        >>> from pysnptools.util import example_file # Download and return local file name
        >>> bgen_file = example_file("pysnptools/examples/2500x100.bgen")
        >>> dist_on_disk = Bgen(bgen_file)
        >>> print(dist_on_disk.sid[:9]) # print the first nine sids
        ['sid_0' 'sid_1' 'sid_2' 'sid_3' 'sid_4' 'sid_5' 'sid_6' 'sid_7' 'sid_8']
        """
        return self.col

    @property
    def sid_count(self):
        """number of sids

        :rtype: integer

        This property (to the degree practical) reads only iid and sid data from the disk, not SNP distribution data. Moreover, the iid and sid data is read from file only once.

        """
        return self.col_count

    #!!document that chr must not be X,Y,M only numbers (as per the PLINK DistNpz format)
    #!!Also what about telling the ref and alt allele? Also, what about tri and quad alleles, etc?
    @property
    def pos(self):
        """A ndarray of the position information for each sid. Each element is a ndarray of three numpy.numbers (chromosome, genetic distance, basepair distance).

        :rtype: ndarray of float64 with shape [:attr:`.sid_count`, 3]

        This property (to the degree practical) reads only iid and sid data from the disk, not SNP distribution data. Moreover, the iid and sid data is read from file only once.

        :Example:

        >>> from pysnptools.distreader import Bgen
        >>> from pysnptools.util import example_file # Download and return local file name
        >>> bgen_file = example_file("pysnptools/examples/2500x100.bgen")
        >>> dist_on_disk = Bgen(bgen_file)
        >>> print(dist_on_disk.pos[:4,].astype('int')) # print position information for the first four sids: #The '...' is for possible space char
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


    def _read(self, iid_index_or_none, sid_index_or_none, order, dtype, force_python_only, view_ok, num_threads):
        raise NotImplementedError
    
    #!!check that views always return contiguous memory by default
    def read(self, order='F', dtype=np.float64, force_python_only=False, view_ok=False, num_threads=None):
        """Reads the SNP values and returns a :class:`.DistData` (with :attr:`DistData.val` property containing a new 3D ndarray of the SNP distribution values).

        :param order: {'F' (default), 'C', 'A'}, optional -- Specify the order of the ndarray. If order is 'F' (default),
            then the array will be in F-contiguous order (iid-index varies the fastest).
            If order is 'C', then the returned array will be in C-contiguous order (sid-index varies the fastest).
            If order is 'A', then the :attr:`DistData.val`
            ndarray may be in any order (either C-, Fortran-contiguous).
        :type order: string or None

        :param dtype: {numpy.float64 (default), numpy.float32}, optional -- The data-type for the :attr:`DistData.val` ndarray.
        :type dtype: data-type

        :param force_python_only: optional -- If False (default), may use outside library code. If True, requests that the read
            be done without outside library code.
        :type force_python_only: bool

        :param view_ok: optional -- If False (default), allocates new memory for the :attr:`DistData.val`'s ndarray. If True,
            if practical and reading from a :class:`DistData`, will return a new 
            :class:`DistData` with a ndarray shares memory with the original :class:`DistData`.
            Typically, you'll also wish to use "order='A'" to increase the chance that sharing will be possible.
            Use these parameters with care because any change to either ndarraywill effect
            the others. Also keep in mind that :meth:`read` relies on ndarray's mechanisms to decide whether to actually
            share memory and so it may ignore your suggestion and allocate a new ndarray anyway.
        :type view_ok: bool

        :param num_threads: optional -- The number of threads with which to read data. Defaults to all available
            processors. Can also be set with these environment variables (listed in priority order):
            'PST_NUM_THREADS', 'NUM_THREADS', 'MKL_NUM_THREADS'.
        :type num_threads: None or int

        :rtype: :class:`.DistData`

        Calling the method again causes the SNP distribution values to be re-read and creates a new in-memory :class:`.DistData` with a new ndarray of SNP values.

        If you request the values for only a subset of the sids or iids, (to the degree practical) only that subset will be read from disk.

        :Example:

        >>> from pysnptools.distreader import Bgen
        >>> from pysnptools.util import example_file # Download and return local file name
        >>> bgen_file = example_file("pysnptools/examples/2500x100.bgen")
        >>> dist_on_disk = Bgen(bgen_file) # Specify SNP data on disk
        >>> distdata1 = dist_on_disk.read() # Read all the SNP data returning a DistData instance
        >>> print(type(distdata1.val).__name__) # The DistData instance contains a ndarray of the data.
        ndarray
        >>> subset_distdata = dist_on_disk[:,::2].read() # From the disk, read SNP values for every other sid
        >>> print(subset_distdata.val[0,0]) # Print the first SNP value in the subset
        [0.466804   0.38812848 0.14506752]
        >>> subsub_distdata = subset_distdata[:10,:].read(order='A',view_ok=True) # Create an in-memory subset of the subset with SNP values for the first ten iids. Share memory if practical.
        >>> import numpy as np
        >>> # print np.may_share_memory(subset_distdata.val, subsub_distdata.val) # Do the two ndarray's share memory? They could. Currently they won't.       
        """
        dtype = np.dtype(dtype)
        val = self._read(None, None, order, dtype, force_python_only, view_ok, num_threads)
        from pysnptools.distreader import DistData
        ret = DistData(self.iid,self.sid,val,pos=self.pos,name=str(self))
        return ret

    def as_snp(self, max_weight=2.0, block_size=None):
        """Returns a :class:`pysnptools.snpreader.SnpReader` such that turns the probability distribution into an expected value.

        For example, if the probability distribution is [0.466804   0.38812848 0.14506752] and the max_weight is 2, then the expected
        value will be (0.466804*0+0.38812848*.5+0.14506752*1)*2, 0.67826352 .

        :param max_weight: optional -- The expected values will range from 0 to *max_weight*. Default of 2.
        :type max_weight: number

        :param block_size: optional -- Default of None (meaning to load all). Suggested number of sids to read into memory at a time.
        :type block_size: int or None

        :rtype: class:`SnpReader`

        :Example:

        >>> from pysnptools.distreader import Bgen
        >>> from pysnptools.util import example_file # Download and return local file name
        >>> bgen_file = example_file("pysnptools/examples/2500x100.bgen")
        >>> distreader = Bgen(bgen_file) # Specify distribution data on disk
        >>> print(distreader[0,0].read().val)
        [[[0.466804   0.38812848 0.14506752]]]
        >>> snpreader = distreader.as_snp(max_weight=2)
        >>> print(snpreader[0,0].read().val)
        [[0.67826352]]
        """
        dist2snp = _Dist2Snp(self,max_weight=max_weight,block_size=block_size)
        return dist2snp

    def iid_to_index(self, list):
        """Takes a list of iids and returns a list of index numbers

        :param list: list of iids
        :type order: list of list of strings

        :rtype: ndarray of int
        
        This method (to the degree practical) reads only iid and sid data from the disk, not SNP value data. Moreover, the iid and sid data is read from file only once.

        :Example:

        >>> from pysnptools.distreader import Bgen
        >>> from pysnptools.util import example_file # Download and return local file name
        >>> bgen_file = example_file("pysnptools/examples/2500x100.bgen")
        >>> dist_on_disk = Bgen(bgen_file) # Specify SNP data on disk
        >>> print(dist_on_disk.iid_to_index([['0','iid_2'],['0','iid_1']])) #Find the indexes for two iids.
        [2 1]
        """
        return self.row_to_index(list)

    def sid_to_index(self, list):
        """Takes a list of sids and returns a list of index numbers

        :param list: list of sids
        :type list: list of strings

        :rtype: ndarray of int
        
        This method (to the degree practical) reads only iid and sid data from the disk, not SNP value data. Moreover, the iid and sid data is read from file only once.

        :Example:

        >>> from pysnptools.distreader import Bgen
        >>> from pysnptools.util import example_file # Download and return local file name
        >>> bgen_file = example_file("pysnptools/examples/2500x100.bgen")
        >>> dist_on_disk = Bgen(bgen_file) # Specify SNP data on disk
        >>> print(dist_on_disk.sid_to_index(['sid_2','sid_9'])) #Find the indexes for two sids.
        [2 9]
        """
        return self.col_to_index(list)

    @property
    def val_shape(self):
        '''
        Tells the shape of value for a given individual and SNP. For DistReaders always returns 3.
        '''
        return 3


    def __getitem__(self, iid_indexer_and_snp_indexer):
        from pysnptools.distreader._subset import _DistSubset
        iid_indexer, snp_indexer = iid_indexer_and_snp_indexer
        return _DistSubset(self, iid_indexer, snp_indexer)

    @staticmethod
    def _as_distdata(distreader, force_python_only, order, dtype, num_threads):
        '''
        Like 'read' except won't read if already a DistData
        '''
        from pysnptools.distreader import DistData #must import here to avoid cycles
        dtype = np.dtype(dtype)

        if hasattr(distreader,'val') and distreader.val.dtype==dtype and (order=="A" or (order=="C" and distreader.val.flags["C_CONTIGUOUS"]) or (order=="F" and distreader.val.flags["F_CONTIGUOUS"])):
            return distreader
        else:
            return distreader.read(order=order,dtype=dtype,view_ok=True, num_threads=num_threads)
    
    def copyinputs(self, copier):
        raise NotImplementedError

    def _assert_iid_sid_pos(self,check_val):
        if check_val:
            assert len(self._val.shape)==3 and self._val.shape[-1]==3, "val should have 3 dimensions and the last dimension should have size 3"
            assert self._val.shape == (len(self._row),len(self._col),3), "val shape should match that of iid_count x sid_count"
        assert self._row.dtype.type is np.str_ and len(self._row.shape)==2 and self._row.shape[1]==2, "iid should be dtype str, have two dimensions, and the second dimension should be size 2"
        assert self._col.dtype.type is np.str_ and len(self._col.shape)==1, "sid should be of dtype of str and one dimensional"




if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    if False:
        from pysnptools.distreader import Bgen
        dist_on_disk = Bgen('../examples/2500x100.bgen')
        print(dist_on_disk.pos[:4,].astype('int')) # print position information for the first three sids: #The '...' is for possible space char

    import doctest
    doctest.testmod(optionflags=doctest.ELLIPSIS|doctest.NORMALIZE_WHITESPACE)
    # There is also a unit test case in 'pysnptools\test.py' that calls this doc t
    print("done")