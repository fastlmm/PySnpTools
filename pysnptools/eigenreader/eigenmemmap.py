import logging
import os
import shutil
import numpy as np
import unittest
import doctest
import pysnptools.util as pstutil
from pysnptools.pstreader import PstData
from pysnptools.pstreader import PstMemMap
from pysnptools.eigenreader import EigenReader, EigenData
from pysnptools.util import log_in_place

#!!!cmk update docs
class EigenMemMap(PstMemMap,EigenData):
    '''
    A :class:`.EigenData` that keeps its data in a memory-mapped file. This allows data large than fits in main memory.

    See :class:`.EigenData` for general examples of using EigenData.

    **Constructor:**
        :Parameters: **filename** (*string*) -- The *\*.eigen.memmap* file to read.
        
        Also see :meth:`.EigenMemMap.empty` and :meth:`.EigenMemMap.write`.

        :Example:

        >>> from pysnptools.eigenreader import EigenMemMap
        >>> from pysnptools.util import example_file # Download and return local file name
        >>> mem_map_file = example_file("pysnptools/examples/tiny.eigen.memmap")
        >>> eigen_mem_map = EigenMemMap(mem_map_file)
        >>> print(eigen_mem_map.val[0,1], eigen_mem_map.iid_count, eigen_mem_map.sid_count)
        [0.43403135 0.28289911 0.28306954] 25 10

    **Methods inherited from** :class:`.EigenData`

        :meth:`.EigenData.allclose`

    **Methods beyond** :class:`.EigenReader`

    '''

    def __init__(self, *args, **kwargs):
        super(EigenMemMap, self).__init__(*args, **kwargs)

    @property
    def vectors(self):
        """The 2D NumPy array of floats that represents the column eigen vectors. You can get or set this property.

        >>> from pysnptools.eigenreader import EigenData # !!!cmk later change this example to read from disk
        >>> eigendata = EigenData(row=[['fam0','iid0'],['fam0','iid1']],
        ...                     values=[0.88196601, 3.11803399],
        ...                     vectors=[[ 0.22975292, -0.97324899],
        ...                              [-0.97324899, -0.22975292]])
        >>> print(eigendata.vectors[:,0]) # print the first eigenvector
        [ 0.22975292 -0.97324899]
        """
        return self._val


    @vectors.setter
    def vectors(self, new_value):
        self._run_once()
        if self._val is new_value:
            return
        raise Exception("EigenMemMap vectors's cannot be set to a different array")


    @property
    def offset(self):
        '''The byte position in the file where the memory-mapped values start.
       
        (The disk space before this is used to store :attr:`EigenReader.iid`, etc. information.
        This property is useful when interfacing with, for example, external Fortran and C matrix libraries.)
        
        '''
        self._run_once()
        return self._offset

    @property
    def filename(self):
        '''The name of the memory-mapped file
        '''
        #Don't need '_run_once'
        return self._filename

    @staticmethod
    def empty(iid, values, filename, eid=None, order="F",dtype=np.float64):
        '''Create an empty :class:`.EigenMemMap` on disk.

        :param iid: The :attr:`EigenReader.iid` information
        :type iid: an array of string pairs

        :param sid: The :attr:`EigenReader.sid` information
        :type sid: an array of strings

        :param filename: name of memory-mapped file to create
        :type filename: string

        :param pos: optional -- The additional :attr:`EigenReader.pos` information associated with each sid. Default: None
        :type pos: an array of numeric triples

        :param order: {'F' (default), 'C'}, optional -- Specify the order of the ndarray.
        :type order: string or None

        :param dtype: {numpy.float64 (default), numpy.float32}, optional -- The data-type for the :attr:`EigenMemMap.val` ndarray.
        :type dtype: data-type

        :rtype: :class:`.EigenMemMap`

        >>> import pysnptools.util as pstutil
        >>> from pysnptools.eigenreader import EigenMemMap
        >>> filename = "tempdir/tiny.eigen.memmap"
        >>> pstutil.create_directory_if_necessary(filename)
        >>> eigen_mem_map = EigenMemMap.empty(iid=[['fam0','iid0'],['fam0','iid1']], sid=['snp334','snp349','snp921'],filename=filename,order="F",dtype=np.float64)
        >>> eigen_mem_map.val[:,:,:] = [[[.5,.5,0],[0,0,1],[.5,.5,0]],
        ...                            [[0,1.,0],[0,.75,.25],[.5,.5,0]]]
        >>> eigen_mem_map.flush()

        '''
        #!!!cmk why is this code so different from EigenData's similar code?
        if eid is None:
            eid = [["", f"eid{eid_index}"] for eid_index in range(len(values))]
        self = EigenMemMap(filename)
        self._empty_inner(row=iid, col=eid, filename=filename, row_property=None, col_property=values,order=order,dtype=dtype,val_shape=None)
        return self

    def flush(self):
        '''Flush :attr:`EigenMemMap.val` to disk and close the file. (If values or properties are accessed again, the file will be reopened.)

        >>> import pysnptools.util as pstutil
        >>> from pysnptools.eigenreader import EigenMemMap
        >>> filename = "tempdir/tiny.eigen.memmap"
        >>> pstutil.create_directory_if_necessary(filename)
        >>> eigen_mem_map = EigenMemMap.empty(iid=[['fam0','iid0'],['fam0','iid1']], sid=['snp334','snp349','snp921'],filename=filename,order="F",dtype=np.float64)
        >>> eigen_mem_map.val[:,:,:] = [[[.5,.5,0],[0,0,1],[.5,.5,0]],
        ...                            [[0,1.,0],[0,.75,.25],[.5,.5,0]]]
        >>> eigen_mem_map.flush()

        '''
        if self._ran_once:
            self.val.flush()
            del self._val
            self._ran_once = False


    @staticmethod
    def write(filename, eigenreader, order='A', dtype=None, block_size=None, num_threads=None):
        """Writes a :class:`EigenReader` to :class:`EigenMemMap` format.

        :param filename: the name of the file to create
        :type filename: string
        :param eigenreader: The data that should be written to disk. It can also be any eigenreader, for example, :class:`.EigenNpz`, :class:`.EigenData`, or
           another :class:`.Bgen`.
        :type eigenreader: :class:`EigenReader`
        :param order: {'A' (default), 'F', 'C'}, optional -- Specify the order of the ndarray. By default, will match the order of the input if knowable; otherwise, 'F'
        :type order: string or None
        :param dtype: {None (default), numpy.float64, numpy.float32}, optional -- The data-type for the :attr:`EigenMemMap.val` ndarray.
             By default, will match the order of the input if knowable; otherwise np.float64.
        :type dtype: data-type
        :param block_size: The number of SNPs to read in a batch from *eigenreader*. Defaults to a *block_size* such that *block_size* \* *iid_count* is about 100,000.
        :type block_size: number
        :param num_threads: optional -- The number of threads with which to write data. Defaults to all available
            processors. Can also be set with these environment variables (listed in priority order):
            'PST_NUM_THREADS', 'NUM_THREADS', 'MKL_NUM_THREADS'.
        :type num_threads: None or int
        :rtype: :class:`.EigenMemMap`

        >>> import pysnptools.util as pstutil
        >>> from pysnptools.eigenreader import Bgen, EigenMemMap
        >>> from pysnptools.util import example_file # Download and return local file name
        >>> bgen_file = example_file("pysnptools/examples/2500x100.bgen")
        >>> eigenreader = Bgen(bgen_file)[:,:10] #Create a reader for the first 10 SNPs
        >>> pstutil.create_directory_if_necessary("tempdir/tiny.eigen.memmap")
        >>> EigenMemMap.write("tempdir/tiny.eigen.memmap",eigenreader)      # Write eigenreader in EigenMemMap format
        EigenMemMap('tempdir/tiny.eigen.memmap')

        """
        block_size = block_size or max((100_000)//max(1,eigenreader.row_count),1)

        if hasattr(eigenreader,'val'):
            order = PstMemMap._order(eigenreader) if order=='A' else order
            dtype = dtype or eigenreader.val.dtype
        else:
            order = 'F' if order=='A' else order
            dtype = dtype or np.float64
        dtype = np.dtype(dtype)

        self = PstMemMap.empty(eigenreader.iid, eigenreader.col, filename+'.temp', row_property=eigenreader.row_property, col_property=eigenreader.col_property,order=order,dtype=dtype, val_shape=None)
        if hasattr(eigenreader,'val'):
            self.val[...] = eigenreader.val
        else:
            start = 0
            with log_in_place("EigenMemMap writing sid_index ", logging.INFO) as updater:
                while start < eigenreader.sid_count:
                    updater('{0} of {1}'.format(start,eigenreader.sid_count))
                    eigendata = eigenreader[:,start:start+block_size].read(order=order,dtype=dtype,num_threads=num_threads)
                    self.val[:,start:start+eigendata.sid_count,:] = eigendata.val
                    start += eigendata.sid_count

        self.flush()
        if os.path.exists(filename):
           os.remove(filename) 
        shutil.move(filename+'.temp',filename)
        logging.debug("Done writing " + filename)
        return EigenMemMap(filename)



    def _run_once(self):
            if (self._ran_once):
                return
            iid_ascii,eid_ascii,vectors,row_property,values = self._run_once_inner()
            iid = np.array(iid_ascii,dtype='str') #!!!avoid this copy when not needed
            eid = np.array(eid_ascii,dtype='str') #!!!avoid this copy when not needed

            EigenData.__init__(self,iid=iid,eid=eid,vectors=vectors,values=values,name="np.memmap('{0}')".format(self._filename))

class TestEigenMemMap(unittest.TestCase):     

    def test1(self):        
        old_dir = os.getcwd()
        os.chdir(os.path.dirname(os.path.realpath(__file__)))

        filename2 = "tempdir/tiny.eigen.memmap"
        pstutil.create_directory_if_necessary(filename2)
        eigenreader2 = EigenMemMap.empty(iid=[['fam0','iid0'],['fam0','iid1']], sid=['snp334','snp349','snp921'],filename=filename2,order="F",dtype=np.float64)
        assert isinstance(eigenreader2.val,np.memmap)
        eigenreader2.val[:,:,:] = [[[.5,.5,0],[0,0,1],[.5,.5,0]],[[0,1.,0],[0,.75,.25],[.5,.5,0]]]
        assert np.array_equal(eigenreader2[[1],[1]].read(view_ok=True).val,np.array([[[0,.75,.25]]]))
        eigenreader2.flush()
        assert isinstance(eigenreader2.val,np.memmap)
        assert np.array_equal(eigenreader2[[1],[1]].read(view_ok=True).val,np.array([[[0,.75,.25]]]))
        eigenreader2.flush()

        eigenreader3 = EigenMemMap(filename2)
        assert np.array_equal(eigenreader3[[1],[1]].read(view_ok=True).val,np.array([[[0,.75,.25]]]))
        assert isinstance(eigenreader3.val,np.memmap)

        logging.info("in TestEigenMemMap test1")
        eigenreader = EigenMemMap('../examples/tiny.eigen.memmap')
        assert eigenreader.iid_count == 25
        assert eigenreader.sid_count == 10
        assert isinstance(eigenreader.val,np.memmap)

        eigendata = eigenreader.read(view_ok=True)
        assert isinstance(eigendata.val,np.memmap)
        os.chdir(old_dir)

    def test2(self):
        from pysnptools.eigenreader import Bgen

        old_dir = os.getcwd()
        os.chdir(os.path.dirname(os.path.realpath(__file__)))

        bgen = Bgen('../examples/example.bgen')
        eigenmemmap = EigenMemMap.write("tempdir/bgentomemmap.eigen.memamp",bgen)
        assert EigenData.allclose(bgen.read(),eigenmemmap.read(),equal_nan=True)
        os.chdir(old_dir)

    def test_doctest(self):
        import pysnptools.eigenreader.eigenmemmap as mod_mm
        import doctest

        old_dir = os.getcwd()
        os.chdir(os.path.dirname(os.path.realpath(__file__)))
        result = doctest.testmod(mod_mm)
        os.chdir(old_dir)
        assert result.failed == 0, "failed doc test: " + __file__

def getTestSuite():
    """
    set up composite test suite
    """
    
    test_suite = unittest.TestSuite([])
    test_suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestEigenMemMap))
    return test_suite


if __name__ == "__main__":
    logging.basicConfig(level=logging.WARN)

    suites = getTestSuite()
    r = unittest.TextTestRunner(failfast=True)
    ret = r.run(suites)
    assert ret.wasSuccessful()


    result = doctest.testmod(optionflags=doctest.ELLIPSIS|doctest.NORMALIZE_WHITESPACE)
    assert result.failed == 0, "failed doc test: " + __file__
