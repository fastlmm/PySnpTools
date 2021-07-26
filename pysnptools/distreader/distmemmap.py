import logging
import os
import shutil
import numpy as np
import unittest
import doctest
import pysnptools.util as pstutil
from pysnptools.pstreader import PstData
from pysnptools.pstreader import PstMemMap
from pysnptools.distreader import DistReader, DistData
from pysnptools.util import log_in_place


class DistMemMap(PstMemMap,DistData):
    '''
    A :class:`.DistData` that keeps its data in a memory-mapped file. This allows data large than fits in main memory.

    See :class:`.DistData` for general examples of using DistData.

    **Constructor:**
        :Parameters: **filename** (*string*) -- The *\*.dist.memmap* file to read.
        
        Also see :meth:`.DistMemMap.empty` and :meth:`.DistMemMap.write`.

        :Example:

        >>> from pysnptools.distreader import DistMemMap
        >>> from pysnptools.util import example_file # Download and return local file name
        >>> mem_map_file = example_file("pysnptools/examples/tiny.dist.memmap")
        >>> dist_mem_map = DistMemMap(mem_map_file)
        >>> print(dist_mem_map.val[0,1], dist_mem_map.iid_count, dist_mem_map.sid_count)
        [0.43403135 0.28289911 0.28306954] 25 10

    **Methods inherited from** :class:`.DistData`

        :meth:`.DistData.allclose`

    **Methods beyond** :class:`.DistReader`

    '''

    def __init__(self, *args, **kwargs):
        super(DistMemMap, self).__init__(*args, **kwargs)

    @property
    def val(self):
        """The 3D NumPy memmap array of floats that represents the distribution of SNP values. You can get this property, but cannot set it (except with itself)


        >>> from pysnptools.distreader import DistMemMap
        >>> from pysnptools.util import example_file # Download and return local file name
        >>> mem_map_file = example_file("pysnptools/examples/tiny.dist.memmap")
        >>> dist_mem_map = DistMemMap(mem_map_file)
        >>> print(dist_mem_map.val[0,1])
        [0.43403135 0.28289911 0.28306954]
        """
        self._run_once()
        return self._val


    @val.setter
    def val(self, new_value):
        self._run_once()
        if self._val is new_value:
            return
        raise Exception("DistMemMap val's cannot be set to a different array")


    @property
    def offset(self):
        '''The byte position in the file where the memory-mapped values start.
       
        (The disk space before this is used to store :attr:`DistReader.iid`, etc. information.
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
    def empty(iid, sid, filename, pos=None,order="F",dtype=np.float64):
        '''Create an empty :class:`.DistMemMap` on disk.

        :param iid: The :attr:`DistReader.iid` information
        :type iid: an array of string pairs

        :param sid: The :attr:`DistReader.sid` information
        :type sid: an array of strings

        :param filename: name of memory-mapped file to create
        :type filename: string

        :param pos: optional -- The additional :attr:`DistReader.pos` information associated with each sid. Default: None
        :type pos: an array of numeric triples

        :param order: {'F' (default), 'C'}, optional -- Specify the order of the ndarray.
        :type order: string or None

        :param dtype: {numpy.float64 (default), numpy.float32}, optional -- The data-type for the :attr:`DistMemMap.val` ndarray.
        :type dtype: data-type

        :rtype: :class:`.DistMemMap`

        >>> import pysnptools.util as pstutil
        >>> from pysnptools.distreader import DistMemMap
        >>> filename = "tempdir/tiny.dist.memmap"
        >>> pstutil.create_directory_if_necessary(filename)
        >>> dist_mem_map = DistMemMap.empty(iid=[['fam0','iid0'],['fam0','iid1']], sid=['snp334','snp349','snp921'],filename=filename,order="F",dtype=np.float64)
        >>> dist_mem_map.val[:,:,:] = [[[.5,.5,0],[0,0,1],[.5,.5,0]],
        ...                            [[0,1.,0],[0,.75,.25],[.5,.5,0]]]
        >>> dist_mem_map.flush()

        '''

        self = DistMemMap(filename)
        self._empty_inner(row=iid, col=sid, filename=filename, row_property=None, col_property=pos,order=order,dtype=dtype,val_shape=3)
        return self

    def flush(self):
        '''Flush :attr:`DistMemMap.val` to disk and close the file. (If values or properties are accessed again, the file will be reopened.)

        >>> import pysnptools.util as pstutil
        >>> from pysnptools.distreader import DistMemMap
        >>> filename = "tempdir/tiny.dist.memmap"
        >>> pstutil.create_directory_if_necessary(filename)
        >>> dist_mem_map = DistMemMap.empty(iid=[['fam0','iid0'],['fam0','iid1']], sid=['snp334','snp349','snp921'],filename=filename,order="F",dtype=np.float64)
        >>> dist_mem_map.val[:,:,:] = [[[.5,.5,0],[0,0,1],[.5,.5,0]],
        ...                            [[0,1.,0],[0,.75,.25],[.5,.5,0]]]
        >>> dist_mem_map.flush()

        '''
        if self._ran_once:
            self.val.flush()
            del self._val
            self._ran_once = False


    @staticmethod
    def write(filename, distreader, order='A', dtype=None, block_size=None, num_threads=None):
        """Writes a :class:`DistReader` to :class:`DistMemMap` format.

        :param filename: the name of the file to create
        :type filename: string
        :param distreader: The data that should be written to disk. It can also be any distreader, for example, :class:`.DistNpz`, :class:`.DistData`, or
           another :class:`.Bgen`.
        :type distreader: :class:`DistReader`
        :param order: {'A' (default), 'F', 'C'}, optional -- Specify the order of the ndarray. By default, will match the order of the input if knowable; otherwise, 'F'
        :type order: string or None
        :param dtype: {None (default), numpy.float64, numpy.float32}, optional -- The data-type for the :attr:`DistMemMap.val` ndarray.
             By default, will match the order of the input if knowable; otherwise np.float64.
        :type dtype: data-type
        :param block_size: The number of SNPs to read in a batch from *distreader*. Defaults to a *block_size* such that *block_size* \* *iid_count* is about 100,000.
        :type block_size: number
        :param num_threads: optional -- The number of threads with which to write data. Defaults to all available
            processors. Can also be set with these environment variables (listed in priority order):
            'PST_NUM_THREADS', 'NUM_THREADS', 'MKL_NUM_THREADS'.
        :type num_threads: None or int
        :rtype: :class:`.DistMemMap`

        >>> import pysnptools.util as pstutil
        >>> from pysnptools.distreader import Bgen, DistMemMap
        >>> from pysnptools.util import example_file # Download and return local file name
        >>> bgen_file = example_file("pysnptools/examples/2500x100.bgen")
        >>> distreader = Bgen(bgen_file)[:,:10] #Create a reader for the first 10 SNPs
        >>> pstutil.create_directory_if_necessary("tempdir/tiny.dist.memmap")
        >>> DistMemMap.write("tempdir/tiny.dist.memmap",distreader)      # Write distreader in DistMemMap format
        DistMemMap('tempdir/tiny.dist.memmap')

        """
        block_size = block_size or max((100_000)//max(1,distreader.row_count),1)

        if hasattr(distreader,'val'):
            order = PstMemMap._order(distreader) if order=='A' else order
            dtype = dtype or distreader.val.dtype
        else:
            order = 'F' if order=='A' else order
            dtype = dtype or np.float64
        dtype = np.dtype(dtype)

        self = PstMemMap.empty(distreader.row, distreader.col, filename+'.temp', row_property=distreader.row_property, col_property=distreader.col_property,order=order,dtype=dtype, val_shape=3)
        if hasattr(distreader,'val'):
            self.val[:,:,:] = distreader.val
        else:
            start = 0
            with log_in_place("DistMemMap writing sid_index ", logging.INFO) as updater:
                while start < distreader.sid_count:
                    updater('{0} of {1}'.format(start,distreader.sid_count))
                    distdata = distreader[:,start:start+block_size].read(order=order,dtype=dtype,num_threads=num_threads)
                    self.val[:,start:start+distdata.sid_count,:] = distdata.val
                    start += distdata.sid_count

        self.flush()
        if os.path.exists(filename):
           os.remove(filename) 
        shutil.move(filename+'.temp',filename)
        logging.debug("Done writing " + filename)
        return DistMemMap(filename)



    def _run_once(self):
            if (self._ran_once):
                return
            row_ascii,col_ascii,val,row_property,col_property = self._run_once_inner()
            row = np.array(row_ascii,dtype='str') #!!!avoid this copy when not needed
            col = np.array(col_ascii,dtype='str') #!!!avoid this copy when not needed

            DistData.__init__(self,iid=row,sid=col,val=val,pos=col_property,name="np.memmap('{0}')".format(self._filename))

class TestDistMemMap(unittest.TestCase):     

    def test1(self):        
        old_dir = os.getcwd()
        os.chdir(os.path.dirname(os.path.realpath(__file__)))

        filename2 = "tempdir/tiny.dist.memmap"
        pstutil.create_directory_if_necessary(filename2)
        distreader2 = DistMemMap.empty(iid=[['fam0','iid0'],['fam0','iid1']], sid=['snp334','snp349','snp921'],filename=filename2,order="F",dtype=np.float64)
        assert isinstance(distreader2.val,np.memmap)
        distreader2.val[:,:,:] = [[[.5,.5,0],[0,0,1],[.5,.5,0]],[[0,1.,0],[0,.75,.25],[.5,.5,0]]]
        assert np.array_equal(distreader2[[1],[1]].read(view_ok=True).val,np.array([[[0,.75,.25]]]))
        distreader2.flush()
        assert isinstance(distreader2.val,np.memmap)
        assert np.array_equal(distreader2[[1],[1]].read(view_ok=True).val,np.array([[[0,.75,.25]]]))
        distreader2.flush()

        distreader3 = DistMemMap(filename2)
        assert np.array_equal(distreader3[[1],[1]].read(view_ok=True).val,np.array([[[0,.75,.25]]]))
        assert isinstance(distreader3.val,np.memmap)

        logging.info("in TestDistMemMap test1")
        distreader = DistMemMap('../examples/tiny.dist.memmap')
        assert distreader.iid_count == 25
        assert distreader.sid_count == 10
        assert isinstance(distreader.val,np.memmap)

        distdata = distreader.read(view_ok=True)
        assert isinstance(distdata.val,np.memmap)
        os.chdir(old_dir)

    def test2(self):
        from pysnptools.distreader import Bgen

        old_dir = os.getcwd()
        os.chdir(os.path.dirname(os.path.realpath(__file__)))

        bgen = Bgen('../examples/example.bgen')
        distmemmap = DistMemMap.write("tempdir/bgentomemmap.dist.memamp",bgen)
        assert DistData.allclose(bgen.read(),distmemmap.read(),equal_nan=True)
        os.chdir(old_dir)

    def test_doctest(self):
        import pysnptools.distreader.distmemmap as mod_mm
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
    test_suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestDistMemMap))
    return test_suite


if __name__ == "__main__":
    logging.basicConfig(level=logging.WARN)

    suites = getTestSuite()
    r = unittest.TextTestRunner(failfast=True)
    ret = r.run(suites)
    assert ret.wasSuccessful()


    result = doctest.testmod(optionflags=doctest.ELLIPSIS|doctest.NORMALIZE_WHITESPACE)
    assert result.failed == 0, "failed doc test: " + __file__
