
import logging
import os
import shutil
import numpy as np
import unittest
import doctest
from pathlib import Path
import pysnptools.util as pstutil
from pysnptools.pstreader import PstData
from pysnptools.pstreader import PstMemMap
from pysnptools.snpreader import SnpReader, SnpData
from pysnptools.standardizer import Identity
from pysnptools.util import log_in_place

class SnpMemMap(PstMemMap,SnpData):
    '''
    A :class:`.SnpData` that keeps its data in a memory-mapped file. This allows data large than fits in main memory.

    See :class:`.SnpData` for general examples of using SnpData.

    **Constructor:**
        :Parameters: **filename** (*string*) -- The *\*.snp.memmap* file to read.
        
        Also see :meth:`.SnpMemMap.empty` and :meth:`.SnpMemMap.write`.

        :Example:

        >>> from pysnptools.snpreader import SnpMemMap
        >>> from pysnptools.util import example_file # Download and return local file name
        >>> mem_map_file = example_file('pysnptools/examples/tiny.snp.memmap')
        >>> snp_mem_map = SnpMemMap(mem_map_file)
        >>> print(snp_mem_map.val[0,1], snp_mem_map.iid_count, snp_mem_map.sid_count)
        2.0 2 3

    **Methods inherited from** :class:`.SnpData`

        :meth:`.SnpData.allclose`, :meth:`.SnpData.standardize`

    **Methods beyond** :class:`.SnpReader`

    '''

    def __init__(self, *args, **kwargs):
        super(SnpMemMap, self).__init__(*args, **kwargs)

    @property
    def val(self):
        """The 2D NumPy memmap array of floats that represents the values. You can get this property, but cannot set it (except with itself)

        >>> from pysnptools.snpreader import SnpMemMap
        >>> from pysnptools.util import example_file # Download and return local file name
        >>> mem_map_file = example_file('pysnptools/examples/tiny.snp.memmap')
        >>> snp_mem_map = SnpMemMap(mem_map_file)
        >>> print(snp_mem_map.val[0,1])
        2.0
        """
        self._run_once()
        return self._val


    @val.setter
    def val(self, new_value):
        self._run_once()
        if self._val is new_value:
            return
        raise Exception("SnpMemMap val's cannot be set to a different array")


    @property
    def offset(self):
        '''The byte position in the file where the memory-mapped values start.
       
        (The disk space before this is used to store :attr:`SnpReader.iid`, etc. information.
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
        '''Create an empty :class:`.SnpMemMap` on disk.

        :param iid: The :attr:`SnpReader.iid` information
        :type iid: an array of string pairs

        :param sid: The :attr:`SnpReader.sid` information
        :type sid: an array of strings

        :param filename: name of memory-mapped file to create
        :type filename: string

        :param pos: optional -- The additional :attr:`SnpReader.pos` information associated with each sid. Default: None
        :type pos: an array of numeric triples

        :param order: {'F' (default), 'C'}, optional -- Specify the order of the ndarray.
        :type order: string or None

        :param dtype: {numpy.float64 (default), numpy.float32}, optional -- The data-type for the :attr:`SnpMemMap.val` ndarray.
        :type dtype: data-type

        :rtype: :class:`.SnpMemMap`

        >>> import pysnptools.util as pstutil
        >>> from pysnptools.snpreader import SnpMemMap
        >>> filename = "tempdir/tiny.snp.memmap"
        >>> pstutil.create_directory_if_necessary(filename)
        >>> snp_mem_map = SnpMemMap.empty(iid=[['fam0','iid0'],['fam0','iid1']], sid=['snp334','snp349','snp921'],filename=filename,order="F",dtype=np.float64)
        >>> snp_mem_map.val[:,:] = [[0.,2.,0.],[0.,1.,2.]]
        >>> snp_mem_map.flush()

        '''

        self = SnpMemMap(filename)
        self._empty_inner(row=iid, col=sid, filename=filename, row_property=None, col_property=pos,order=order,dtype=dtype,val_shape=None)
        return self

    def flush(self):
        '''Flush :attr:`SnpMemMap.val` to disk and close the file. (If values or properties are accessed again, the file will be reopened.)

        >>> import pysnptools.util as pstutil
        >>> from pysnptools.snpreader import SnpMemMap
        >>> filename = "tempdir/tiny.snp.memmap"
        >>> pstutil.create_directory_if_necessary(filename)
        >>> snp_mem_map = SnpMemMap.empty(iid=[['fam0','iid0'],['fam0','iid1']], sid=['snp334','snp349','snp921'],filename=filename,order="F",dtype=np.float64)
        >>> snp_mem_map.val[:,:] = [[0.,2.,0.],[0.,1.,2.]]
        >>> snp_mem_map.flush()

        '''
        if self._ran_once:
            self.val.flush()
            del self._val
            self._ran_once = False


    @staticmethod
    def write(filename, snpreader, standardizer=Identity(), order='A', dtype=None, block_size=None, num_threads=None):
        """Writes a :class:`SnpReader` to :class:`SnpMemMap` format.

        :param filename: the name of the file to create
        :type filename: string
        :param snpreader: The data that should be written to disk.
        :type snpreader: :class:`SnpReader`
        :rtype: :class:`.SnpMemMap`

        >>> import pysnptools.util as pstutil
        >>> from pysnptools.util import example_file # Download and return local file name
        >>> from pysnptools.snpreader import Bed, SnpMemMap
        >>> bed_file = example_file("pysnptools/examples/toydata.5chrom.*","*.bed")
        >>> bed = Bed(bed_file)
        >>> pstutil.create_directory_if_necessary("tempdir/toydata.5chrom.snp.memmap") #LATER should we just promise to create directories?
        >>> SnpMemMap.write("tempdir/toydata.5chrom.snp.memmap",bed)      # Write bed in SnpMemMap format
        SnpMemMap('tempdir/toydata.5chrom.snp.memmap')
        """
        block_size = block_size or max((100_000)//max(1,snpreader.row_count),1)

        if hasattr(snpreader,'val'):
            order = PstMemMap._order(snpreader) if order=='A' else order
            dtype = dtype or snpreader.val.dtype
        else:
            order = 'F' if order=='A' else order
            dtype = dtype or np.float64
        dtype = np.dtype(dtype)

        snpmemmap = SnpMemMap.empty(iid=snpreader.iid, sid=snpreader.sid, filename=filename+'.temp', pos=snpreader.col_property,order=order,dtype=dtype)
        if hasattr(snpreader,'val'):
            standardizer.standardize(snpreader,num_threads=num_threads)
            snpmemmap.val[:,:] = snpreader.val
        else:
            with log_in_place("SnpMemMap write sid_index ", logging.INFO) as updater:
                for start in range(0,snpreader.sid_count,block_size):
                    updater('{0} of {1}'.format(start,snpreader.sid_count))
                    snpdata = snpreader[:,start:start+block_size].read(order=order,dtype=dtype,num_threads=num_threads)
                    standardizer.standardize(snpdata,num_threads=num_threads)
                    snpmemmap.val[:,start:start+snpdata.sid_count] = snpdata.val

        snpmemmap.flush()
        if os.path.exists(filename):
           os.remove(filename) 
        shutil.move(filename+'.temp',filename)
        logging.debug("Done writing " + filename)
        return SnpMemMap(filename)


    def _run_once(self):
            if (self._ran_once):
                return
            row_ascii,col_ascii,val,row_property,col_property = self._run_once_inner()
            row = np.array(row_ascii,dtype='str') #!!!avoid this copy when not needed
            col = np.array(col_ascii,dtype='str') #!!!avoid this copy when not needed

            SnpData.__init__(self,iid=row,sid=col,val=val,pos=col_property,name="np.memmap('{0}')".format(self._filename))

class TestSnpMemMap(unittest.TestCase):     

    def test1(self):
        from pysnptools.snpreader import Bed, SnpMemMap
        from pysnptools.util import example_file # Download and return local file name

        old_dir = os.getcwd()
        os.chdir(os.path.dirname(os.path.realpath(__file__)))

        filename2 = "tempdir/tiny.snp.memmap"
        pstutil.create_directory_if_necessary(filename2)
        snpreader2 = SnpMemMap.empty(iid=[['fam0','iid0'],['fam0','iid1']], sid=['snp334','snp349','snp921'],filename=filename2,order="F",dtype=np.float64)
        assert isinstance(snpreader2.val,np.memmap)
        snpreader2.val[:,:] = [[0.,2.,0.],[0.,1.,2.]]
        assert np.array_equal(snpreader2[[1],[1]].read(view_ok=True).val,np.array([[1.]]))
        snpreader2.flush()
        assert isinstance(snpreader2.val,np.memmap)
        assert np.array_equal(snpreader2[[1],[1]].read(view_ok=True).val,np.array([[1.]]))
        snpreader2.flush()

        snpreader3 = SnpMemMap(filename2)
        assert np.array_equal(snpreader3[[1],[1]].read(view_ok=True).val,np.array([[1.]]))
        assert isinstance(snpreader3.val,np.memmap)

        logging.info("in TestSnpMemMap test1")
        snpreader = SnpMemMap('tempdir/tiny.snp.memmap')
        assert snpreader.iid_count == 2
        assert snpreader.sid_count == 3
        assert isinstance(snpreader.val,np.memmap)

        snpdata = snpreader.read(view_ok=True)
        assert isinstance(snpdata.val,np.memmap)

        bed_file = example_file("pysnptools/examples/toydata.5chrom.*","*.bed")
        bed = Bed(bed_file)
        pstutil.create_directory_if_necessary("tempdir/toydata.5chrom.snp.memmap") #LATER should we just promise to create directories?
        SnpMemMap.write("tempdir/toydata.5chrom.snp.memmap",bed)      # Write bed in SnpMemMap format
        SnpMemMap.write("tempdir/toydata.5chromsnpdata.snp.memmap",bed[:,::2].read())  # Write snpdata in SnpMemMap format

        os.chdir(old_dir)


def getTestSuite():
    """
    set up composite test suite
    """
    
    test_suite = unittest.TestSuite([])
    test_suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestSnpMemMap))
    return test_suite


if __name__ == "__main__":
    logging.basicConfig(level=logging.WARN)

    if False:
        from pysnptools.snpreader import Bed, _MergeSIDs, SnpMemMap

        # Data will appear in the memory mapped file in the order given.
        # The *.fam and *.bim files must be in the same order as the bed files.
        bed_file_list = []
        fam_file_list = []
        bim_file_list = []
        for piece in range(25):
             bed_file_list += [r"M:\deldir\testsnps_1_10_250000_10000\chrom10.piece{0}of25.bed".format(piece)]
             fam_file_list += [r"M:\deldir\testsnps_1_10_250000_10000\chrom10.piece{0}of25.fam".format(piece)]
             bim_file_list += [r"M:\deldir\testsnps_1_10_250000_10000\chrom10.piece{0}of25.bim".format(piece)]
        #for chrom in range(21,23):
        #    bed_file_list += [r"d:\deldir\genbgen\merged_487400x220000.{0}.bed".format(chrom)]
        #    fam_file_list += [r"m:\deldir\genbgen\merged_487400x220000.{0}.fam".format(chrom)]
        #    bim_file_list += [r"m:\deldir\genbgen\merged_487400x220000.{0}.bim".format(chrom)]

        memmap_file = r"D:\deldir\memmap1.snp.memmap"

        #######
        # For this demo, erase the memmap output file
        ######

        if Path(memmap_file).exists():
            Path(memmap_file).unlink()

    
        #######
        # Merge the input files
        ######
        merge = _MergeSIDs(
            [Bed(bed_file,fam_filename=fam_file,bim_filename=bim_file,count_A1=True,skip_format_check=True)
             for bed_file, fam_file, bim_file in zip(bed_file_list,fam_file_list,bim_file_list)]
            )


        # memmap = _bed_to_memmap2(merge,memmap_file=memmap_file,dtype='float32',step=10)
        from pysnptools.standardizer import Unit
        memmap = SnpMemMap.write(memmap_file, merge, standardizer=Unit(),dtype='float32')
        memmap

    suites = getTestSuite()
    r = unittest.TextTestRunner(failfast=True)
    ret = r.run(suites)
    assert ret.wasSuccessful()

    result = doctest.testmod(optionflags=doctest.ELLIPSIS)
    assert result.failed == 0, "failed doc test: " + __file__
