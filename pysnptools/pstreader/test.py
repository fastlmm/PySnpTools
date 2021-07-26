import numpy as np
import scipy as sp
import logging
import doctest
import unittest
import os.path
import time
import sys
import tempfile
from pysnptools.pstreader import PstData, PstNpz, PstHdf5, PstMemMap
from pysnptools.util import create_directory_if_necessary
from pysnptools.kernelreader.test import _fortesting_JustCheckExists


class TestPstReader(unittest.TestCase):     

    def test_big_npz(self):
        logging.info("in test_big_npz")
        n = 1000
        pstdata = PstData(row=range(n-1),col=range(n+1),val=np.zeros([n-1,n+1]))
        output = "tempdir/pstreader/big.npz"
        create_directory_if_necessary(output)
        PstNpz.write(output,pstdata)
        pstnpz = PstNpz(output)
        pstdata1 = pstnpz[::2,::4].read()
        pstdata2 = pstnpz.read(order='A')
        assert pstdata2.val.flags['C_CONTIGUOUS']

        pstdata = PstData(row=range(n-1),col=range(n+1),val=np.zeros([n-1,n+1],order='F'))
        PstNpz.write(output,pstdata)
        pstnpz = PstNpz(output)
        pstdata2 = pstnpz.read(order='A')
        pstdata2.val.flags['F_CONTIGUOUS']

    def test_big_npz_3d(self):
        logging.info("in test_big_npz_3d")
        n = 1000
        pstdata = PstData(row=range(n-1),col=range(n+1),val=np.zeros([n-1,n+1,3]))
        output = "tempdir/pstreader/big3d.npz"
        create_directory_if_necessary(output)
        PstNpz.write(output,pstdata)
        pstnpz = PstNpz(output)
        pstdata1 = pstnpz[::2,::4].read()
        pstdata2 = pstnpz.read(order='A')
        assert pstdata2.val.flags['C_CONTIGUOUS']

        pstdata = PstData(row=range(n-1),col=range(n+1),val=np.zeros([n-1,n+1,3],order='F'))
        PstNpz.write(output,pstdata)
        pstnpz = PstNpz(output)
        pstdata2 = pstnpz.read(order='A')
        pstdata2.val.flags['F_CONTIGUOUS']

    def test_writes(self):
        #===================================
        #    Defining sub functions
        #===================================
        def _oned_int(c):
            return range(c)
        def _oned_str(c):
            return [str(i).encode('ascii') for i in range(c)]
        def _twooned_int(c):
            return [[i] for i in range(c)]
        def _twooned_str(c):
            return [[str(i).encode('ascii')] for i in range(c)]
        def _twotwod_int(c):
            return [[i,i] for i in range(c)]
        def _twotwod_str(c):
            return [[str(i).encode('ascii'),b"hello"] for i in range(c)]
        #def _twotwod_U(c):
        #    return [[str(i).encode('UTF-8'),u"hello"] for i in range(c)]
        def _none(c):
            return None
        def _zero(c):
            return np.empty([c,0],dtype='S')
        #===================================
        #    Starting main function
        #===================================
        logging.info("starting 'test_writes'")
        np.random.seed(0)
        temp_dir = tempfile.TemporaryDirectory("pstreader")
        output_template = temp_dir.name + '/writes.{0}.{1}'
        i = 0
        for row_count in [5,2,1,0]:
            for col_count in [4,2,1,0]:
                for val_shape in [3,None,1]:
                    val = np.random.normal(.5,2,size=(row_count,col_count)) if val_shape is None else np.random.normal(.5,2,size=(row_count,col_count,val_shape))
                    for row_or_col_gen in [_oned_int,_oned_str,_twooned_int,_twooned_str,_twotwod_int,_twotwod_str]:#!!!,_twotwod_U can't roundtrop Unicode in hdf5
                        row = row_or_col_gen(row_count)
                        col = row_or_col_gen(col_count)
                        for prop_gen in [_none,_oned_str,_oned_int,_twooned_int,_twooned_str,_twotwod_int,_twotwod_str,_zero]: #!!!_twotwod_U can't round trip Unicode because Hdf5 doesn't like it.
                            row_prop = prop_gen(row_count)
                            col_prop = prop_gen(col_count)
                            pstdata = PstData(row,col,val,row_prop,col_prop,str(i))
                            for the_class,suffix in [(PstMemMap,"memmap"),(PstHdf5,"hdf5"),(PstNpz,"npz")]:
                                filename = output_template.format(i,suffix)
                                logging.info(filename)
                                i += 1
                                the_class.write(filename,pstdata)
                                reader = the_class(filename) if suffix!='hdf5' else the_class(filename,block_size=3)
                                _fortesting_JustCheckExists().input(reader)
                                for subsetter in [None, np.s_[::2,::3]]:
                                    subreader = reader if subsetter is None else reader[subsetter[0],subsetter[1]]
                                    expected = pstdata if subsetter is None else pstdata[subsetter[0],subsetter[1]].read()
                                    for order in ['C','F','A']:
                                        for force_python_only in [True,False]:
                                            readdata = subreader.read(order=order,force_python_only=force_python_only)
                                            assert np.array_equal(readdata.val,expected.val)
                                            assert np.array_equal(readdata.row,expected.row)
                                            assert np.array_equal(readdata.col,expected.col)
                                            assert np.array_equal(readdata.row_property,expected.row_property) or (readdata.row_property.shape[1]==0 and expected.row_property.shape[1]==0)
                                            assert np.array_equal(readdata.col_property,expected.col_property) or (readdata.col_property.shape[1]==0 and expected.col_property.shape[1]==0)
                                if suffix in {'memmap','hdf5'}:
                                    reader.flush()
                                os.remove(filename)
        temp_dir.cleanup()
        logging.info("done with 'test_writes'")

    def test_every_read(self):
        for order_from in ['F','C']:
            for order_to in ['F','C']:
                for dtype_from in [np.float32,np.float64]:
                    for dtype_to in [np.float32,np.float64]:
                        for val_shape in [None,1,3]:
                            for force_python_only in [True,False]:
                                np.random.seed(0)
                                val0 = np.random.normal(.5,2,size=(3,2)) if val_shape is None else np.random.normal(.5,2,size=(3,2,val_shape))
                                val = np.array(val0,order=order_from,dtype=dtype_from)
                                pstdata = PstData(val=val,row=list(range(3)),col=list(range(2)))
                                expected = np.array(val[::-2,:][:,::-1],order=order_to,dtype=dtype_to)
                                result = pstdata[::-2,::-1].read(order=order_to,dtype=dtype_to,force_python_only=force_python_only)
                                assert result.val.dtype==dtype_to
                                assert (order_to == 'F' and result.val.flags['F_CONTIGUOUS']) or (order_to == 'C' and result.val.flags['C_CONTIGUOUS'])
                                assert np.array_equal(result.val,expected)

    def test_repr_test(self):
        np.random.seed(0)
        row_property=np.array([[1.0,2,2.5],[3,4,4.5],[5,6,6.5]])
        col_property=np.array([[1.0,2,2.5,1],[3,4,4.5,3]])
        pstdata = PstData(row=np.array([[1.0,2],[3,4],[5,6]]),
                          col=np.array([("A","a"),("B","b")]),
                          val = np.random.normal(.5,2,size=(3,2)),
                          row_property=row_property,
                          col_property=col_property)
        assert pstdata.col_to_index([("B","b")])[0] == 1
        s = str(pstdata)

    def test_read(self):
        np.random.seed(0)
        row_property=np.array([[1.0,2,2.5],[3,4,4.5],[5,6,6.5]])
        col_property=np.array([[1.0,2,2.5,1],[3,4,4.5,3]])
        pstdata = PstData(row=np.array([[1.0,2],[3,4],[5,6]]),
                          col=np.array([["A","a"],["B","b"]]),
                          val = np.random.normal(.5,2,size=(3,2)),
                          row_property=row_property,
                          col_property=col_property,
                          name="test_read")

        assert pstdata.row_to_index([np.array([3.0,4])])[0] == 1
        assert pstdata.col_to_index([np.array(["A","a"])])[0] == 0
        assert np.array_equal(pstdata[1:,:2].row_property,row_property[1:])
        assert np.array_equal(pstdata[1:,:2].col_property,col_property[:2])


        pstdata2 = pstdata[:2,:2].read()
        from pysnptools.kernelreader.test import _fortesting_JustCheckExists
        _fortesting_JustCheckExists().input(pstdata)
        _fortesting_JustCheckExists().input(pstdata2)

        np.testing.assert_array_almost_equal(pstdata2.val, pstdata.val[:2,:2], decimal=10)
        pstdata3 = pstdata[[],:].read()
        assert pstdata3.val.shape[0] == 0 and pstdata3.val.shape[1]==2
        pstdata.val = pstdata.val.copy(order='F')
        pstdata2 = pstdata[:2,:2].read()
        np.testing.assert_array_almost_equal(pstdata2.val, pstdata.val[:2,:2], decimal=10)
        pstdata2 = pstdata[:2,:2].read(order='F')
        np.testing.assert_array_almost_equal(pstdata2.val, pstdata.val[:2,:2], decimal=10)
        pstdata2 = pstdata[:2,:2].read(order='A')
        np.testing.assert_array_almost_equal(pstdata2.val, pstdata.val[:2,:2], decimal=10)
        pstdata2 = pstdata[:2,:2].read(force_python_only=True,dtype=None,order='C')
        np.testing.assert_array_almost_equal(pstdata2.val, pstdata.val[:2,:2], decimal=10)
        pstdata2 = pstdata[:2,:2].read(force_python_only=True,dtype='float32',order='C')
        np.testing.assert_array_almost_equal(pstdata2.val, pstdata.val[:2,:2].astype(dtype='float32'), decimal=10)
        pstdata2 = pstdata[:2,:2].read(force_python_only=True,dtype='float32',order=None)
        np.testing.assert_array_almost_equal(pstdata2.val, pstdata.val[:2,:2].astype(dtype='float32'), decimal=10)
        pstdata2 = pstdata[:2,:2].read(force_python_only=True,dtype=None,order='F')
        np.testing.assert_array_almost_equal(pstdata2.val, pstdata.val[:2,:2], decimal=10)
        pstdata4 = pstdata[::,::].read(force_python_only=True)
        np.testing.assert_array_almost_equal(pstdata4.val, pstdata.val, decimal=10)


        logging.info("done with test")


    def test_inputs(self):
        from pysnptools.pstreader import PstData
        np.random.seed(0)
        row_property=np.array([1.0,2,2.5])
        col_property=np.array([1,2])
        pstdata = PstData(row=np.array([1.0,3,6]),
                          col=np.array(["Aa","Bb"]),
                          val = np.random.normal(.5,2,size=(3,2)),
                          row_property=row_property,
                          col_property=col_property,
                          name="test_read")

        assert pstdata.row_to_index([3])[0] == 1
        assert pstdata.col_to_index(["Aa"])[0] == 0
        assert np.array_equal(pstdata[1:,:2].row_property,row_property[1:])
        assert np.array_equal(pstdata[1:,:2].col_property,col_property[:2])
        logging.info("done with test")


    def test_inputs2(self):
        from pysnptools.pstreader import PstData
        np.random.seed(0)
        row_property=None
        col_property=None
        pstdata = PstData(row=np.array([1.0,3,6]),
                          col=np.array(["Aa","Bb"]),
                          val = np.random.normal(.5,2,size=(3,2)),
                          row_property=row_property,
                          col_property=col_property,
                          name="test_read")

        assert pstdata.row_to_index([3])[0] == 1
        assert pstdata.col_to_index(["Aa"])[0] == 0
        assert np.array_equal(pstdata[1:,:2].row_property,pstdata.row_property[1:])
        assert np.array_equal(pstdata[1:,:2].col_property,pstdata.col_property[:2])
        logging.info("done with test")

    def test_inputs3(self):
        from pysnptools.pstreader import PstData
        np.random.seed(0)
        row_property=None
        col_property=None
        pstdata = PstData(row=[[1.0,2.0],[3,4],[6,7]],
                          col=np.array([]),
                          val = [[],[],[]],
                          row_property=row_property,
                          col_property=col_property,
                          name="test_read")

        assert pstdata.row_to_index([[3,4]])[0] == 1
        assert np.array_equal(pstdata[1:,:2].row_property,pstdata.row_property[1:])
        assert np.array_equal(pstdata[1:,:2].col_property,pstdata.col_property[:2])
        logging.info("done with test")

    def test_flush(self):

        previous_wd = os.getcwd()
        os.chdir(os.path.dirname(os.path.realpath(__file__)))

        reader = PstHdf5('../examples/toydata.kernel.hdf5')
        val1 = reader[0,0].read()
        reader.flush()
        val2 = reader[0,0].read()
        assert val1 == val2

    def test_respect_read_inputs(self):
        from pysnptools.pstreader import _MergeRows,_MergeCols

        previous_wd = os.getcwd()
        os.chdir(os.path.dirname(os.path.realpath(__file__)))

        for pstreader in [
                           PstNpz('../examples/toydata10.snp.npz'),
                           _MergeRows([PstHdf5('../examples/toydata.snpmajor.snp.hdf5')[:5,:].read(),PstHdf5('../examples/toydata.snpmajor.snp.hdf5')[5:,:].read()]),
                           _MergeCols([PstHdf5('../examples/toydata.snpmajor.snp.hdf5')[:,:5].read(),PstHdf5('../examples/toydata.snpmajor.snp.hdf5')[:,5:].read()]),
                           PstHdf5('../examples/toydata.snpmajor.snp.hdf5')[::2,::2],
                           PstHdf5('../examples/toydata.snpmajor.dist.hdf5').read(),
                           PstHdf5('../examples/toydata.kernel.hdf5'),
                           PstMemMap('../examples/tiny.pst.memmap')
                          ]:
            logging.info(str(pstreader))
            for order in ['F','C','A']:
                for dtype in [np.float32,np.float64]:
                    for force_python_only in [True,False]:
                        for view_ok in [True,False]:
                            val = pstreader.read(order=order,dtype=dtype,force_python_only=force_python_only,view_ok=view_ok).val
                            has_right_order = order=="A" or (order=="C" and val.flags["C_CONTIGUOUS"]) or (order=="F" and val.flags["F_CONTIGUOUS"])
                            if hasattr(pstreader,'val') and not view_ok:
                                assert pstreader.val is not val
                            if (hasattr(pstreader,'val') and view_ok and pstreader.val is not val and
                                (order == 'A' or (order == 'F' and pstreader.val.flags['F_CONTIGUOUS']) or (order == 'C' and pstreader.val.flags['C_CONTIGUOUS'])) and
                                (dtype is None or  pstreader.val.dtype == dtype)):
                                logging.info("{0} could have read a view, but didn't".format(pstreader))
                            assert val.dtype == dtype and has_right_order

        os.chdir(previous_wd)

    def test_inputs4(self):
        from pysnptools.pstreader import PstData
        pstdata = PstData(row=None,
                          col=None,
                          val = None,
                          row_property=None,
                          col_property=None,
                          name="test_read")

        assert pstdata.row_count == 0 and pstdata.col_count == 0 and pstdata.val.shape[0] == 0 and pstdata.val.shape[1]==0 and len(pstdata.row_property)==0 and len(pstdata.col_property)==0 
        logging.info("done with test")

# We do it this way instead of using doctest.DocTestSuite because doctest.DocTestSuite requires modules to be pickled, which python doesn't allow.
# We need tests to be pickleable so that they can be run on a cluster.
class TestPstDocStrings(unittest.TestCase):
    pass

    def test_pstdata(self):
        import pysnptools.pstreader.pstdata
        old_dir = os.getcwd()
        os.chdir(os.path.dirname(os.path.realpath(__file__)))
        result = doctest.testmod(pysnptools.pstreader.pstdata,optionflags=doctest.ELLIPSIS)
        os.chdir(old_dir)
        assert result.failed == 0, "failed doc test: " + __file__

    def test_psthdf5(self):
        import pysnptools.pstreader.psthdf5
        old_dir = os.getcwd()
        os.chdir(os.path.dirname(os.path.realpath(__file__)))
        result = doctest.testmod(pysnptools.pstreader.psthdf5,optionflags=doctest.ELLIPSIS)
        os.chdir(old_dir)
        assert result.failed == 0, "failed doc test: " + __file__

    def test_pstmemmap(self):
        import pysnptools.pstreader.pstmemmap
        old_dir = os.getcwd()
        os.chdir(os.path.dirname(os.path.realpath(__file__)))
        result = doctest.testmod(pysnptools.pstreader.pstmemmap,optionflags=doctest.ELLIPSIS)
        os.chdir(old_dir)
        assert result.failed == 0, "failed doc test: " + __file__

    def test_pstnpz(self):
        import pysnptools.pstreader.pstnpz
        old_dir = os.getcwd()
        os.chdir(os.path.dirname(os.path.realpath(__file__)))
        result = doctest.testmod(pysnptools.pstreader.pstnpz,optionflags=doctest.ELLIPSIS)
        os.chdir(old_dir)
        assert result.failed == 0, "failed doc test: " + __file__

    def test_pstreader(self):
        import pysnptools.pstreader.pstreader
        old_dir = os.getcwd()
        os.chdir(os.path.dirname(os.path.realpath(__file__)))
        result = doctest.testmod(pysnptools.pstreader.pstreader,optionflags=doctest.ELLIPSIS|doctest.NORMALIZE_WHITESPACE)
        os.chdir(old_dir)
        assert result.failed == 0, "failed doc test: " + __file__

    def test_subset(self):
        np.random.seed(0)
        row_property=np.array([[1.0,2,2.5],[3,4,4.5],[5,6,6.5]])
        col_property=np.array([[1.0,2,2.5,1],[3,4,4.5,3]])
        val = np.random.normal(.5,2,size=(3,2))
        pstdata = PstData(row=np.array([[1.0,2],[3,4],[5,6]]),
                          col=np.array([["A","a"],["B","b"]]),
                          val = val,
                          row_property=row_property,
                          col_property=col_property,
                          name="test_read")

        assert np.array_equal(pstdata[-1:0:-1,:].read().val, val[-1:0:-1,:])
        assert pstdata[-1,-1].read().val[0,0] == val[-1,-1]
        assert np.array_equal(pstdata[[-1,0],[-1,0]].read().val,val[[-1,0],:][:,[-1,0]])
        assert np.array_equal(pstdata[[True,False,True],[False,True]].read().val,val[[True,False,True],[False,True]].reshape(2,1))
        assert pstdata[0,0].read().val[0,0] == val[0,0]
        assert np.array_equal(pstdata[1::2,1::2].read().val,val[1::2,1::2])

        logging.info("done with test")


def getTestSuite():
    """
    set up composite test suite
    """
    
    test_suite = unittest.TestSuite([])
    test_suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestPstDocStrings))
    test_suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestPstReader))
    return test_suite

if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)

    suites = getTestSuite()
    r = unittest.TextTestRunner(failfast=False)
    ret = r.run(suites)
    assert ret.wasSuccessful()

