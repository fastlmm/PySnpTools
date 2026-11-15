import os
import numpy as np
import logging
import doctest
import unittest

class TestUtilTools(unittest.TestCase):

    def test_util_mapreduce1_testmod(self):
        import pysnptools.util.mapreduce1
        import pysnptools.util.pheno
        import pysnptools.util.mapreduce1.examples
        old_dir = os.getcwd()
        os.chdir(os.path.dirname(os.path.realpath(__file__))+"/mapreduce1")
        for mod in [
                    pysnptools.util,
                    pysnptools.util.pheno,
                    pysnptools.util.mapreduce1.examples,
                    pysnptools.util.mapreduce1.mapreduce,
                    pysnptools.util.mapreduce1.runner.local,
                    pysnptools.util.mapreduce1.runner.localinparts,
                    pysnptools.util.mapreduce1.runner.localmultiproc,
                    pysnptools.util.mapreduce1.runner.localmultithread,
                    ]:
            result = doctest.testmod(mod,optionflags=doctest.ELLIPSIS|doctest.NORMALIZE_WHITESPACE)
            assert result.failed == 0, "failed doc test: " + __file__
        os.chdir(old_dir)

    def test_local(self):
        from pysnptools.util.mapreduce1 import map_reduce
        from pysnptools.util.mapreduce1.runner import Local
        def holder1(n,runner):
            def mapper1(x):
                return x*x
            def reducer1(sequence):
               return sum(sequence)
            return map_reduce(range(n),mapper=mapper1,reducer=reducer1,runner=runner)
        assert 328350 == holder1(100,Local())

    def test_localinparts(self):
        from pysnptools.util.mapreduce1 import map_reduce
        from pysnptools.util.mapreduce1.runner import LocalInParts
        def holder1(n,runner):
            def mapper1(x):
                return x*x
            def reducer1(sequence):
               return sum(sequence)
            return map_reduce(range(n),mapper=mapper1,reducer=reducer1,runner=runner)
        holder1(100,LocalInParts(0,4)) #Run part 0 of 4 and save partial results to disk as '0.4.p'.
        holder1(100,LocalInParts(1,4)) #Run part 1 of 4 and save partial results to disk as '1.4.p'.
        holder1(100,LocalInParts(2,4)) #Run part 2 of 4 and save partial results to disk as '2.4.p'.
        holder1(100,LocalInParts(3,4)) #Run part 3 of 4 and save partial results to disk as '3.4.p'.
        assert 328350 == holder1(100,LocalInParts(4,4)) #Read the all partial results and then apply the reduce function and return the result.

    def test_localmultiproc(self):
        from pysnptools.util.mapreduce1 import map_reduce
        from pysnptools.util.mapreduce1.runner import LocalMultiProc
        def holder1(n,runner):
            def mapper1(x):
                return x*x
            def reducer1(sequence):
               return sum(sequence)
            return map_reduce(range(n),mapper=mapper1,reducer=reducer1,runner=runner)
        assert 328350 == holder1(100,LocalMultiProc(4))

    def test_localmultithread(self):
        from pysnptools.util.mapreduce1 import map_reduce
        from pysnptools.util.mapreduce1.runner import LocalMultiThread
        def holder1(n,runner):
            def mapper1(x):
                return x*x
            def reducer1(sequence):
               return sum(sequence)
            return map_reduce(range(n),mapper=mapper1,reducer=reducer1,runner=runner)
        assert 328350 == holder1(100,LocalMultiThread(4))

    def test_localinparts_with_weights(self):
        from pysnptools.util.mapreduce1 import map_reduce
        from pysnptools.util.mapreduce1.runner import LocalInParts
        def holder1(n,runner):
            def mapper1(x):
                return int(os.environ['TEST_ENVIRON'])
            def reducer1(sequence):
               return sum(sequence)+int(os.environ['TEST_ENVIRON'])
            return map_reduce(range(n),mapper=mapper1,reducer=reducer1,runner=runner)
        weights = [1,97,1,1]
        for i in range(5):
            environ = {'TEST_ENVIRON':str(i)}
            assert 'TEST_ENVIRON' not in os.environ
            result = holder1(100,LocalInParts(i,4,weights=weights,environ=environ)) #ignore all but the last result
            assert 'TEST_ENVIRON' not in os.environ
        assert result == 0*1+1*97+2*1+3*1+4

    def test_localinmultiproc_with_weights(self):
        from pysnptools.util.mapreduce1 import map_reduce
        from pysnptools.util.mapreduce1.runner import LocalMultiProc
        def holder1(n,runner):
            def mapper1(x):
                return int(os.environ['TEST_ENVIRON'])
            def reducer1(sequence):
               return sum(sequence)+int(os.environ['TEST_ENVIRON'])
            return map_reduce(range(n),mapper=mapper1,reducer=reducer1,runner=runner)
        weights = [1,97,1,1]
        def taskindex_to_environ(taskindex):
            return {'TEST_ENVIRON': str(taskindex)}
        runner = LocalMultiProc(4,weights=weights,taskindex_to_environ=taskindex_to_environ)
        assert 'TEST_ENVIRON' not in os.environ
        result = holder1(100,runner)
        assert result == 0*1+1*97+2*1+3*1+4

    def test_sub_matrix(self):
        import pysnptools.util as pstutil
        np.random.seed(0) # set seed so that results are deterministic
        matrix = np.random.rand(12,7) # create a 12 x 7 ndarray
        submatrix = pstutil.sub_matrix(matrix,[0,2,11],[6,5,4,3,2,1,0])
        assert  matrix[2,0] == submatrix[1,6] #The row # 2 is now #1, the column #0 is now #6.

        np.random.seed(0) # set seed so that results are deterministic
        matrix = np.random.rand(12,7,3) # create a 12 x 7 ndarray
        submatrix = pstutil.sub_matrix(matrix,[0,2,11],[6,5,4,3,2,1,0])
        assert  matrix[2,0,1] == submatrix[1,6,1] #The row # 2 is now #1, the column #0 is now #6.

    def test_sub_matrix_cupy(self):
        # cupy is optional; this test needs it and a GPU.
        try:
            import cupy as cp
        except ModuleNotFoundError:
            raise unittest.SkipTest("cupy is not installed")
        import warnings
        import pysnptools.util as pstutil
        rng = np.random.default_rng(0)
        row_index_list = [0, 2, 11, 2]
        col_index_list = [6, 5, 4, 3, 2, 1, 0]
        for shape in [(12, 7), (12, 7, 3)]:
            for in_order in ["F", "C"]:
                for in_dtype in [np.float64, np.float32]:
                    matrix = np.asarray(rng.random(shape), dtype=in_dtype, order=in_order)
                    for out_order in ["A", "F", "C"]:
                        for out_dtype in [np.float64, np.float32]:
                            with warnings.catch_warnings():
                                warnings.simplefilter("ignore", UserWarning)
                                expected = pstutil.sub_matrix(matrix, row_index_list, col_index_list, order=out_order, dtype=out_dtype)
                                actual = pstutil.sub_matrix(cp.asarray(matrix, order=in_order), row_index_list, col_index_list, order=out_order, dtype=out_dtype)
                            assert isinstance(actual, cp.ndarray)
                            assert actual.dtype == expected.dtype
                            assert actual.flags["F_CONTIGUOUS"] == expected.flags["F_CONTIGUOUS"]
                            assert actual.flags["C_CONTIGUOUS"] == expected.flags["C_CONTIGUOUS"]
                            np.testing.assert_array_equal(cp.asnumpy(actual), expected)

    def test_snpdata_subset_cupy(self):
        # End to end: SnpData created under ARRAY_MODULE=cupy holds a cupy val and must subset.
        try:
            import cupy as cp
        except ModuleNotFoundError:
            raise unittest.SkipTest("cupy is not installed")
        from unittest.mock import patch
        from pysnptools.snpreader import SnpData
        val = np.arange(24, dtype=np.float64).reshape(6, 4)
        with patch.dict("os.environ", {"ARRAY_MODULE": "cupy"}):
            snpdata = SnpData(iid=[["f", str(i)] for i in range(6)], sid=[f"s{j}" for j in range(4)], val=val)
            assert isinstance(snpdata.val, cp.ndarray)
            subset = snpdata[np.array([True, False, True, True, False, True]), [3, 1]].read()
        np.testing.assert_array_equal(cp.asnumpy(subset.val), val[[0, 2, 3, 5]][:, [3, 1]])

    def test_snpdata_cupy_keeps_float32(self):
        # A numpy float32 val moved to cupy must stay float32 (not be upcast to float64).
        try:
            import cupy as cp
        except ModuleNotFoundError:
            raise unittest.SkipTest("cupy is not installed")
        from unittest.mock import patch
        from pysnptools.snpreader import SnpData
        iid = [["f", str(i)] for i in range(3)]
        sid = ["s0", "s1"]
        with patch.dict("os.environ", {"ARRAY_MODULE": "cupy"}):
            for in_dtype, out_dtype in [(np.float32, np.float32), (np.float64, np.float64), (np.int8, np.float64)]:
                snpdata = SnpData(iid=iid, sid=sid, val=np.ones((3, 2), dtype=in_dtype))
                assert isinstance(snpdata.val, cp.ndarray)
                assert snpdata.val.dtype == out_dtype, (in_dtype, snpdata.val.dtype)



def getTestSuite():
    """
    set up composite test suite
    """

    test_suite = unittest.TestSuite([])
    test_suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestUtilTools))
    return test_suite

if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)

    suites = getTestSuite()
    r = unittest.TextTestRunner(failfast=False)
    ret = r.run(suites)
    assert ret.wasSuccessful()
