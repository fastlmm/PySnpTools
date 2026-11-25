import os
import numpy as np
import logging
import doctest
import unittest
import warnings
from pysnptools.util.gpu import GpuFloat64Warning

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

    def test_sub_matrix_indices(self):
        self._check_sub_matrix_indices(np)

    def test_sub_matrix_indices_cupy(self):
        try:
            import cupy as cp
        except ModuleNotFoundError:
            raise unittest.SkipTest("cupy is not installed")
        self._check_sub_matrix_indices(cp)

    def _check_sub_matrix_indices(self, xp):
        import pysnptools.util as pstutil

        invalid_indices = [
            [0.9], [-0.9], [1.0], [0, 0.9],
            [np.nan], [np.inf], [-np.inf], [1 + 0j], ["1"],
            [True], [False, 1], (np.bool_(True), 1),
            np.array([True, False]), np.array([0.9], dtype=object),
            [-1], [3], np.array([2**64 - 1], dtype=np.uint64),
            [2**64], [2**100], [-2**100],
            0, [[0]], np.empty((0, 1), dtype=np.int64),
            # The vectorized list/tuple/range path must reject these too.
            range(-1, 2), range(4), (0, True), [0, 2**63],
        ]
        valid_indices = [
            [2, 0, 2], (2, 0, 2), range(2, -1, -1),
            np.array([2, 0, 2], dtype=np.int32),
            np.array([2, 0, 2], dtype=np.uint64),
            np.array([2, 0, 2], dtype=object),
            [np.int64(2), np.uint64(0)], [], range(0),
        ]
        for shape in [(3, 3), (3, 3, 2)]:
            for order in ["C", "F"]:
                for dtype in [np.float32, np.float64]:
                    values = np.arange(np.prod(shape), dtype=dtype).reshape(shape)
                    matrix = xp.array(values, order=order)
                    for indices in invalid_indices:
                        for axis in ["row", "column"]:
                            with self.subTest(shape=shape, order=order, dtype=dtype,
                                              indices=repr(indices), axis=axis):
                                rows, columns = (indices, [0]) if axis == "row" else ([0], indices)
                                with self.assertRaisesRegex(IndexError, axis):
                                    pstutil.sub_matrix(matrix, rows, columns, dtype=dtype)
                    for indices in valid_indices:
                        # Independent NumPy indexing supplies the expected Cartesian product.
                        expected_indices = np.asarray(list(indices), dtype=np.intp)
                        for rows, columns in [(indices, [2, 0]), ([2, 0], indices)]:
                            with self.subTest(shape=shape, order=order, dtype=dtype,
                                              rows=repr(rows), columns=repr(columns)):
                                actual = pstutil.sub_matrix(matrix, rows, columns,
                                                            order=order, dtype=dtype)
                                expected_rows = expected_indices if rows is indices else [2, 0]
                                expected_columns = expected_indices if columns is indices else [2, 0]
                                expected = values[np.ix_(expected_rows, expected_columns)]
                                self.assertIsInstance(actual, xp.ndarray)
                                self.assertEqual(actual.dtype, dtype)
                                self.assertTrue(actual.flags[order + "_CONTIGUOUS"])
                                np.testing.assert_array_equal(pstutil.asnumpy(actual), expected)

    def test_sub_matrix_cupy(self):
        # cupy is optional; this test needs it and a GPU.
        try:
            import cupy as cp
        except ModuleNotFoundError:
            raise unittest.SkipTest("cupy is not installed")
        import warnings
        import pysnptools.util as pstutil

        def sub_matrix_and_warnings(val, out_order, out_dtype):
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                result = pstutil.sub_matrix(val, row_index_list, col_index_list, order=out_order, dtype=out_dtype)
            return result, caught

        rng = np.random.default_rng(0)
        row_index_list = [0, 2, 11, 2]
        col_index_list = [6, 5, 4, 3, 2, 1, 0]
        for shape in [(12, 7), (12, 7, 3)]:
            for in_order in ["F", "C"]:
                for in_dtype in [np.float64, np.float32]:
                    matrix = np.asarray(rng.random(shape), dtype=in_dtype, order=in_order)
                    for out_order in ["A", "F", "C"]:
                        for out_dtype in [np.float64, np.float32]:
                            expected, expected_warnings = sub_matrix_and_warnings(matrix, out_order, out_dtype)
                            actual, actual_warnings = sub_matrix_and_warnings(cp.asarray(matrix, order=in_order), out_order, out_dtype)
                            # Both paths report narrowing; GPU float64 output additionally
                            # reports its performance cost. Check every notice rather than
                            # filtering the new category or weakening the precision check.
                            warning_type = GpuFloat64Warning
                            narrowing = in_dtype == np.float64 and out_dtype == np.float32
                            precision_warning = ["Converting float64 to float32 can cause loss of information"]
                            expected_messages = [str(warning.message) for warning in expected_warnings]
                            assert expected_messages == (precision_warning if narrowing else []), expected_messages
                            assert len(actual_warnings) == len(expected_warnings) + (out_dtype == np.float64)
                            for warning in actual_warnings:
                                if str(warning.message) not in expected_messages:
                                    assert out_dtype == np.float64
                                    assert issubclass(warning.category, warning_type)
                                    assert "consumer GPUs" in str(warning.message)
                            assert isinstance(actual, cp.ndarray)
                            assert actual.dtype == expected.dtype
                            assert actual.flags["F_CONTIGUOUS"] == expected.flags["F_CONTIGUOUS"]
                            assert actual.flags["C_CONTIGUOUS"] == expected.flags["C_CONTIGUOUS"]
                            np.testing.assert_array_equal(cp.asnumpy(actual), expected)
        # cupy indexing would silently wrap these, so sub_matrix must reject them.
        matrix = cp.asarray(rng.random((12, 7)))
        for bad_rows, bad_cols in [([12], [0]), ([-1], [0]), (np.array([-1]), [0]), ([0], [7]), ([0], np.array([-7]))]:
            with self.assertRaises(IndexError):
                pstutil.sub_matrix(matrix, bad_rows, bad_cols)

    def test_asarray_warns_when_narrowing_any_float64_array(self):
        # Like sub_matrix, the internal allocator warns for float64-to-float32
        # narrowing of a NumPy array as well as a CuPy array; untyped lists do not warn.
        from pysnptools.util import _array

        with self.assertWarnsRegex(UserWarning, "Converting float64 to float32"):
            result = _array.asarray(np.ones(3), dtype=np.float32, order="C", xp=np)
        self.assertEqual(result.dtype, np.float32)
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            _array.asarray([1.0, 2.0], dtype=np.float32, order="C", xp=np)
            _array.asarray(np.ones(3, dtype=np.float32), dtype=np.float32, order="C", xp=np)

    def test_snpdata_subset_cupy(self):
        # End to end: SnpData created under ARRAY_MODULE=cupy holds a cupy val and must subset.
        # An existing NumPy array needs an explicit xp; xp=None never copies it to the GPU.
        try:
            import cupy as cp
        except ModuleNotFoundError:
            raise unittest.SkipTest("cupy is not installed")
        from unittest.mock import patch
        from pysnptools.snpreader import SnpData
        val = np.arange(24, dtype=np.float64).reshape(6, 4)
        with patch.dict("os.environ", {"ARRAY_MODULE": "cupy"}):
            iid, sid = [["f", str(i)] for i in range(6)], [f"s{j}" for j in range(4)]
            with self.assertRaisesRegex(ValueError, "xp=None selects cupy"):
                SnpData(iid=iid, sid=sid, val=val)
            snpdata = SnpData(iid=iid, sid=sid, val=val, xp="cupy")
            assert isinstance(snpdata.val, cp.ndarray)
            subset = snpdata[np.array([True, False, True, True, False, True]), [3, 1]].read()
        # The subset stays on the GPU with its dtype (cp.asnumpy would also accept numpy).
        assert isinstance(subset.val, cp.ndarray) and subset.val.dtype == np.float64
        np.testing.assert_array_equal(cp.asnumpy(subset.val), val[[0, 2, 3, 5]][:, [3, 1]])

    def test_snpdata_cupy_keeps_dtype(self):
        # A val moved between numpy and cupy (explicit xp; xp=None raises) keeps the
        # dtype that same-module input keeps.
        try:
            import cupy as cp
        except ModuleNotFoundError:
            raise unittest.SkipTest("cupy is not installed")
        from unittest.mock import patch
        from pysnptools.snpreader import SnpData
        iid = [["f", str(i)] for i in range(3)]
        sid = ["s0", "s1"]
        cases = [  # (input dtype, _require_float32_64, expected dtype)
            (np.float32, True, np.float32),
            (np.float64, True, np.float64),
            (np.int8, True, np.float64),
            (np.int8, False, np.int8),
        ]
        for array_module, target in [("cupy", cp), ("numpy", np)]:
            source = np if target is cp else cp
            with patch.dict("os.environ", {"ARRAY_MODULE": array_module}):
                for in_dtype, require, out_dtype in cases:
                    val = source.ones((3, 2), dtype=in_dtype)
                    with self.assertRaisesRegex(ValueError, f"xp=None selects {target.__name__}"):
                        SnpData(iid=iid, sid=sid, val=val, _require_float32_64=require)
                    snpdata = SnpData(iid=iid, sid=sid, val=val, _require_float32_64=require, xp=target)
                    assert isinstance(snpdata.val, target.ndarray), (array_module, in_dtype)
                    assert snpdata.val.dtype == out_dtype, (array_module, in_dtype, require, snpdata.val.dtype)
                    # Same-module input keeps the same dtype.
                    same = SnpData(iid=iid, sid=sid, val=target.ones((3, 2), dtype=in_dtype), _require_float32_64=require)
                    assert same.val.dtype == out_dtype, (array_module, in_dtype, require, same.val.dtype)
                lists = SnpData(iid=iid, sid=sid, val=[[1, 2], [3, 4], [5, 6]])
                assert isinstance(lists.val, target.ndarray) and lists.val.dtype == np.float64


    def test_standardize_cupy_matches_rust(self):
        # cupy arrays are standardized in Python; results and errors must match the Rust path.
        try:
            import cupy as cp
        except ModuleNotFoundError:
            raise unittest.SkipTest("cupy is not installed")
        from pysnptools.standardizer import Unit, Beta
        rng = np.random.default_rng(0)
        val = rng.binomial(2, 0.3, (40, 5)).astype(np.float64)
        val[[1, 5, 9], 0] = np.nan  # some missing
        val[:, 1] = 1.0  # constant
        val[:, 2] = np.nan  # one observed value
        val[0, 2] = 2.0
        for standardizer in [Unit(), Beta(1, 25)]:
            for dtype in [np.float64, np.float32]:
                expected = standardizer.standardize(np.asfortranarray(val.astype(dtype)))
                actual = cp.asnumpy(standardizer.standardize(cp.asarray(val.astype(dtype), order="F")))
                assert actual.dtype == expected.dtype
                np.testing.assert_allclose(actual, expected, rtol=0, atol=8 * np.finfo(dtype).eps)
            all_missing = val.copy()
            all_missing[:, 3] = np.nan
            for array in [np.asfortranarray(all_missing), cp.asarray(all_missing, order="F")]:
                with self.assertRaisesRegex(ValueError, "No individual observed"):
                    standardizer.standardize(array)


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
