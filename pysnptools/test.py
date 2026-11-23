from pysnptools.test_optional_dependencies import _has_hdf5, _requires_hdf5
import numpy as np
import logging
import doctest
import shutil
import cloudpickle as pickle

from pysnptools.snpreader import Bed
from pysnptools.snpreader import SnpHdf5, SnpNpz
from pysnptools.snpreader import Dat
from pysnptools.snpreader import Dense
from pysnptools.snpreader import Pheno
from pysnptools.snpreader import Ped
from pysnptools.snpreader import DistributedBed
from pysnptools.standardizer import Unit
from pysnptools.standardizer import Beta
from pysnptools.util import create_directory_if_necessary
from pysnptools.kernelreader.test import TestKernelReader
from pysnptools.kernelreader.test import TestKrDocStrings
from pysnptools.test_optional_dependencies import TestOptionalHdf5, TestOptionalBgen
from pysnptools.pstreader.test import TestPstReader
from pysnptools.pstreader.test import TestPstDocStrings
from pysnptools.pstreader.pstmemmap import TestPstMemMap
from pysnptools.snpreader.snpmemmap import TestSnpMemMap
from pysnptools.snpreader.snpgen import TestSnpGen
from pysnptools.snpreader.distributedbed import TestDistributedBed
from pysnptools.util.generate import TestGenerate
from pysnptools.util._example_file import TestExampleFile
from pysnptools.kernelreader.test import _fortesting_JustCheckExists
from pysnptools.util.intrangeset import TestIntRangeSet
from pysnptools.util.test import TestUtilTools
from pysnptools.util.filecache.test import TestFileCache
from pysnptools.util.filecache.test_hashdown_retry import TestHashdownRetry
from pysnptools.util.mapreduce1.test_threads import TestRunnerThreads
from pysnptools.distreader.test import TestDistReaders

import unittest
import os.path
import time


class TestBetaPython(unittest.TestCase):
    def test_empty_data_factory(self):
        from pysnptools.snpreader import SnpData
        from pysnptools.kernelreader import KernelData

        source = SnpData(iid=[], sid=["s0", "s1"], val=None, xp=np)
        kernel = KernelData(iid=[], val=None, xp=np)
        for target, shape in ((source, (0, 2)), (kernel, (0, 0))):
            self.assertIsInstance(target.val, np.ndarray)
            self.assertEqual(target.val.dtype, np.float64)
            self.assertEqual(target.val.shape, shape)

    def test_trained_statistics_across_dtypes(self):
        from pysnptools.snpreader import SnpData

        iid = [["family", str(sample)] for sample in range(3)]
        sid = ["variable", "dosage", "constant"]
        training_values = [[0, 0.25, 1], [1, 1, 1], [2, 1.75, np.nan]]
        test_values = [[0.5, np.nan, 0], [1.75, 0.25, 2], [np.nan, 1.25, 1]]
        for training_dtype, target_dtype in ((np.float32, np.float64), (np.float64, np.float32)):
            for order in ("C", "F"):
                for standardizer in (Unit(), Beta(2, 3)):
                    training = SnpData(iid=iid, sid=sid, xp=np,
                                       val=np.array(training_values, dtype=training_dtype, order=order))
                    _, trained = training.standardize(standardizer, return_trained=True)
                    original_stats = trained.stats.copy()
                    # Stored moments are rounded to the target's supported dtype.
                    means, deviations = trained.stats.astype(target_dtype).astype(np.float64).T
                    if isinstance(standardizer, Unit):
                        factors = 1 / deviations
                    else:
                        frequencies = np.minimum(means / 2, 1 - means / 2)
                        factors = 12 * frequencies * (1 - frequencies) ** 2
                        factors[np.isinf(deviations)] = 0
                    expected = (np.array(test_values) - means) * factors
                    expected[np.isnan(test_values)] = 0
                    for force_python_only in (False, True):
                        with self.subTest(training_dtype=training_dtype, target_dtype=target_dtype,
                                          order=order, standardizer=standardizer,
                                          force_python_only=force_python_only):
                            target = SnpData(iid=iid, sid=sid, xp=np,
                                             val=np.array(test_values, dtype=target_dtype, order=order))
                            self.assertIs(target.standardize(trained, force_python_only=force_python_only), target)
                            self.assertEqual(target.val.dtype, target_dtype)
                            # Small bounded inputs; allow arithmetic and density roundoff.
                            tolerance = 32 * np.finfo(target_dtype).eps
                            np.testing.assert_allclose(target.val, expected, rtol=tolerance, atol=tolerance)
                            np.testing.assert_array_equal(trained.stats, original_stats)
                            self.assertEqual(trained.stats.dtype, training_dtype)

    def test_constant_fractional_snps(self):
        # bed-reader's native CPU standardizer (force_python_only=False) is not
        # covered: it leaves small nonzero residues for these inputs, tracked
        # for bed-reader rather than fixed here.
        self._check_constant_fractional_snps(np, force_python_only=True)

    def test_constant_fractional_snps_cupy(self):
        try:
            import cupy
        except ModuleNotFoundError:
            self.skipTest("requires CuPy and a working CUDA GPU")
        self._check_constant_fractional_snps(cupy, force_python_only=False)

    def _check_constant_fractional_snps(self, xp, force_python_only):
        # A constant dosage such as 0.1 has a rounded mean (0.1 averages to
        # 0.09999999999999999), so its computed variance can be tiny but
        # nonzero. It must still be treated as constant: standardized to zero,
        # with an infinite trained deviation, not to +/-1.
        from pysnptools.snpreader import SnpData
        from pysnptools.util import asnumpy

        # Each case standardized a constant column to +/-1 before the fix.
        for dtype, value, count in ((np.float64, 0.1, 3), (np.float64, 0.3, 7), (np.float64, 1.7, 7),
                                    (np.float32, 0.1, 7), (np.float32, 1.7, 3)):
            values = np.empty((count + 1, 2), dtype=dtype)
            values[:, 0] = value
            values[-1, 0] = np.nan  # missing values do not break constancy
            values[:, 1] = np.arange(count + 1) % 3
            for order in ("C", "F"):
                for standardizer in (Unit(), Beta(1, 25)):
                    with self.subTest(dtype=dtype, value=value, count=count, order=order,
                                      standardizer=standardizer):
                        data = SnpData(iid=[["f", str(row)] for row in range(count + 1)], sid=["constant", "varying"],
                                       val=xp.array(values, order=order), xp=xp)
                        _, trained = data.standardize(standardizer, return_trained=True,
                                                      force_python_only=force_python_only)
                        np.testing.assert_array_equal(asnumpy(data.val)[:, 0], 0)
                        self.assertTrue(np.isinf(asnumpy(trained.stats)[0, 1]))
                        self.assertTrue(np.isfinite(asnumpy(trained.stats)[1, 1]))

    def test_blocked_trained_kernel(self):
        self._check_blocked_trained_kernel(np)

    def test_blocked_trained_kernel_cupy(self):
        try:
            import cupy
        except ModuleNotFoundError:
            self.skipTest("requires CuPy and a working CUDA GPU")
        self._check_blocked_trained_kernel(cupy)

    def _check_blocked_trained_kernel(self, xp):
        from scipy.stats import beta as beta_distribution
        from pysnptools.kernelreader import SnpKernel
        from pysnptools.snpreader import SnpData
        from pysnptools.standardizer import BetaTrained
        from pysnptools.util import asnumpy

        training_values = np.array([
            [0, 2, 1, 0, 0.25, 0, 2],
            [0, 2, 1, 0, 0.5, 1, 2],
            [0, 2, 1, 1, 1.75, 1, 1],
            [np.nan, 2, 1, 1, 1.5, 2, np.nan],
            [0, np.nan, 1, 2, 0.75, np.nan, 2],
            [0, 2, 1, np.nan, 1, 1, 2],
        ])
        # Test means differ from training means; constant training SNPs vary here.
        values = np.array([
            [0, 1, 2, 2, 1.75, 2, 0],
            [1, 2, 0, 1, np.nan, 2, 1],
            [2, 0, 1, 2, 1.5, 0, 0],
            [np.nan, 1, 0, np.nan, 1, 1, np.nan],
        ])
        # Repeated IDs have distinct moments, so positional correspondence matters.
        sid = ["constant0", "constant2", "constant1", "repeated", "dosage", "other", "repeated"]
        iid = [["family", str(sample)] for sample in range(len(values))]
        for dtype in (np.float32, np.float64):
            for order in ("C", "F"):
                for beta_a, beta_b in ((0.5, 0.5), (1, 25)):
                    training = SnpData(iid=[["train", str(sample)] for sample in range(len(training_values))],
                                       sid=sid, val=np.array(training_values, dtype=dtype, order=order), xp=np)
                    _, trained = training.standardize(Beta(beta_a, beta_b), return_trained=True)
                    original_stats = trained.stats.copy()
                    standardized = np.zeros(values.shape, dtype=dtype)
                    for column in range(3, len(sid)):
                        mean = dtype(np.nanmean(training_values[:, column]))
                        minor_frequency = min(float(mean) / 2, 1 - float(mean) / 2)
                        weight = dtype(beta_distribution.pdf(minor_frequency, beta_a, beta_b))
                        observed = ~np.isnan(values[:, column])
                        standardized[observed, column] = (values[observed, column].astype(dtype) - mean) * weight
                    expected = standardized @ standardized.T
                    reader = SnpData(iid=iid, sid=sid, val=xp.array(values, dtype=dtype, order=order), xp=xp)
                    # Seven SNP products and block additions contribute rounding.
                    tolerance = 64 * np.finfo(dtype).eps
                    unblocked = reader.read_kernel(trained, dtype=dtype, order=order, xp=xp).val
                    # The native statrs Beta(1,25) density differs from SciPy's
                    # by about 1.5e-13 relatively here, even without blocking.
                    # Keep the block comparison at floating-point roundoff;
                    # allow that library difference only in the independent check.
                    reference_tolerance = max(tolerance, 1e-12) if xp is np and beta_a == 1 else tolerance
                    np.testing.assert_allclose(asnumpy(unblocked), expected,
                                               rtol=reference_tolerance, atol=reference_tolerance * len(sid))
                    for block_size in (1, 3, 9):
                        with self.subTest(backend=xp.__name__, dtype=dtype, order=order,
                                          beta_a=beta_a, beta_b=beta_b, block_size=block_size):
                            kernel = reader.read_kernel(trained, block_size=block_size, dtype=dtype, order=order, xp=xp)
                            self.assertIsInstance(kernel.val, xp.ndarray)
                            self.assertEqual(kernel.val.dtype, dtype)
                            self.assertTrue(kernel.val.flags["C_CONTIGUOUS" if order == "C" else "F_CONTIGUOUS"])
                            np.testing.assert_allclose(asnumpy(kernel.val), asnumpy(unblocked), rtol=tolerance, atol=tolerance * len(sid))
                            specified = SnpKernel(reader, trained, block_size=block_size).read(dtype=dtype, order=order, xp=xp)
                            np.testing.assert_allclose(asnumpy(specified.val), asnumpy(unblocked), rtol=tolerance, atol=tolerance * len(sid))
                            result, returned = reader._read_kernel(trained, block_size=block_size, dtype=dtype,
                                                                  order=order, xp=xp, return_trained=True)
                            self.assertIs(returned, trained)
                            np.testing.assert_allclose(asnumpy(result), asnumpy(unblocked), rtol=tolerance, atol=tolerance * len(sid))
                    np.testing.assert_array_equal(trained.stats, original_stats)
                    np.testing.assert_array_equal(trained.sid, sid)
                    # Blocking must not turn an invalid full-cohort alignment into
                    # an accepted subset, permutation, or unknown SNP selection.
                    for selection in (reader[:, :-1], reader[:, ::-1],
                                      SnpData(iid=iid, sid=["unknown"] + sid[1:], val=values, xp=np)):
                        with self.assertRaisesRegex(AssertionError, "sid in training and use"):
                            selection.read_kernel(trained, block_size=3, dtype=dtype, xp=xp)
                    wrong_stats = BetaTrained(beta_a, beta_b, trained.sid,
                                             np.concatenate([trained.stats, trained.stats[:1]]))
                    with self.assertRaisesRegex(AssertionError, "stats must have size"):
                        reader.read_kernel(wrong_stats, block_size=3, dtype=dtype, xp=xp)

    def test_trained_checks_survive_optimize(self):
        # Under python -O, plain assert statements vanish. Run each user-input
        # check in a subprocess with -O and require its AssertionError.
        import subprocess
        import sys

        script = r"""
import numpy as np
from pysnptools.snpreader import SnpData
from pysnptools.standardizer import Beta, BetaTrained
assert False, "expected -O to strip this assert"
"""
        cases = {
            "sid in training and use": """
values = np.array([[0.0, 1.0], [1.0, 2.0], [2.0, 0.0]])
data = SnpData(iid=[["f", str(i)] for i in range(3)], sid=["s0", "s1"], val=values)
_, trained = data.read().standardize(Beta(1, 25), return_trained=True)
data[:, ::-1].read().standardize(trained)
""",
            "stats must have size": """
values = np.array([[0.0, 1.0], [1.0, 2.0], [2.0, 0.0]])
data = SnpData(iid=[["f", str(i)] for i in range(3)], sid=["s0", "s1"], val=values)
BetaTrained(1, 25, data.sid, np.ones((3, 2)))._standardize_unit_and_beta(
    data.val, is_beta=True, a=1, b=25, apply_in_place=True, use_stats=True,
    stats=np.ones((3, 2)), num_threads=1)
""",
            "blocked sid": """
values = np.array([[0.0, 1.0], [1.0, 2.0], [2.0, 0.0]])
data = SnpData(iid=[["f", str(i)] for i in range(3)], sid=["s0", "s1"], val=values)
_, trained = data.read().standardize(Beta(1, 25), return_trained=True)
data[:, ::-1].read_kernel(trained, block_size=1)
""",
            "blocked stats": """
values = np.array([[0.0, 1.0], [1.0, 2.0], [2.0, 0.0]])
data = SnpData(iid=[["f", str(i)] for i in range(3)], sid=["s0", "s1"], val=values)
data.read_kernel(BetaTrained(1, 25, data.sid, np.ones((3, 2))), block_size=1)
""",
            "same 'a'": """
Beta(1, 25)._merge_trained([BetaTrained(1, 25, ["s0"], np.ones((1, 2))),
                            BetaTrained(2, 25, ["s1"], np.ones((1, 2)))])
""",
            "same 'b'": """
Beta(1, 25)._merge_trained([BetaTrained(1, 25, ["s0"], np.ones((1, 2))),
                            BetaTrained(1, 26, ["s1"], np.ones((1, 2)))])
""",
        }
        messages = {"blocked sid": "sid in training and use", "blocked stats": "stats must have size"}
        for name, body in cases.items():
            with self.subTest(check=name):
                completed = subprocess.run([sys.executable, "-O", "-c", script + body],
                                           capture_output=True, text=True, check=False)
                self.assertNotEqual(completed.returncode, 0, completed.stdout)
                last_line = completed.stderr.strip().splitlines()[-1]
                self.assertTrue(last_line.startswith("AssertionError: "), completed.stderr)
                self.assertIn(messages.get(name, name), last_line)

    def test_constant_columns(self):
        from scipy.stats import beta as beta_distribution
        from pysnptools.snpreader import SnpData

        values = np.array([
            [0, 2, 1, 0, 0.25],
            [0, 2, 1, 1, 0.5],
            [np.nan, 2, 1, 2, 1.75],
            [0, np.nan, 1, np.nan, 1.5],
        ])
        iid = [["family", str(sample)] for sample in range(len(values))]
        sid = ["constant0", "constant2", "constant1", "genotype", "dosage"]
        for dtype in (np.float32, np.float64):
            for order in ("C", "F"):
                for beta_a, beta_b in ((0.5, 0.5), (0.5, 2), (2, 0.5), (1, 25)):
                    with self.subTest(dtype=dtype, order=order, beta_a=beta_a, beta_b=beta_b):
                        expected = np.zeros(values.shape, dtype=dtype, order=order)
                        for column in (3, 4):
                            observed = ~np.isnan(values[:, column])
                            mean = values[observed, column].mean()
                            minor_frequency = min(mean / 2, 1 - mean / 2)
                            weight = beta_distribution.pdf(minor_frequency, beta_a, beta_b)
                            expected[observed, column] = (values[observed, column] - mean) * weight
                        # These well-scaled columns need only a few rounded operations.
                        tolerance = 32 * np.finfo(dtype).eps
                        standardizer = Beta(beta_a, beta_b)
                        for python_only in (False, True):
                            data = SnpData(iid=iid, sid=sid, val=np.array(values, dtype=dtype, order=order), xp=np)
                            with np.errstate(invalid="raise"):
                                result, trained = data.standardize(
                                    standardizer, return_trained=True, force_python_only=python_only
                                )
                            self.assertIs(result, data)
                            self.assertEqual(data.val.dtype, dtype)
                            self.assertTrue(data.val.flags["C_CONTIGUOUS" if order == "C" else "F_CONTIGUOUS"])
                            np.testing.assert_array_equal(data.val[:, :3], 0)
                            np.testing.assert_allclose(data.val, expected, rtol=tolerance, atol=tolerance)
                            np.testing.assert_array_equal(trained.stats[:3, 0], [0, 2, 1])
                            self.assertTrue(np.isinf(trained.stats[:3, 1]).all())

                            changed = np.array(values, dtype=dtype, order=order)
                            changed[:, :3] = [[0, 1, 2], [1, 2, 0], [2, 0, 1], [np.nan, 1, 0]]
                            fresh = SnpData(iid=iid, sid=sid, val=changed, xp=np)
                            with np.errstate(invalid="raise"):
                                fresh.standardize(trained, force_python_only=True)
                            np.testing.assert_allclose(fresh.val, expected, rtol=tolerance, atol=tolerance)

                        # A contiguous view of a larger allocation automatically uses
                        # the Python path even when force_python_only is False.
                        shape = (5, 5) if order == "C" else (4, 6)
                        backing = np.empty(shape, dtype=dtype, order=order)
                        view = backing[:4, :5]
                        view[:] = values
                        self.assertTrue(view.flags["C_CONTIGUOUS" if order == "C" else "F_CONTIGUOUS"])
                        fallback = SnpData(iid=iid, sid=sid, val=view, xp=np)
                        with np.errstate(invalid="raise"):
                            fallback.standardize(standardizer)
                        np.testing.assert_allclose(fallback.val, expected, rtol=tolerance, atol=tolerance)


class TestDiagKtoN(unittest.TestCase):
    def test_diagKtoN_invalid_denominator(self):
        self._check_diagKtoN_denominator(np)

    def test_diagKtoN_invalid_denominator_cupy(self):
        try:
            import cupy
        except ModuleNotFoundError:
            self.skipTest("requires CuPy and a working CUDA GPU")
        self._check_diagKtoN_denominator(cupy)

    def _check_diagKtoN_denominator(self, xp):
        from pysnptools.kernelreader import KernelData
        from pysnptools.snpreader import SnpData
        from pysnptools.standardizer import DiagKtoN
        from pysnptools.util import asnumpy

        def data(values, dtype, order, kernel):
            iid = [["family", str(sample)] for sample in range(values.shape[0])]
            val = xp.array(values, dtype=dtype, order=order)
            if kernel:
                return KernelData(iid=iid, val=val, xp=xp)
            return SnpData(iid=iid, sid=[str(marker) for marker in range(values.shape[1])], val=val, xp=xp)

        for dtype in (np.float32, np.float64):
            for order in ("C", "F"):
                for kernel in (False, True):
                    invalid = [np.zeros((2, 2)), np.array([[np.nan, 0], [0, 1]]),
                               np.array([[np.inf, 0], [0, 1]]), np.empty((0, 0))]
                    if kernel:
                        invalid.append(-np.eye(2))
                    else:
                        invalid.extend([np.empty((0, 2)), np.empty((2, 0))])
                    for values in invalid:
                        for return_trained in (False, True):
                            with self.subTest(backend=xp.__name__, dtype=dtype, order=order,
                                              kernel=kernel, values=values, return_trained=return_trained):
                                target = data(values, dtype, order, kernel)
                                with self.assertRaisesRegex(ValueError, "finite, positive normalization denominator"):
                                    target.standardize(DiagKtoN(), return_trained=return_trained)
                                np.testing.assert_array_equal(asnumpy(target.val), values)

                    # Independently known scales: kernel trace 8 -> factor 1/4;
                    # SNP sum of squares 4 -> factor 1/2, applied as sqrt(1/2).
                    values = np.array([[4., 2.], [2., 4.]]) if kernel else np.array([[-1., 1.], [1., -1.]])
                    factor = 0.25 if kernel else 0.5
                    expected = values * (factor if kernel else np.sqrt(factor))
                    target = data(values, dtype, order, kernel)
                    result, trained = target.standardize(DiagKtoN(), return_trained=True)
                    self.assertIs(result, target)
                    self.assertEqual(trained.factor, factor)
                    self.assertEqual(target.val.dtype, dtype)
                    self.assertTrue(target.val.flags["C_CONTIGUOUS" if order == "C" else "F_CONTIGUOUS"])
                    # One scalar multiplication; allow sqrt and multiplication roundoff.
                    np.testing.assert_allclose(asnumpy(target.val), expected,
                                               rtol=8 * np.finfo(dtype).eps, atol=0)
                    # Applying a fitted scale to zero test data does not refit it.
                    zero = data(np.zeros((2, 2)), dtype, order, kernel)
                    zero.standardize(trained)
                    np.testing.assert_array_equal(asnumpy(zero.val), 0)


class TestPySnpTools(unittest.TestCase):
    def xtest_aaa_hdf5_speed(self):  # !!too slow to use all the time
        # currentFolder + "/examples/toydata.5chrom.bed"
        # currentFolder + "/examples/delme.hdf5"
        bedFileName = r"c:\Source\carlk\fastlmm\tests\datasets\all_chr.maf0.001.N300"  # !! local paths
        hdf5Pattern = r"c:\Source\carlk\fastlmm\tests\datasets\del.{0}.hdf5"  # !!
        tt0 = time.time()
        snpreader_bed = Bed(bedFileName, count_A1=False)

        S0 = snpreader_bed.sid_count
        snp_index_list0 = range(min(S0, 15000))

        hdf5FileName = hdf5Pattern.format(len(snp_index_list0))

        snpDataBed = snpreader_bed[:, snp_index_list0].read()
        tt1 = time.time()
        logging.info("Read bed %.2f seconds" % (tt1 - tt0))
        SnpHdf5.write(hdf5FileName, snpDataBed)
        tt2 = time.time()
        logging.info("write SnpHdf5 bed %.2f seconds" % (tt2 - tt1))

        snpreader_hdf5 = SnpHdf5(hdf5FileName)
        assert snpreader_hdf5.iid_count == snpreader_bed.iid_count
        S = snpreader_hdf5.sid_count
        N_original = snpreader_hdf5.iid_count
        iid_index_list = sorted(range(N_original - 1, 0, -2))

        snp_index_list = sorted(range(S - 1, 0, -2))  # !!
        # snp_index_list = range(S//2)

        snpreader_hdf5 = snpreader_hdf5[iid_index_list, :]
        _ = snpreader_hdf5[:, snp_index_list].read()
        tt3 = time.time()
        logging.info(
            "read SnpHdf5 with reversed indexes bed %.2f seconds" % (tt3 - tt2)
        )

        _ = snpreader_hdf5[:, snp_index_list].read(order="C")
        tt4 = time.time()
        logging.info(
            "read SnpHdf5 C with reversed indexes bed %.2f seconds" % (tt4 - tt3)
        )

        print("the end")

    @classmethod
    def setUpClass(self):
        self.currentFolder = os.path.dirname(os.path.realpath(__file__))
        # TODO: get data set with NANs!
        snpreader = Bed(
            self.currentFolder + "/examples/toydata.5chrom.bed", count_A1=False
        )
        assert snpreader.pos[0, 1] != snpreader.pos[0, 1]  # real assert
        self.pheno_fn = self.currentFolder + "/examples/toydata.phe"
        self.snpdata = snpreader.read(order="F", force_python_only=True)
        self.snps = self.snpdata.val

    # def test_bug_bed_reader_27(self):
    #     from pysnptools.snpreader import Bed, SnpData
    #     import pysnptools.util as pstutil
    #     from pysnptools.util import example_file

    #     bed_fn = example_file("pysnptools/examples/toydata.*", "*.bed")
    #     bed = Bed(bed_fn, count_A1=False)
    #     bed.pos[:, :] = np.nan

    #     file_name = "tempdir/bug_bed_reader_27.bed"
    #     pstutil.create_directory_if_necessary(file_name)
    #     Bed.write(file_name, bed.read(), count_A1=False)
    #     bed2 = Bed(file_name, count_A1=False).read()
    #     assert np.all(np.isnan(bed2.pos))

    def test_write_bad_value_and_good(self):
        from pysnptools.snpreader import Pheno, Bed, SnpData
        import pysnptools.util as pstutil
        from pysnptools.util import example_file  # Download and return local file name

        pheno_fn = example_file("pysnptools/examples/toydata.phe")
        snpdata = Pheno(pheno_fn).read()  # Read data from Pheno format
        pstutil.create_directory_if_necessary("tempdir/toydata.5chrom.bed")
        error_seen = False
        try:
            Bed.write(
                "tempdir/toydata.5chrom.bed", snpdata, count_A1=False
            )  # Write data in Bed format
        except Exception:
            error_seen = True
        assert error_seen
        # Can write from an int8 array, too.
        snpdata_int = SnpData(
            val=np.int_(snpdata.val).astype("int8"),
            iid=snpdata.iid,
            sid=snpdata.sid,
            pos=snpdata.pos,
            _require_float32_64=False,
        )
        assert snpdata_int.val.dtype == "int8"
        error_seen = False
        try:
            Bed.write(
                "tempdir/toydata.5chrom.bed",
                snpdata_int,
                count_A1=False,
                _require_float32_64=False,
            )
        except Exception:
            error_seen = True
        assert error_seen

        bed_fn = example_file("pysnptools/examples/toydata.5chrom.*", "*.bed")
        snpdata = Bed(bed_fn)[:, ::2].read()  # Read every-other SNP
        pstutil.create_directory_if_necessary("tempdir/everyother.bed")
        Bed.write(
            "tempdir/everyother.bed", snpdata, count_A1=False
        )  # Write data in Bed format
        # Bed('tempdir/everyother.bed',count_A1=False)
        # Can write from an int8 array, too.
        snpdata_int = SnpData(
            val=np.int_(snpdata.val).astype("int8"),
            iid=snpdata.iid,
            sid=snpdata.sid,
            pos=snpdata.pos,
            _require_float32_64=False,
        )
        snpdata_int.val.dtype
        assert snpdata_int.val.dtype == "int8"
        Bed.write(
            "tempdir/everyother.bed",
            snpdata_int,
            count_A1=False,
            _require_float32_64=False,
        )
        # Bed('tempdir/everyother',count_A1=False)

    def test_val_assign(self):
        from pysnptools.snpreader import SnpData
        from pysnptools.kernelreader import KernelData

        iid_count = 3
        sid_count = 2
        iid = [[str(i), str(i)] for i in range(iid_count)]
        sid = [str(s) for s in range(sid_count)]
        np.random.seed(0)
        vali = np.random.randint(10, size=[iid_count, sid_count])
        vali2 = np.random.randint(10, size=[iid_count, sid_count])
        snpdata = SnpData(iid=iid, sid=sid, val=vali)
        snpdata.val = vali2
        assert snpdata.val.dtype == np.float64

        valk = np.random.randint(10, size=[iid_count, iid_count])
        valk2 = np.random.randint(10, size=[iid_count, iid_count])
        kerneldata = KernelData(iid=iid, val=valk)
        assert np.float64 == kerneldata.val.dtype
        kerneldata.val = valk2
        assert np.float64 == kerneldata.val.dtype

        valstr = np.array(vali, dtype=np.str_)
        valstr[0, 0] = "cannot convert"
        error_seen = False
        try:
            snpdata.val = valstr
        except Exception:
            error_seen = True
        assert error_seen
        error_seen = False
        try:
            kerneldata.val = vali
        except Exception:
            error_seen = True
        assert error_seen
        error_seen = False

    def test_diagKtoN(self):
        """
        make sure standardization on SNPs results in sum(diag(K))=N
        """

        np.random.seed(42)
        m = np.random.random((100, 1000))
        from pysnptools.standardizer import DiagKtoN

        s = DiagKtoN()
        s.standardize(m)
        K = m.dot(m.T)
        sum_diag = np.sum(np.diag(K))

        np.testing.assert_almost_equal(100, sum_diag)

    def test_c_reader_bed(self):
        snpreader = Bed(
            self.currentFolder + "/examples/toydata.5chrom.bed", count_A1=False
        )
        self.c_reader(snpreader)

    def test_c_reader_bed_count_A1(self):
        snpreader = Bed(
            self.currentFolder + "/examples/toydata.5chrom.bed", count_A1=True
        )
        snpdata = snpreader.read()
        snpdata.val = 2 - snpdata.val
        self.c_reader(snpdata)

    def test_p_reader_bed(self):
        snpreader = Bed(
            self.currentFolder + "/examples/toydata.5chrom.bed", count_A1=False
        ).read(force_python_only=True)
        self.c_reader(snpreader)

    def test_p_reader_bed_count_A1(self):
        snpreader = Bed(
            self.currentFolder + "/examples/toydata.5chrom.bed", count_A1=True
        )
        snpdata = snpreader.read(force_python_only=True)
        snpdata.val = 2 - snpdata.val
        self.c_reader(snpdata)

    def test_bed_2021(self):
        snpreader = Bed(
            self.currentFolder + "/examples/toydata.5chrom.bed",
            count_A1=True,
            fam_filename=self.currentFolder + "/examples/toydata.5chrom.2.fam",
            bim_filename=self.currentFolder + "/examples/toydata.5chrom.XYMT.bim",
        )
        snpdata = snpreader.read()

        from pysnptools.snpreader import reverse_plink_chrom_map
        import pysnptools.util as pstutil

        pstutil.create_directory_if_necessary("tempdir/test_reverse_chrom_map.bed")
        Bed.write(
            "tempdir/test_reverse_chrom_map.bed",
            snpdata,
            count_A1=False,
            reverse_chrom_map=reverse_plink_chrom_map,
        )

        see_exception = False
        try:
            snpreader = Bed(
                self.currentFolder + "/examples/toydata.5chrom.bed",
                count_A1=True,
                bim_filename=self.currentFolder + "/examples/toydata.5chrom.bad.bim",
            )
            snpdata = snpreader.read()
        except ValueError:
            see_exception = True
        assert see_exception

        snpreader = Bed(
            self.currentFolder + "/examples/toydata.5chrom.bed",
            count_A1=True,
            bim_filename=self.currentFolder + "/examples/toydata.5chrom.bad.bim",
            chrom_map={"B": 1, "A": 2, "D": 3},
        )
        snpdata = snpreader.read()

    def test_bed_int8(self):
        snpreader = Bed(
            self.currentFolder + "/../tests/datasets/distributed_bed_test1_X.bed",
            count_A1=True,
        )
        ref = snpreader.read()
        ref.val[ref.val != ref.val] = -127
        ref.val = ref.val.astype("int8")
        for force_python_only in [False, True]:
            for order in ["F", "C"]:
                snpdata = snpreader.read(
                    dtype="int8",
                    _require_float32_64=False,
                    force_python_only=force_python_only,
                    order=order,
                )
                assert snpdata.val.dtype == "int8"
                assert np.allclose(snpdata.val, ref.val, equal_nan=True)
                output = "tempdir/snpreader/int8.bed"
                create_directory_if_necessary(output)
                for count_A1 in [False, True]:
                    bed2 = Bed.write(
                        output,
                        snpdata,
                        count_A1=count_A1,
                        _require_float32_64=False,
                        force_python_only=force_python_only,
                    )
                    assert np.allclose(
                        bed2.read(
                            dtype="int8",
                            _require_float32_64=False,
                            force_python_only=force_python_only,
                        ).val,
                        ref.val,
                        equal_nan=True,
                    )

    def test_scalar_index(self):
        snpreader = Bed(
            self.currentFolder + "/examples/toydata.5chrom.bed", count_A1=False
        )
        arr = np.int64(1)
        snpreader[arr, arr]

    @_requires_hdf5
    def test_c_reader_hdf5(self):
        snpreader = SnpHdf5(self.currentFolder + "/examples/toydata.snpmajor.snp.hdf5")
        self.c_reader(snpreader)

    @_requires_hdf5
    def test_hdf5_case3(self):
        snpreader1 = SnpHdf5(
            self.currentFolder + "/examples/toydata.snpmajor.snp.hdf5"
        )[::2, :]
        snpreader2 = Bed(
            self.currentFolder + "/examples/toydata.5chrom.bed", count_A1=False
        )[::2, :]
        self.assertTrue(
            np.allclose(
                snpreader1.read().val, snpreader2.read().val, rtol=1e-05, atol=1e-05
            )
        )

    @_requires_hdf5
    def test_standardize_hdf5(self):
        snpreader = SnpHdf5(self.currentFolder + "/examples/toydata.iidmajor.snp.hdf5")
        self.standardize(snpreader)

    def test_c_reader_npz(self):
        snpreader = SnpNpz(self.currentFolder + "/examples/toydata10.snp.npz")
        snpdata = snpreader.read(order="F", force_python_only=False)
        snp_c = snpdata.val

        self.assertEqual(np.float64, snp_c.dtype)
        self.assertTrue(np.allclose(self.snps[:, :10], snp_c, rtol=1e-05, atol=1e-05))

        snpreader1 = SnpNpz(self.currentFolder + "/examples/toydata10.snp.npz")
        snpreader2 = Bed(
            self.currentFolder + "/examples/toydata.5chrom.bed", count_A1=False
        )[:, :10]
        self.assertTrue(
            np.allclose(
                snpreader1.read().val, snpreader2.read().val, rtol=1e-05, atol=1e-05
            )
        )

        snpdata.val[1, 2] = (
            np.nan
        )  # Inject a missing value to test writing and reading missing values
        output = "tempdir/snpreader/toydata10.snp.npz"
        create_directory_if_necessary(output)
        SnpNpz.write(output, snpdata)
        snpdata2 = SnpNpz(output).read()
        np.testing.assert_array_almost_equal(snpdata.val, snpdata2.val, decimal=10)

        snpdata3 = snpdata[:, 0:0].read()  # create snpdata with no sids
        output = "tempdir/snpreader/toydata0.snp.npz"
        SnpNpz.write(output, snpdata3)
        snpdata4 = SnpNpz(output).read()
        assert snpdata3 == snpdata4

    def test_standardize_npz(self):
        snpreader = SnpNpz(self.currentFolder + "/examples/toydata10.snp.npz")
        self.standardize(snpreader)

    def test_c_reader_dat(self):
        snpreader = Dat(self.currentFolder + "/examples/toydata.dat")[:, ::100]
        _fortesting_JustCheckExists().input(snpreader)

        snpdata1 = snpreader.read()
        self.assertEqual(np.float64, snpdata1.val.dtype)
        self.assertTrue(
            np.allclose(self.snps[:, ::100], snpdata1.val, rtol=1e-05, atol=1e-05)
        )

        snpdata1.val[1, 2] = (
            np.nan
        )  # Inject a missing value to test writing and reading missing values
        output = "tempdir/snpreader/toydata.dat"
        create_directory_if_necessary(output)
        Dat.write(output, snpdata1)
        snpdata2 = Dat(output).read()
        np.testing.assert_array_almost_equal(snpdata1.val, snpdata2.val, decimal=10)

        snpdata3 = snpdata1[:, 0:0].read()  # create snpdata with no sids
        output = "tempdir/snpreader/toydata3.dat"
        Dat.write(output, snpdata3)
        snpdata4 = Dat(output).read()
        assert snpdata3 == snpdata4

    @staticmethod
    def assert_match_012_210(snpdata1, snpdata2):
        for sid_index in range(
            snpdata1.sid_count
        ):  # Check that every row matches (except OK if 0,1,2 can be 2,1,0)
            goods1 = (
                snpdata1.val[:, sid_index] == snpdata1.val[:, sid_index]
            )  # find non-missing
            goods2 = (
                snpdata2.val[:, sid_index] == snpdata2.val[:, sid_index]
            )  # find non-missing
            assert (goods1 == goods2).all()  # assert that they agree on non-missing
            is_ok = (
                snpdata1.val[goods1, sid_index] == snpdata2.val[goods2, sid_index]
            ).all() or (
                snpdata1.val[goods1, sid_index]
                == snpdata2.val[goods2, sid_index] * -1 + 2
            ).all()
            assert is_ok

    def test_c_reader_ped(self):
        if False:  # Too slow for routine testing
            snpdata1 = Ped(self.currentFolder + "/examples/toydata.ped")[
                ::25, ::1000
            ].read()
            self.assertEqual(np.float64, snpdata1.val.dtype)
            TestPySnpTools.assert_match_012_210(
                self.snpdata[::25, ::1000].read(), snpdata1
            )
        else:
            snpdata1 = self.snpdata[::25, ::1000].read()

        output = "tempdir/snpreader/toydata.ped"
        create_directory_if_necessary(output)

        snpdata1.val[1, 2] = (
            np.nan
        )  # Inject a missing value to test writing and reading missing values
        Ped.write(output, snpdata1)
        snpreader = Ped(output)
        _fortesting_JustCheckExists().input(snpreader)
        _ = str(snpreader)
        snpdata2 = snpreader.read()
        TestPySnpTools.assert_match_012_210(snpdata1, snpdata2)

    def test_c_reader_pheno(self):
        snpdata1 = Pheno(self.currentFolder + "/examples/toydata.phe").read()

        self.assertEqual(np.float64, snpdata1.val.dtype)

        snpdata1.val[1, 0] = (
            np.nan
        )  # Inject a missing value to test writing and reading missing values
        output = "tempdir/snpreader/toydata.phe"
        create_directory_if_necessary(output)
        Pheno.write(output, snpdata1)
        snpreader = Pheno(output)
        _fortesting_JustCheckExists().input(snpreader)
        _ = str(snpreader)
        snpdata2 = snpreader.read()
        np.testing.assert_array_almost_equal(snpdata1.val, snpdata2.val, decimal=10)

        snpdata1 = Pheno(self.currentFolder + "/examples/toydata.phe").read()
        import pysnptools.util.pheno as pstpheno

        dict = pstpheno.loadOnePhen(
            self.currentFolder + "/examples/toydata.phe", missing=""
        )
        snpdata3 = Pheno(dict).read()
        np.testing.assert_array_almost_equal(snpdata1.val, snpdata3.val, decimal=10)

        dict = pstpheno.loadOnePhen(
            self.currentFolder + "/examples/toydata.phe", missing="", vectorize=True
        )
        assert len(dict["vals"].shape) == 1, "test 1-d array of values"
        snpdata3 = Pheno(dict).read()
        np.testing.assert_array_almost_equal(snpdata1.val, snpdata3.val, decimal=10)

        snpdata4 = Pheno(None, iid_if_none=snpdata1.iid)
        assert (snpdata4.row == snpdata1.row).all() and snpdata4.col_count == 0

        snpdata5 = Pheno(self.currentFolder + "/examples/toydata.id.phe").read()
        np.testing.assert_array_almost_equal(snpdata1.val, snpdata5.val, decimal=10)
        snpdata6 = Pheno(self.currentFolder + "/examples/toydata.fid.phe").read()
        np.testing.assert_array_almost_equal(snpdata1.val, snpdata6.val, decimal=10)

    def test_c_reader_dense(self):
        snpdata1 = self.snpdata[:, ::100].read()
        snpdata1.val[1, 2] = (
            np.nan
        )  # Inject a missing value to test writing and reading missing values
        output = "tempdir/snpreader/toydata.dense.txt"
        create_directory_if_necessary(output)
        Dense.write(output, snpdata1)
        snpreader = Dense(output)
        _fortesting_JustCheckExists().input(snpreader)
        snpdata2 = snpreader.read()
        np.testing.assert_array_almost_equal(snpdata1.val, snpdata2.val, decimal=10)

    def test_c_reader_distributedbed(self):
        from pysnptools.util.filecache import LocalCache

        snpdata1 = self.snpdata[:, ::100].read()
        snpdata1.val[1, 2] = (
            np.nan
        )  # Inject a missing value to test writing and reading missing values
        output = "tempdir/snpreader/toydata.distributedbed"
        LocalCache(output).rmtree()
        DistributedBed.write(output, snpdata1, piece_per_chrom_count=5)
        snpreader = DistributedBed(output)
        _fortesting_JustCheckExists().input(snpreader)
        snpdata2 = snpreader.read()
        np.testing.assert_array_almost_equal(snpdata1.val, snpdata2.val, decimal=10)

    def test_some_std(self):
        k0 = self.snpdata.read_kernel(standardizer=Unit()).val

        k1 = self.snpdata.read_kernel(standardizer=Unit())
        np.testing.assert_array_almost_equal(k0, k1.val, decimal=10)

        from pysnptools.snpreader import SnpData

        snpdata2 = SnpData(
            iid=self.snpdata.iid,
            sid=self.snpdata.sid,
            pos=self.snpdata.pos,
            val=np.array(self.snpdata.val),
        )
        _ = str(snpdata2)
        snpdata2.standardize()
        _ = str(snpdata2)

        snpreader = Bed(
            self.currentFolder + "/examples/toydata.5chrom.bed", count_A1=False
        )
        k2 = snpreader.read_kernel(standardizer=Unit(), block_size=500).val
        np.testing.assert_array_almost_equal(k0, k2, decimal=10)

        from pysnptools.standardizer.identity import Identity
        from pysnptools.standardizer.diag_K_to_N import DiagKtoN

        for dtype in [np.float64, np.float32]:
            for std in [Unit(), Beta(1, 25), Identity(), DiagKtoN()]:
                _ = str(std)
                np.random.seed(0)
                x = np.array(np.random.randint(3, size=[60, 100]), dtype=dtype)
                x2 = x[:, ::2]
                _ = np.array(x2)
                # LATER what's this about? It doesn't do non-contiguous?
                # assert not x2.flags['C_CONTIGUOUS'] and not x2.flags['F_CONTIGUOUS'] #set up to test non contiguous
                # assert x2b.flags['C_CONTIGUOUS'] or x2b.flags['F_CONTIGUOUS'] #set up to test non contiguous
                # a,b = std.standardize(x2b),std.standardize(x2)
                # np.testing.assert_array_almost_equal(a,b)

    def c_reader(self, snpreader):
        """
        make sure c-reader yields same result
        """
        snpdata = snpreader.read(order="F", force_python_only=False)
        snp_c = snpdata.val

        self.assertEqual(np.float64, snp_c.dtype)
        self.assertTrue(np.allclose(self.snps, snp_c, rtol=1e-05, atol=1e-05))
        return snpdata

    def test_standardize_bed(self):
        snpreader = Bed(
            self.currentFolder + "/examples/toydata.5chrom.bed", count_A1=False
        )
        self.standardize(snpreader)

    def test_standardize_dat(self):
        snpreader = Dat(self.currentFolder + "/examples/toydata.dat")
        self.standardize(snpreader)

    def test_standardize_ped(self):
        snpreader = Ped(self.currentFolder + "/examples/toydata")
        self.standardize(snpreader)

    def standardize(self, snpreader):
        """
        make sure blocked standardize yields same result as regular standardize
        """

        for dtype in [np.float64, np.float32]:
            snps = snpreader.read(order="F", force_python_only=True, dtype=dtype).val
            self.assertEqual(dtype, snps.dtype)

            snp_s1 = Unit().standardize(snps.copy(), force_python_only=True)
            snp_s2 = Unit().standardize(
                snps.copy(), block_size=100, force_python_only=True
            )
            snps_F = np.array(snps, dtype=dtype, order="F")
            snp_s3 = Unit().standardize(snps_F)
            snps_C = np.array(snps, dtype=dtype, order="C")
            snp_s4 = Unit().standardize(snps_C)

            self.assertEqual(snp_s1.shape[0], snp_s2.shape[0])
            self.assertEqual(snp_s1.shape[1], snp_s2.shape[1])

            self.assertEqual(snp_s1.shape[0], snp_s3.shape[0])
            self.assertEqual(snp_s1.shape[1], snp_s3.shape[1])

            self.assertEqual(snp_s1.shape[0], snp_s4.shape[0])
            self.assertEqual(snp_s1.shape[1], snp_s4.shape[1])

            self.assertTrue(np.allclose(snp_s1, snp_s2, rtol=1e-05, atol=1e-05))
            self.assertTrue(np.allclose(snp_s1, snp_s3, rtol=1e-05, atol=1e-05))
            self.assertTrue(np.allclose(snp_s1, snp_s4, rtol=1e-05, atol=1e-05))

            snp_beta1 = Beta(1, 25).standardize(snps.copy(), force_python_only=True)
            snps_F = np.array(snps, dtype=dtype, order="F")
            snp_beta2 = Beta(1, 25).standardize(snps_F)
            snps_C = np.array(snps, dtype=dtype, order="C")
            snp_beta3 = Beta(1, 25).standardize(snps_C)

            self.assertEqual(snp_beta1.shape[0], snp_beta2.shape[0])
            self.assertEqual(snp_beta1.shape[1], snp_beta2.shape[1])
            self.assertEqual(snp_beta1.shape[0], snp_beta3.shape[0])
            self.assertEqual(snp_beta1.shape[1], snp_beta3.shape[1])

            self.assertTrue(np.allclose(snp_beta1, snp_beta2, rtol=1e-05, atol=1e-05))
            self.assertTrue(np.allclose(snp_beta1, snp_beta3, rtol=1e-05, atol=1e-05))

    def test_load_and_standardize_bed(self):
        snpreader2 = Bed(
            self.currentFolder + "/examples/toydata.5chrom.bed", count_A1=False
        )
        self.load_and_standardize(snpreader2, snpreader2)

        snpreader3 = Bed(
            self.currentFolder + "/examples/toydata.5chrom", count_A1=False
        )
        self.load_and_standardize(snpreader3, snpreader3)

    def too_slow_test_write_bedbig(self):
        iid_count = 100000
        sid_count = 50000
        from pysnptools.snpreader import SnpData

        iid = np.array([[str(i), str(i)] for i in range(iid_count)])
        sid = np.array(["sid_{0}".format(i) for i in range(sid_count)])
        pos = np.array([[i, i, i] for i in range(sid_count)])
        np.random.seed(0)
        snpdata = SnpData(
            iid, sid, np.zeros((iid_count, sid_count)), pos=pos
        )  # random.choice((0.0,1.0,2.0,float("nan")),size=(iid_count,sid_count)))
        output = "tempdir/bedbig.{0}.{1}".format(iid_count, sid_count)
        create_directory_if_necessary(output)
        Bed.write(output, snpdata, count_A1=False)
        snpdata2 = Bed(output, count_A1=False).read()
        np.testing.assert_array_almost_equal(snpdata.val, snpdata2.val, decimal=10)

    def test_write_bed_f64cpp_0(self):
        snpreader = Bed(
            self.currentFolder + "/examples/toydata.5chrom.bed", count_A1=False
        )
        iid_index = 0
        logging.info("iid={0}".format(iid_index))
        # if snpreader.iid_count % 4 == 0: # divisible by 4 isn't a good test
        #    snpreader = snpreader[0:-1,:]
        # assert snpreader.iid_count % 4 != 0
        snpdata = snpreader[0:iid_index, :].read(order="F", dtype=np.float64)
        if snpdata.iid_count > 0:
            snpdata.val[-1, 0] = float("NAN")
        output = "tempdir/toydata.F64cpp.{0}.bed".format(iid_index)
        create_directory_if_necessary(output)
        Bed.write(output, snpdata, count_A1=False)
        snpdata2 = Bed(output, count_A1=False).read()
        np.testing.assert_array_almost_equal(snpdata.val, snpdata2.val, decimal=10)

    def test_write_bed_f64cpp_1(self):
        snpreader = Bed(
            self.currentFolder + "/examples/toydata.5chrom.bed", count_A1=False
        )
        iid_index = 1
        logging.info("iid={0}".format(iid_index))
        # if snpreader.iid_count % 4 == 0: # divisible by 4 isn't a good test
        #    snpreader = snpreader[0:-1,:]
        # assert snpreader.iid_count % 4 != 0
        snpdata = snpreader[0:iid_index, :].read(order="F", dtype=np.float64)
        if snpdata.iid_count > 0:
            snpdata.val[-1, 0] = float("NAN")
        output = "tempdir/toydata.F64cpp.{0}".format(iid_index)
        create_directory_if_necessary(output)
        Bed.write(output, snpdata, count_A1=False)
        snpdata2 = Bed(output, count_A1=False).read()
        np.testing.assert_array_almost_equal(snpdata.val, snpdata2.val, decimal=10)

    def test_write_bed_f64cpp_5(self):
        snpreader = Bed(
            self.currentFolder + "/examples/toydata.5chrom.bed", count_A1=False
        )

        from pysnptools.kernelreader.test import _fortesting_JustCheckExists

        _fortesting_JustCheckExists().input(snpreader)

        iid_index = 5
        logging.info("iid={0}".format(iid_index))
        # if snpreader.iid_count % 4 == 0: # divisible by 4 isn't a good test
        #    snpreader = snpreader[0:-1,:]
        # assert snpreader.iid_count % 4 != 0
        snpdata = snpreader[0:iid_index, :].read(order="F", dtype=np.float64)
        if snpdata.iid_count > 0:
            snpdata.val[-1, 0] = float("NAN")
        output = "tempdir/toydata.F64cpp.{0}".format(iid_index)
        create_directory_if_necessary(output)
        Bed.write(output, snpdata, count_A1=False)  # ,force_python_only=True)
        snpdata2 = Bed(output, count_A1=False).read()
        np.testing.assert_array_almost_equal(snpdata.val, snpdata2.val, decimal=10)

    def test_write_bed_f64cpp_5_python(self):
        snpreader = Bed(
            self.currentFolder + "/examples/toydata.5chrom.bed", count_A1=False
        )
        iid_index = 5
        logging.info("iid={0}".format(iid_index))
        # if snpreader.iid_count % 4 == 0: # divisible by 4 isn't a good test
        #    snpreader = snpreader[0:-1,:]
        # assert snpreader.iid_count % 4 != 0
        snpdata = snpreader[0:iid_index, :].read(order="F", dtype=np.float64)
        if snpdata.iid_count > 0:
            snpdata.val[-1, 0] = float("NAN")
        output = "tempdir/toydata.F64python.{0}".format(iid_index)
        create_directory_if_necessary(output)
        Bed.write(output, snpdata, force_python_only=True)
        snpdata2 = Bed(output, count_A1=False).read()
        np.testing.assert_array_almost_equal(snpdata.val, snpdata2.val, decimal=10)

    def test_write_x_x_cpp(self):
        for count_A1 in [False, True]:
            snpreader = Bed(
                self.currentFolder + "/examples/toydata.5chrom.bed", count_A1=count_A1
            )
            for order in ["C", "F"]:
                for dtype in [np.float32, np.float64]:
                    snpdata = snpreader.read(order=order, dtype=dtype)
                    snpdata.val[-1, 0] = float("NAN")
                    output = "tempdir/toydata.{0}{1}.cpp".format(
                        order, "32" if dtype == np.float32 else "64"
                    )
                    create_directory_if_necessary(output)
                    Bed.write(output, snpdata, count_A1=count_A1)
                    snpdata2 = Bed(output, count_A1=count_A1).read()
                    np.testing.assert_array_almost_equal(
                        snpdata.val, snpdata2.val, decimal=10
                    )

    def test_subset_view(self):
        snpreader2 = Bed(
            self.currentFolder + "/examples/toydata.5chrom.bed", count_A1=False
        )[:, :]
        result = snpreader2.read(view_ok=True)
        self.assertFalse(snpreader2 is result)
        result2 = result[:, :].read()
        self.assertFalse(np.may_share_memory(result2.val, result.val))
        result3 = result[:, :].read(view_ok=True)
        self.assertTrue(np.may_share_memory(result3.val, result.val))
        result4 = result3.read()
        self.assertFalse(np.may_share_memory(result4.val, result3.val))
        result5 = result4.read(view_ok=True)
        self.assertTrue(np.may_share_memory(result4.val, result5.val))

    @_requires_hdf5
    def test_load_and_standardize_hdf5(self):
        snpreader2 = SnpHdf5(self.currentFolder + "/examples/toydata.snpmajor.snp.hdf5")
        snpreader3 = SnpHdf5(self.currentFolder + "/examples/toydata.iidmajor.snp.hdf5")
        self.load_and_standardize(snpreader2, snpreader3)
        snpreaderref = Bed(
            self.currentFolder + "/examples/toydata.5chrom.bed", count_A1=False
        )
        self.load_and_standardize(snpreader2, snpreaderref)

    def test_load_and_standardize_dat(self):
        snpreader2 = Dat(self.currentFolder + "/examples/toydata.dat")
        self.load_and_standardize(snpreader2, snpreader2)
        # snpreaderref = Bed(self.currentFolder + "/examples/toydata.5chrom.bed",count_A1=False)
        # self.load_and_standardize(snpreader2, snpreaderref)

    def test_load_and_standardize_ped(self):
        # !!Ped columns can be ambiguous
        # ##Creating Ped data ...
        # currentFolder = os.path.dirname(os.path.realpath(__file__))
        # snpData = Bed(currentFolder + "/examples/toydata.5chrom.bed",count_A1=False).read()
        # #Ped.write(snpData, currentFolder + "/examples/toydata.ped")
        # fromPed = Ped(currentFolder + "/examples/toydata").read()
        # self.assertTrue(np.allclose(snpData.val, fromPed.val, rtol=1e-05, atol=1e-05))

        snpreader2 = Ped(self.currentFolder + "/examples/toydata")
        self.load_and_standardize(snpreader2, snpreader2)
        # snpreaderref = Bed(self.currentFolder + "/examples/toydata.5chrom.bed",count_A1=False)
        # self.load_and_standardize(snpreader2, snpreaderref)

    def load_and_standardize(self, snpreader2, snpreader3):
        """
        test c-version of load and standardize
        """

        S = snpreader2.sid_count
        N_original = snpreader2.iid_count

        iid_index_list = range(N_original - 1, 0, -2)
        snpreader3 = snpreader3[iid_index_list, :]

        for dtype in [np.float64, np.float32]:
            G2 = snpreader2.read(order="F", force_python_only=True).val
            G2 = Unit().standardize(G2, block_size=10000, force_python_only=True)

            SNPs_floatF = snpreader2.read(
                order="F", dtype=dtype, force_python_only=False
            ).val
            GF = Unit().standardize(SNPs_floatF)

            SNPs_floatC = snpreader2.read(
                order="C", dtype=dtype, force_python_only=False
            ).val
            GC = Unit().standardize(SNPs_floatC)

            self.assertTrue(np.allclose(GF, G2, rtol=1e-05, atol=1e-05))
            self.assertTrue(np.allclose(GF, GC, rtol=1e-05, atol=1e-05))

            # testing selecting a subset of snps and iids
            snp_index_list = range(S - 1, 0, -2)

            G2x = snpreader2.read(order="F", force_python_only=True).val
            G2x = G2x[iid_index_list, :][:, snp_index_list]
            G2x = Unit().standardize(G2x, block_size=10000, force_python_only=True)

            SNPs_floatFx = (
                snpreader3[:, snp_index_list]
                .read(order="F", dtype=dtype, force_python_only=False)
                .val
            )
            GFx = Unit().standardize(SNPs_floatFx)
            self.assertTrue(np.allclose(GFx, G2x, rtol=1e-05, atol=1e-05))

            SNPs_floatCx = (
                snpreader3[:, snp_index_list]
                .read(order="C", dtype=dtype, force_python_only=False)
                .val
            )
            _ = Unit().standardize(SNPs_floatCx)
            self.assertTrue(np.allclose(GFx, G2x, rtol=1e-05, atol=1e-05))

    def test_val_is_float(self):
        from pysnptools.snpreader import SnpData
        from pysnptools.kernelreader import KernelData

        iid_count = 3
        sid_count = 2
        iid = [[str(i), str(i)] for i in range(iid_count)]
        sid = [str(s) for s in range(sid_count)]
        np.random.seed(0)
        vali = np.random.randint(10, size=[iid_count, sid_count])
        snpdata = SnpData(iid=iid, sid=sid, val=vali)
        assert snpdata.val.dtype == np.float64
        snpdata = SnpData(iid=iid, sid=sid, val=np.array(vali, dtype=np.float32))
        assert snpdata.val.dtype == np.float32
        snpdata = SnpData(iid=iid, sid=sid, val=np.array(vali, dtype=np.float64))
        assert snpdata.val.dtype == np.float64
        snpdata = SnpData(iid=iid, sid=sid, val=np.array(vali, dtype=np.float16))
        assert snpdata.val.dtype == np.float64
        assert (
            np.float64
            == SnpData(iid=iid, sid=sid, val=np.array(vali, dtype=np.str_)).val.dtype
        )

        valstr = np.array(vali, dtype=np.str_)
        valstr[0, 0] = "cannot convert"
        error_seen = False
        try:
            SnpData(iid=iid, sid=sid, val=valstr)  # expect error
        except Exception:
            error_seen = True
        assert error_seen
        error_seen = False
        try:
            KernelData(iid=iid, val=valstr)  # expect error
        except Exception:
            error_seen = True
        assert error_seen

    def test_read_dtype(self):
        previous_wd = os.getcwd()
        os.chdir(os.path.dirname(os.path.realpath(__file__)))

        bed = Bed("examples/toydata.5chrom.bed", count_A1=True)
        bed[:3, :2].read(dtype="float32").val

        os.chdir(previous_wd)

    def test_respect_read_inputs(self):
        from pysnptools.snpreader import _MergeIIDs, _MergeSIDs, SnpGen, SnpMemMap

        previous_wd = os.getcwd()
        os.chdir(os.path.dirname(os.path.realpath(__file__)))

        snpreader_list = [
            Bed("examples/toydata.5chrom.bed", count_A1=True),
            SnpGen(seed=0, iid_count=500, sid_count=50),
            *([SnpHdf5("examples/toydata.snpmajor.snp.hdf5")] if _has_hdf5 else []),
            SnpMemMap("examples/tiny.snp.memmap"),
            SnpNpz("examples/toydata10.snp.npz"),
            Bed("examples/toydata.5chrom.bed", count_A1=True)[::2, ::2],
            _MergeSIDs(
                [
                    Bed("examples/toydata.5chrom.bed", count_A1=True)[:, :5].read(),
                    Bed("examples/toydata.5chrom.bed", count_A1=True)[:, 5:].read(),
                ]
            ),
            Dat("examples/toydata.dat"),
            Dense("examples/toydata100.dense.txt"),
            DistributedBed("examples/toydataSkip10.distributedbed"),
            Ped("examples/toydata.ped"),
            Pheno("examples/toydata.phe"),
            Bed("examples/toydata.5chrom.bed", count_A1=True).read(),
            _MergeIIDs(
                [
                    Bed("examples/toydata.5chrom.bed", count_A1=True)[:5, :].read(),
                    Bed("examples/toydata.5chrom.bed", count_A1=True)[5:, :].read(),
                ]
            ),
        ]

        output_p = "tempdir/snpreader/respect.p"
        create_directory_if_necessary(output_p)

        for snpreader in snpreader_list:
            logging.info(str(snpreader))
            for order in ["F", "C", "A"]:
                for dtype in [np.float32, np.float64]:
                    for force_python_only in [True, False]:
                        for view_ok in [True, False]:
                            val = snpreader.read(
                                order=order,
                                dtype=dtype,
                                force_python_only=force_python_only,
                                view_ok=view_ok,
                            ).val
                            has_right_order = (
                                order == "A"
                                or (order == "C" and val.flags["C_CONTIGUOUS"])
                                or (order == "F" and val.flags["F_CONTIGUOUS"])
                            )
                            if hasattr(snpreader, "val") and not view_ok:
                                assert snpreader.val is not val
                            if (
                                hasattr(snpreader, "val")
                                and view_ok
                                and snpreader.val is not val
                                and (
                                    order == "A"
                                    or (
                                        order == "F"
                                        and snpreader.val.flags["F_CONTIGUOUS"]
                                    )
                                    or (
                                        order == "C"
                                        and snpreader.val.flags["C_CONTIGUOUS"]
                                    )
                                )
                                and (dtype is None or snpreader.val.dtype == dtype)
                            ):
                                logging.info(
                                    "{0} could have read a view, but didn't".format(
                                        snpreader
                                    )
                                )
                            if not force_python_only:
                                # Don't check this when force_python_only -- Bed is know to get the order wrong,
                                # but it doesn't matter
                                assert val.dtype == dtype and has_right_order
            with open(output_p, "wb") as f:
                pickle.dump(snpreader, f)
            with open(output_p, "rb") as f:
                snpreader_p = pickle.load(f)
            val_p = snpreader_p.read(
                order=order,
                dtype=dtype,
                force_python_only=force_python_only,
                view_ok=view_ok,
            ).val
            np.allclose(val, val_p, equal_nan=True)

        os.chdir(previous_wd)

    def test_writes(self):
        from pysnptools.snpreader import SnpData, SnpMemMap

        the_class_and_suffix_list = [(SnpMemMap, "memmap", None, None)]

        # ===================================
        #    Starting main function
        # ===================================
        logging.info("starting 'test_writes'")
        np.random.seed(0)
        output_template = "tempdir/snpreader/writes.{0}.{1}"
        create_directory_if_necessary(output_template.format(0, "npz"))
        for row_count in [1]:
            for col_count in [4]:
                val = np.random.randint(0, 4, size=(row_count, col_count)) * 1.0
                val[val == 3] = np.nan
                row = [("0", "0"), ("1", "1"), ("2", "2"), ("3", "3"), ("4", "4")][
                    :row_count
                ]
                col = ["s0", "s1", "s2", "s3", "s4"][:col_count]
                for is_none in [True]:
                    col_prop = None
                    snpdata = SnpData(iid=row, sid=col, val=val, pos=col_prop)
                    for (
                        the_class,
                        suffix,
                        constructor,
                        writer,
                    ) in the_class_and_suffix_list:

                        def constructor(filename):
                            return the_class(filename)

                        def writer(filename, _data):
                            return the_class.write(filename, _data)

                        filename = output_template.format("debug", suffix)
                        logging.info(filename)
                        ret = writer(filename, snpdata)
                        assert ret is not None
                        for subsetter in [np.s_[::2, ::3]]:
                            reader = constructor(filename)
                            _fortesting_JustCheckExists().input(reader)
                            subreader = (
                                reader
                                if subsetter is None
                                else reader[subsetter[0], subsetter[1]]
                            )
                            _ = subreader.read(order="C")
                        try:
                            os.remove(filename)
                        except Exception:
                            pass
        logging.info("done with 'test_writes'")

    def test_writes2(self):
        from pysnptools.snpreader import SnpData, SnpHdf5, SnpNpz, SnpMemMap

        the_class_and_suffix_list = [
            (DistributedBed, "distributed_bed", None, None),
            (Dense, "dense", None, None),
            (Bed, "bed", lambda filename: Bed(filename, count_A1=False), None),
            (Dat, "dat", None, None),
            (Ped, "ped", None, None),
            (Pheno, "pheno", None, None),
            *([(SnpHdf5, "hdf5", None, None)] if _has_hdf5 else []),
            (SnpNpz, "npz", None, None),
            (SnpMemMap, "memmap", None, None),
        ]
        cant_do_col_prop_none_set = {"dense", "distributed_bed"}
        cant_do_col_len_0_set = {"distributed_bed"}
        cant_do_row_count_zero_set = {"dense", "ped", "pheno"}
        can_swap_0_2_set = {"ped"}
        can_change_col_names_set = {"pheno"}
        ignore_fam_id_set = {"dense"}
        ignore_pos_set = {"dense", "pheno"}
        erase_any_write_dir = {"distributed_bed"}

        # ===================================
        #    Starting main function
        # ===================================
        logging.info("starting 'test_writes'")
        np.random.seed(0)
        output_template = "tempdir/snpreader/writes.{0}.{1}"
        create_directory_if_necessary(output_template.format(0, "npz"))
        i = 0
        for row_count in [0, 5, 2, 1]:
            for col_count in [4, 2, 1, 0]:
                val = np.random.randint(0, 4, size=(row_count, col_count)) * 1.0
                val[val == 3] = np.nan
                row = [("0", "0"), ("1", "1"), ("2", "2"), ("3", "3"), ("4", "4")][
                    :row_count
                ]
                col = ["s0", "s1", "s2", "s3", "s4"][:col_count]
                for is_none in [True, False]:
                    if is_none:
                        col_prop = None
                    else:
                        col_prop = [(x, x, x) for x in range(1, 6)][:col_count]
                    snpdata = SnpData(
                        iid=row, sid=col, val=val, pos=col_prop, name=str(i)
                    )
                    for (
                        the_class,
                        suffix,
                        constructor,
                        writer,
                    ) in the_class_and_suffix_list:
                        constructor = constructor or (
                            lambda filename: the_class(filename)
                        )
                        writer = writer or (
                            lambda filename, _data: the_class.write(filename, _data)
                        )

                        if col_count == 0 and suffix in cant_do_col_len_0_set:
                            continue
                        if col_prop is None and suffix in cant_do_col_prop_none_set:
                            continue
                        if row_count == 0 and suffix in cant_do_row_count_zero_set:
                            continue
                        filename = output_template.format(i, suffix)
                        logging.info(filename)
                        i += 1
                        if suffix in erase_any_write_dir and os.path.exists(filename):
                            shutil.rmtree(filename)
                        ret = writer(filename, snpdata)
                        assert ret is not None
                        for subsetter in [None, np.s_[::2, ::3]]:
                            reader = constructor(filename)
                            _fortesting_JustCheckExists().input(reader)
                            subreader = (
                                reader
                                if subsetter is None
                                else reader[subsetter[0], subsetter[1]]
                            )
                            readdata = subreader.read(order="C")
                            expected = (
                                snpdata
                                if subsetter is None
                                else snpdata[subsetter[0], subsetter[1]].read()
                            )
                            if suffix not in can_swap_0_2_set:
                                assert np.allclose(
                                    readdata.val, expected.val, equal_nan=True
                                )
                            else:
                                for col_index in range(readdata.col_count):
                                    assert np.allclose(
                                        readdata.val[:, col_index],
                                        expected.val[:, col_index],
                                        equal_nan=True,
                                    ) or np.allclose(
                                        readdata.val[:, col_index] * -1 + 2,
                                        expected.val[:, col_index],
                                        equal_nan=True,
                                    )
                            if suffix not in ignore_fam_id_set:
                                assert np.array_equal(readdata.row, expected.row)
                            else:
                                assert np.array_equal(
                                    readdata.row[:, 1], expected.row[:, 1]
                                )
                            if suffix not in can_change_col_names_set:
                                assert np.array_equal(readdata.col, expected.col)
                            else:
                                assert readdata.col_count == expected.col_count
                            assert np.array_equal(
                                readdata.row_property, expected.row_property
                            ) or (
                                readdata.row_property.shape[1] == 0
                                and expected.row_property.shape[1] == 0
                            )

                            if suffix not in ignore_pos_set:
                                assert np.allclose(
                                    readdata.col_property,
                                    expected.col_property,
                                    equal_nan=True,
                                ) or (
                                    readdata.col_property.shape[1] == 0
                                    and expected.col_property.shape[1] == 0
                                )
                            else:
                                assert len(readdata.col_property) == len(
                                    expected.col_property
                                )
                        try:
                            os.remove(filename)
                        except Exception:
                            pass
        logging.info("done with 'test_writes'")


class NaNCNCTestCases(unittest.TestCase):
    def __init__(
        self,
        iid_index_list,
        snp_index_list,
        standardizer,
        snpreader,
        dtype,
        order,
        force_python_only,
        reference_snps,
        reference_dtype,
    ):
        self.iid_index_list = iid_index_list
        self.snp_index_list = snp_index_list
        self.standardizer = standardizer
        self.snpreader = snpreader
        self.dtype = dtype
        self.order = order
        self.force_python_only = force_python_only
        self.reference_snps = reference_snps
        self.reference_dtype = reference_dtype

    _testMethodName = "runTest"
    _testMethodDoc = None

    @staticmethod
    def factory_iterator():
        def snp_reader_factory_bed():
            return Bed("examples/toydata.5chrom.bed", count_A1=False)

        def snp_reader_factory_snpmajor_hdf5():
            return SnpHdf5("examples/toydata.snpmajor.snp.hdf5")

        def snp_reader_factory_iidmajor_hdf5():
            return SnpHdf5("examples/toydata.iidmajor.snp.hdf5")

        def snp_reader_factory_dat():
            return Dat("examples/toydata.dat")

        previous_wd = os.getcwd()
        os.chdir(os.path.dirname(os.path.realpath(__file__)))

        snpreader0 = snp_reader_factory_bed()
        S_original = snpreader0.sid_count
        N_original = snpreader0.iid_count

        snps_to_read_count = min(S_original, 100)

        for iid_index_list in [
            range(N_original),
            range(N_original // 2),
            range(N_original - 1, 0, -2),
        ]:
            for snp_index_list in [
                range(snps_to_read_count),
                range(snps_to_read_count // 2),
                range(snps_to_read_count - 1, 0, -2),
            ]:
                for standardizer in [Unit(), Beta(1, 25)]:
                    reference_snps, reference_dtype = NaNCNCTestCases(
                        iid_index_list,
                        snp_index_list,
                        standardizer,
                        snp_reader_factory_bed(),
                        np.float64,
                        "C",
                        "False",
                        None,
                        None,
                    ).read_and_standardize()
                    for snpreader_factory in [
                        snp_reader_factory_bed,
                        *([snp_reader_factory_snpmajor_hdf5,
                           snp_reader_factory_iidmajor_hdf5] if _has_hdf5 else []),
                        snp_reader_factory_dat,
                    ]:
                        for dtype in [np.float64, np.float32]:
                            for order in ["C", "F"]:
                                for force_python_only in [False, True]:
                                    snpreader = snpreader_factory()
                                    test_case = NaNCNCTestCases(
                                        iid_index_list,
                                        snp_index_list,
                                        standardizer,
                                        snpreader,
                                        dtype,
                                        order,
                                        force_python_only,
                                        reference_snps,
                                        reference_dtype,
                                    )
                                    yield test_case
        os.chdir(previous_wd)

    def __str__(self):
        iid_index_list = self.iid_index_list
        snp_index_list = self.snp_index_list
        standardizer = self.standardizer
        snpreader = self.snpreader
        dtype = self.dtype
        order = self.order
        force_python_only = self.force_python_only
        return (
            "{0}(iid_index_list=[{1}], snp_index_list=[{2}], standardizer={3}, snpreader={4}, "
            + "dtype={5}, order='{6}', force_python_only=={7})"
        ).format(
            self.__class__.__name__,
            (
                ",".join([str(i) for i in iid_index_list])
                if len(iid_index_list) < 10
                else ",".join([str(i) for i in iid_index_list[0:10]]) + ",..."
            ),
            (
                ",".join([str(i) for i in snp_index_list])
                if len(snp_index_list) < 10
                else ",".join([str(i) for i in snp_index_list[0:10]]) + ",..."
            ),
            standardizer,
            snpreader,
            dtype,
            order,
            force_python_only,
        )

    def read_and_standardize(self):
        previous_wd = os.getcwd()
        os.chdir(os.path.dirname(os.path.realpath(__file__)))
        iid_index_list = self.iid_index_list
        snp_index_list = self.snp_index_list
        standardizer = self.standardizer
        snpreader = self.snpreader
        dtype = self.dtype
        order = self.order
        force_python_only = self.force_python_only

        snps = (
            snpreader[iid_index_list, snp_index_list]
            .read(order=order, dtype=dtype, force_python_only=force_python_only)
            .val
        )
        snps[0, 0] = np.nan  # add a NaN
        snps[:, 1] = 2  # make a SNC
        snps = standardizer.standardize(snps, force_python_only=force_python_only)
        os.chdir(previous_wd)
        return snps, dtype

    def doCleanups(self):
        pass
        # return super(NaNCNCTestCases, self).doCleanups()

    def runTest(self, result=None):
        snps, dtype = self.read_and_standardize()
        self.assertTrue(snps[0, 0] == 0)
        self.assertTrue(np.all(snps[:, 1] == 0))
        if self.reference_snps is not None:
            self.assertTrue(
                np.allclose(
                    self.reference_snps,
                    snps,
                    rtol=(
                        1e-04
                        if dtype == np.float32 or self.reference_dtype == np.float32
                        else 1e-12
                    ),
                )
            )


# We do it this way instead of using doctest.DocTestSuite because doctest.DocTestSuite requires modules to be pickled,
# which python doesn't allow. We need tests to be pickleable so that they can be run on a cluster.
class TestSnpDocStrings(unittest.TestCase):
    def _run_doctest_without_progress_logs(self, module):
        root_logger = logging.getLogger()
        old_level = root_logger.level
        root_logger.setLevel(logging.WARN)
        try:
            return doctest.testmod(module, optionflags=doctest.ELLIPSIS)
        finally:
            root_logger.setLevel(old_level)

    def test_bed(self):
        import pysnptools.snpreader.bed

        old_dir = os.getcwd()
        os.chdir(os.path.dirname(os.path.realpath(__file__)) + "/snpreader")
        result = doctest.testmod(pysnptools.snpreader.bed, optionflags=doctest.ELLIPSIS)
        os.chdir(old_dir)
        assert result.failed == 0, "failed doc test: " + __file__

    def test_dat(self):
        import pysnptools.snpreader.snpdata

        old_dir = os.getcwd()
        os.chdir(os.path.dirname(os.path.realpath(__file__)) + "/snpreader")
        result = doctest.testmod(pysnptools.snpreader.dat, optionflags=doctest.ELLIPSIS)
        os.chdir(old_dir)
        assert result.failed == 0, "failed doc test: " + __file__

    def test_dense(self):
        import pysnptools.snpreader.snpdata

        old_dir = os.getcwd()
        os.chdir(os.path.dirname(os.path.realpath(__file__)) + "/snpreader")
        result = doctest.testmod(
            pysnptools.snpreader.dense, optionflags=doctest.ELLIPSIS
        )
        os.chdir(old_dir)
        assert result.failed == 0, "failed doc test: " + __file__

    def test_pairs(self):
        old_dir = os.getcwd()
        os.chdir(os.path.dirname(os.path.realpath(__file__)) + "/snpreader")
        import pysnptools.snpreader.pairs

        result = doctest.testmod(
            pysnptools.snpreader.pairs, optionflags=doctest.ELLIPSIS
        )
        os.chdir(old_dir)
        assert result.failed == 0, "failed doc test: " + __file__

    def test_ped(self):
        old_dir = os.getcwd()
        os.chdir(os.path.dirname(os.path.realpath(__file__)) + "/snpreader")
        import pysnptools.snpreader.ped

        result = doctest.testmod(pysnptools.snpreader.ped, optionflags=doctest.ELLIPSIS)
        os.chdir(old_dir)
        assert result.failed == 0, "failed doc test: " + __file__

    def test_pheno(self):
        old_dir = os.getcwd()
        os.chdir(os.path.dirname(os.path.realpath(__file__)) + "/snpreader")
        import pysnptools.snpreader.pheno

        result = doctest.testmod(
            pysnptools.snpreader.pheno, optionflags=doctest.ELLIPSIS
        )
        os.chdir(old_dir)
        assert result.failed == 0, "failed doc test: " + __file__

    def test_snpdata(self):
        old_dir = os.getcwd()
        os.chdir(os.path.dirname(os.path.realpath(__file__)) + "/snpreader")
        import pysnptools.snpreader.snpdata

        result = doctest.testmod(
            pysnptools.snpreader.snpdata, optionflags=doctest.ELLIPSIS
        )
        os.chdir(old_dir)
        assert result.failed == 0, "failed doc test: " + __file__

    @_requires_hdf5
    def test_snphdf5(self):
        old_dir = os.getcwd()
        os.chdir(os.path.dirname(os.path.realpath(__file__)) + "/snpreader")
        import pysnptools.snpreader.snphdf5

        result = doctest.testmod(
            pysnptools.snpreader.snphdf5, optionflags=doctest.ELLIPSIS
        )
        os.chdir(old_dir)
        assert result.failed == 0, "failed doc test: " + __file__

    def test_snpmemmap(self):
        old_dir = os.getcwd()
        os.chdir(os.path.dirname(os.path.realpath(__file__)) + "/snpreader")
        import pysnptools.snpreader.snpmemmap

        result = self._run_doctest_without_progress_logs(pysnptools.snpreader.snpmemmap)
        os.chdir(old_dir)
        assert result.failed == 0, "failed doc test: " + __file__

    def test_snpgen(self):
        old_dir = os.getcwd()
        os.chdir(os.path.dirname(os.path.realpath(__file__)) + "/util")
        import pysnptools.snpreader.snpgen

        result = self._run_doctest_without_progress_logs(pysnptools.snpreader.snpgen)
        os.chdir(old_dir)
        assert result.failed == 0, "failed doc test: " + __file__

    def test_snpnpz(self):
        old_dir = os.getcwd()
        os.chdir(os.path.dirname(os.path.realpath(__file__)) + "/snpreader")
        import pysnptools.snpreader.snpnpz

        result = doctest.testmod(
            pysnptools.snpreader.snpnpz, optionflags=doctest.ELLIPSIS
        )
        os.chdir(old_dir)
        assert result.failed == 0, "failed doc test: " + __file__

    def test_distributedbed(self):
        import pysnptools.snpreader.distributedbed

        old_dir = os.getcwd()
        os.chdir(os.path.dirname(os.path.realpath(__file__)) + "/snpreader")
        result = doctest.testmod(
            pysnptools.snpreader.distributedbed, optionflags=doctest.ELLIPSIS
        )
        os.chdir(old_dir)
        assert result.failed == 0, "failed doc test: " + __file__

    def test_snpreader(self):
        old_dir = os.getcwd()
        os.chdir(os.path.dirname(os.path.realpath(__file__)) + "/snpreader")
        import pysnptools.snpreader.snpreader

        result = doctest.testmod(
            pysnptools.snpreader.snpreader,
            optionflags=doctest.ELLIPSIS | doctest.NORMALIZE_WHITESPACE,
        )
        os.chdir(old_dir)
        assert result.failed == 0, "failed doc test: " + __file__

    def test_util(self):
        import pysnptools.util

        old_dir = os.getcwd()
        os.chdir(os.path.dirname(os.path.realpath(__file__)) + "/util")

        result = doctest.testmod(pysnptools.util, optionflags=doctest.ELLIPSIS)
        os.chdir(old_dir)
        assert result.failed == 0, "failed doc test: " + __file__

    def test_standardize_testmod(self):
        import pysnptools.standardizer

        old_dir = os.getcwd()
        os.chdir(os.path.dirname(os.path.realpath(__file__)) + "/standardizer")
        for mod in [
            pysnptools.standardizer.beta,
            pysnptools.standardizer.betatrained,
            pysnptools.standardizer.bysidcount,
            pysnptools.standardizer.bysqrtsidcount,
            pysnptools.standardizer.diag_K_to_N,
            pysnptools.standardizer.identity,
            pysnptools.standardizer.standardizer,
            pysnptools.standardizer.unit,
            pysnptools.standardizer.unittrained,
        ]:
            result = doctest.testmod(
                mod, optionflags=doctest.ELLIPSIS | doctest.NORMALIZE_WHITESPACE
            )
            assert result.failed == 0, "failed doc test: " + __file__
        os.chdir(old_dir)

    def test_kernelstandardize_testmod(self):
        import pysnptools.kernelstandardizer

        old_dir = os.getcwd()
        os.chdir(os.path.dirname(os.path.realpath(__file__)) + "/kernelstandardizer")
        for mod in [
            pysnptools.kernelstandardizer,
        ]:
            result = doctest.testmod(mod)
            assert result.failed == 0, "failed doc test: " + __file__
        os.chdir(old_dir)


def getTestSuite():
    """
    set up composite test suite
    """
    test_suite = unittest.TestSuite([])

    test_suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestBetaPython))
    test_suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestDiagKtoN))
    test_suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestPySnpTools))
    test_suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestDistReaders))
    test_suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestDistributedBed))
    test_suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestFileCache))
    test_suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestHashdownRetry))
    test_suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestRunnerThreads))
    test_suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestUtilTools))
    test_suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestIntRangeSet))
    test_suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestSnpDocStrings))
    test_suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestPstDocStrings))
    test_suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestKrDocStrings))
    test_suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestSnpGen))
    test_suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestGenerate))
    test_suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestExampleFile))
    test_suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestPstMemMap))
    test_suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestSnpMemMap))
    test_suite.addTests(NaNCNCTestCases.factory_iterator())
    test_suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestOptionalHdf5))
    test_suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestOptionalBgen))
    test_suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestPstReader))
    test_suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestKernelReader))

    return test_suite


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    suites = getTestSuite()
    r = unittest.TextTestRunner(failfast=False)
    ret = r.run(suites)
    assert ret.wasSuccessful()
