"""Exercise runner thread limits and BLAS pools visible to threadpoolctl."""
import json
import os
import pickle
import platform
from inspect import signature
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

import bed_reader
import numpy as np
from scipy import linalg
from threadpoolctl import threadpool_info, threadpool_limits

from pysnptools.util.mapreduce1 import map_reduce
from pysnptools.util.mapreduce1.mapreduce import _MapReduce
from pysnptools.util.mapreduce1.runner import LocalInParts, LocalMultiProc
from pysnptools.util.mapreduce1.runner import _threads


class _RunnerWithSlot(LocalMultiProc):
    __slots__ = ("label",)


_THREAD_VARIABLES = (
    "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS",
    "PST_NUM_THREADS", "BLIS_NUM_THREADS",
)


def _snapshot():
    # Load both libraries before inspecting their actual runtime thread counts.
    matrix = np.eye(3)
    np.dot(matrix, matrix)
    linalg.eigh(matrix)
    return {
        "pools": [(pool["filepath"], pool["num_threads"])
                  for pool in threadpool_info() if pool["user_api"] == "blas"],
        "environment": {name: os.environ.get(name) for name in _THREAD_VARIABLES},
        "bed_threads": bed_reader.get_num_threads(),
    }


def _mapper(value):
    return _snapshot()


def _reducer(sequence):
    return {"workers": list(sequence), "reducer": _snapshot()}


def _failing_mapper(value):
    raise RuntimeError("mapper failed while thread limits were active")


def _failing_reducer(sequence):
    list(sequence)
    raise RuntimeError("reducer failed while thread limits were active")


class TestRunnerThreads(unittest.TestCase):
    def setUp(self):
        temporary = TemporaryDirectory(prefix="pysnptools-runner-threads-")
        self.addCleanup(temporary.cleanup)
        previous_directory = os.getcwd()
        os.chdir(temporary.name)
        self.addCleanup(os.chdir, previous_directory)
        self.directory = Path(temporary.name)
        self.before = _snapshot()
        self.assertBlasPools(self.before["pools"])

    def assertBlasPools(self, pools):
        # macOS wheels use Accelerate, which threadpoolctl cannot inspect.
        # Windows ARM64 wheels can also lack a discoverable BLAS pool. Keep
        # requiring real pools on Linux and Windows x86, and check every pool
        # that is visible on any platform; an empty list is not a test skip.
        system = platform.system()
        windows_arm = system == "Windows" and platform.machine().lower() in ("arm64", "aarch64")
        if system != "Darwin" and not windows_arm:
            self.assertTrue(pools, "Requires an actual NumPy or SciPy BLAS pool")

    def assertLimited(self, snapshot, count):
        self.assertBlasPools(snapshot["pools"])
        self.assertTrue(all(threads == count for path, threads in snapshot["pools"]), snapshot)
        self.assertEqual(snapshot["bed_threads"], count)
        self.assertEqual(snapshot["environment"], {name: str(count) for name in _THREAD_VARIABLES})

    def test_platform_specific_pool_checks(self):
        for system, machine, requires_pool in (
            ("Darwin", "arm64", False), ("Darwin", "x86_64", False),
            ("Windows", "ARM64", False), ("Windows", "aarch64", False),
            ("Windows", "AMD64", True), ("Windows", "x86", True),
            ("Linux", "x86_64", True), ("Linux", "aarch64", True),
        ):
            with self.subTest(system=system, machine=machine):
                with patch.object(platform, "system", return_value=system), patch.object(platform, "machine", return_value=machine):
                    snapshot = {"pools": [], "bed_threads": 2,
                                "environment": {name: "2" for name in _THREAD_VARIABLES}}
                    if requires_pool:
                        with self.assertRaisesRegex(AssertionError, "Requires an actual"):
                            self.assertLimited(snapshot, 2)
                    else:
                        self.assertLimited(snapshot, 2)
                        with self.assertRaises(AssertionError):
                            self.assertLimited(dict(snapshot, bed_threads=1), 2)
                        with self.assertRaises(AssertionError):
                            self.assertLimited(dict(snapshot, environment={}), 2)
                    self.assertLimited(dict(snapshot, pools=[("blas", 2)]), 2)
                    with self.assertRaises(AssertionError):
                        self.assertLimited(dict(snapshot, pools=[("blas", 1)]), 2)

    def test_canonical_signature_and_positional_compatibility(self):
        for constructor, positional in ((LocalMultiProc, (2,)), (LocalInParts, (0, 2))):
            with self.subTest(constructor=constructor.__name__):
                self.assertIsNone(signature(constructor).parameters["num_threads"].default)
                self.assertEqual(constructor(*positional, 2).num_threads, 2)
                self.assertIsNone(constructor(*positional).num_threads)
                self.assertIsNone(constructor(*positional, num_threads=None).num_threads)
        # Parameters after the historical thread-count slot keep their positions.
        self.assertEqual(LocalMultiProc(2, 1, [1, 1]).weights, [1, 1])
        self.assertEqual(LocalInParts(0, 2, 1, [1, 1]).weights, [1, 1])

    def test_deprecated_constructor_keyword(self):
        for constructor, positional in ((LocalMultiProc, (1,)), (LocalInParts, (0, 1))):
            for count in (None, 1, 2):
                with self.subTest(constructor=constructor.__name__, count=count):
                    with self.assertWarnsRegex(DeprecationWarning, "mkl_num_threads.*num_threads"):
                        runner = constructor(*positional, mkl_num_threads=count)
                    self.assertEqual(runner.num_threads, count)
        with self.assertWarns(DeprecationWarning):
            runner = LocalMultiProc(1, mkl_num_threads=2)
        result = map_reduce(range(1), mapper=_mapper, reducer=_reducer, runner=runner)
        self.assertLimited(result["workers"][0], 2)
        self.assertEqual(_snapshot(), self.before)

    def test_constructor_rejects_both_names_even_with_none(self):
        for constructor, positional in ((LocalMultiProc, (1,)), (LocalInParts, (0, 1))):
            for canonical, legacy in ((None, None), (None, 1), (1, None), (1, 1), (1, 2)):
                with self.subTest(constructor=constructor.__name__, canonical=canonical, legacy=legacy):
                    with self.assertRaisesRegex(TypeError, "Specify only num_threads"):
                        constructor(*positional, num_threads=canonical, mkl_num_threads=legacy)
                    with self.assertRaisesRegex(TypeError, "Specify only num_threads"):
                        constructor(*positional, canonical, mkl_num_threads=legacy)

    def test_deprecated_attribute_remains_synchronized(self):
        for runner in (LocalMultiProc(1, num_threads=1), LocalInParts(0, 1, num_threads=1)):
            with self.subTest(runner=type(runner).__name__):
                with self.assertWarns(DeprecationWarning):
                    self.assertEqual(runner.mkl_num_threads, 1)
                with self.assertWarns(DeprecationWarning):
                    runner.mkl_num_threads = 2
                self.assertEqual(runner.num_threads, 2)
                runner.num_threads = None
                with self.assertWarns(DeprecationWarning):
                    self.assertIsNone(runner.mkl_num_threads)

    def test_legacy_and_current_pickle_states(self):
        for runner in (LocalMultiProc(1, num_threads=2), LocalInParts(0, 1, num_threads=2)):
            with self.subTest(runner=type(runner).__name__):
                self.assertEqual(runner.__getstate__()["mkl_num_threads"], 2)
                self.assertNotIn("num_threads", runner.__getstate__())
                self.assertEqual(pickle.loads(pickle.dumps(runner)).num_threads, 2)
                # The native pickle state keeps the field used by older releases.
                self.assertEqual(runner.__dict__["mkl_num_threads"], 2)
                restored = pickle.loads(pickle.dumps(runner))
                self.assertEqual(restored.num_threads, 2)
                self.assertNotIn("num_threads", restored.__dict__)
                result = map_reduce(range(1), mapper=_mapper, reducer=_reducer, runner=restored)
                if result is not None:
                    self.assertLimited(result["workers"][0], 2)
                self.assertEqual(_snapshot(), self.before)
    def test_pickle_preserves_subclass_slots(self):
        runner = _RunnerWithSlot(1, num_threads=2)
        runner.label = "preserved"
        restored = pickle.loads(pickle.dumps(runner))
        self.assertEqual(restored.num_threads, 2)
        self.assertEqual(restored.label, runner.label)

    def test_shared_limits_restore_loaded_pools_and_environment(self):
        for count in (1, 2):
            with self.subTest(count=count):
                with _threads.limits(count):
                    self.assertLimited(_snapshot(), count)
                self.assertEqual(_snapshot(), self.before)
        with self.assertRaisesRegex(RuntimeError, "sentinel"):
            with _threads.limits(1):
                self.assertLimited(_snapshot(), 1)
                raise RuntimeError("sentinel")
        self.assertEqual(_snapshot(), self.before)

    def test_none_is_strict_noop(self):
        # None must neither change settings nor create a runtime controller.
        with patch.object(_threads, "threadpool_limits", side_effect=AssertionError("None created a controller")):
            with _threads.limits(None):
                self.assertEqual(_snapshot(), self.before)
            self.assertEqual(_threads.environment(None), {})
        self.assertEqual(_snapshot(), self.before)

    def test_validate_thread_counts(self):
        for count in (True, False, np.bool_(True), np.bool_(False), 1.5, "2", object()):
            with self.subTest(count=count):
                with self.assertRaises(TypeError):
                    _threads.environment(count)
                with self.assertRaises(TypeError):
                    with _threads.limits(count):
                        self.fail("Invalid thread limit was accepted")
        for count in (0, -1):
            with self.subTest(count=count):
                with self.assertRaises(ValueError):
                    _threads.environment(count)
        self.assertEqual(_threads.environment(np.int64(2)), {name: "2" for name in _THREAD_VARIABLES})
        self.assertEqual(_snapshot(), self.before)

    def test_multiprocess_and_debug_routes(self):
        for just_one_process in (False, True):
            for count in (1, 2):
                with self.subTest(just_one_process=just_one_process, count=count):
                    result = map_reduce(
                        range(2), mapper=_mapper, reducer=_reducer,
                        runner=LocalMultiProc(2, num_threads=count,
                                              just_one_process=just_one_process),
                    )
                    self.assertEqual(len(result["workers"]), 2)
                    for snapshot in result["workers"] + [result["reducer"]]:
                        self.assertLimited(snapshot, count)
                    self.assertEqual(_snapshot(), self.before)

    def test_local_in_parts_limits_mapper_and_reducer(self):
        self.assertIsNone(map_reduce(
            range(1), mapper=_mapper, reducer=_reducer,
            runner=LocalInParts(0, 1, num_threads=1),
        ))
        self.assertEqual(_snapshot(), self.before)
        result = map_reduce(
            range(1), mapper=_mapper, reducer=_reducer,
            runner=LocalInParts(1, 1, num_threads=1),
        )
        self.assertLimited(result["workers"][0], 1)
        self.assertLimited(result["reducer"], 1)
        self.assertEqual(_snapshot(), self.before)

    def test_none_preserves_current_debug_and_local_in_parts_settings(self):
        with patch.dict(os.environ, {name: "2" for name in _THREAD_VARIABLES}), threadpool_limits(limits=2, user_api="blas"):
            before = _snapshot()
            for runner in (LocalMultiProc(1), LocalMultiProc(1, just_one_process=True), LocalInParts(0, 1)):
                with self.subTest(runner=type(runner).__name__):
                    result = map_reduce(range(1), mapper=_mapper, reducer=_reducer, runner=runner)
                    if result is not None:
                        self.assertLimited(result["workers"][0], 2)
                        self.assertLimited(result["reducer"], 2)
                    self.assertEqual(_snapshot(), before)
        self.assertEqual(_snapshot(), self.before)

    def test_mapper_failures_restore_parent(self):
        for runner in (LocalMultiProc(1, num_threads=1, just_one_process=True), LocalInParts(0, 1, num_threads=1)):
            with self.subTest(runner=type(runner).__name__):
                with self.assertRaisesRegex(RuntimeError, "mapper failed"):
                    map_reduce(range(1), mapper=_failing_mapper, runner=runner)
                self.assertEqual(_snapshot(), self.before)

    def test_reducer_failures_restore_parent(self):
        for just_one_process in (False, True):
            with self.subTest(just_one_process=just_one_process):
                with self.assertRaisesRegex(RuntimeError, "reducer failed"):
                    map_reduce(range(1), mapper=_mapper, reducer=_failing_reducer,
                               runner=LocalMultiProc(1, num_threads=1,
                                                     just_one_process=just_one_process))
                self.assertEqual(_snapshot(), self.before)
        # Exercise the independent LocalInParts reduction route too.
        distributable = _MapReduce(range(1), _mapper, None, _failing_reducer)
        LocalInParts(0, 1, num_threads=1).run(distributable)
        with self.assertRaisesRegex(RuntimeError, "reducer failed"):
            LocalInParts(1, 1, num_threads=1).run(distributable)
        self.assertEqual(_snapshot(), self.before)

    def test_subprocess_startup_limits_precede_numpy_scipy_imports(self):
        startup = self.directory / "startup"
        startup.mkdir()
        script = '''import json, os
from pathlib import Path
import numpy as np
from scipy import linalg
from threadpoolctl import threadpool_info
matrix = np.eye(3)
np.dot(matrix, matrix)
linalg.eigh(matrix)
record = {"environment": {name: os.environ.get(name) for name in %r},
          "pools": [(pool["filepath"], pool["num_threads"]) for pool in threadpool_info() if pool["user_api"] == "blas"]}
Path(os.environ["PYSNPTOOLS_THREAD_STARTUP"], str(os.getpid()) + ".json").write_text(json.dumps(record))
''' % (_THREAD_VARIABLES,)
        (startup / "sitecustomize.py").write_text(script)
        pythonpath = str(startup)
        if os.environ.get("PYTHONPATH"):
            pythonpath += os.pathsep + os.environ["PYTHONPATH"]
        with patch.dict(os.environ, {"PYTHONPATH": pythonpath, "PYSNPTOOLS_THREAD_STARTUP": str(startup)}):
            map_reduce(range(2), mapper=_mapper, reducer=_reducer,
                       runner=LocalMultiProc(2, num_threads=1))
        records = [json.loads(path.read_text()) for path in startup.glob("*.json")]
        self.assertEqual(len(records), 2, records)
        for record in records:
            self.assertBlasPools(record["pools"])
            self.assertTrue(all(threads == 1 for path, threads in record["pools"]), record)
            self.assertEqual(record["environment"], {name: "1" for name in _THREAD_VARIABLES})
        self.assertEqual(_snapshot(), self.before)
