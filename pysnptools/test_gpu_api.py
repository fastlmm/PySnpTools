"""GPU API regressions. Requires the CUDA bed-reader build and a working CuPy GPU."""
import numpy as np
import pytest
from unittest.mock import patch

from bed_reader import sample_file
from pysnptools.snpreader import Bed, SnpData
from pysnptools.kernelreader import SnpKernel
from pysnptools.standardizer import Unit, Beta


@pytest.fixture
def cupy():
    # Optional GPU qualification is explicitly reported on CPU-only installations.
    from bed_reader import has_gpu_support
    if not has_gpu_support():
        pytest.skip("requires bed-reader built with the cubecl-cuda Cargo feature")
    return pytest.importorskip("cupy")


@pytest.mark.parametrize("dtype", [np.float32, np.float64])
@pytest.mark.parametrize("standardizer", [Unit(), Beta(1, 25)])
def test_subset_fused_and_blocked(cupy, dtype, standardizer):
    reader = Bed(sample_file("toydata.5chrom.bed"), count_A1=True)[::-3, ::97][::-1, ::-1]
    expected, trained = reader._read_with_standardizing(
        standardizer, dtype=dtype, xp=np, return_trained=True)
    with patch("bed_reader.bed_reader.standardize_gpu", wraps=__import__(
            "bed_reader.bed_reader", fromlist=["standardize_gpu"]).standardize_gpu) as native:
        actual, gpu_trained = reader._read_with_standardizing(
            standardizer, dtype=dtype, xp=cupy, return_trained=True)
        assert native.call_count == 1
        np.testing.assert_array_equal(cupy.asnumpy(actual.val), expected.val)
        np.testing.assert_array_equal(gpu_trained.stats, trained.stats)
        kernel = reader.read_kernel(standardizer, block_size=7, dtype=dtype, xp=cupy)
        assert native.call_count > 2  # the requested blocks were actually used
    # BLAS accumulation order differs between blocked GPU and whole CPU products.
    tolerance = 32 * np.finfo(dtype).eps
    np.testing.assert_allclose(cupy.asnumpy(kernel.val), expected.val @ expected.val.T,
                               atol=tolerance * reader.sid_count, rtol=tolerance)
    for to_kerneldata in [False, True]:
        result = SnpKernel(reader, standardizer, block_size=7)._read_with_standardizing(
            to_kerneldata, dtype=dtype, xp=cupy)
        val = result.val if to_kerneldata else result.snpreader.val
        assert val.dtype == dtype and isinstance(val, cupy.ndarray)


@pytest.mark.parametrize("dtype", [np.float32, np.float64])
@pytest.mark.parametrize("order", ["C", "F"])
@pytest.mark.parametrize("standardizer", [Beta(1, 25), Beta(0.5, 0.5)])
def test_blocked_beta_trained_bed(cupy, tmp_path, dtype, order, standardizer):
    from bed_reader import to_bed

    rng = np.random.default_rng(37)
    values = rng.integers(0, 3, size=(13, 9)).astype(dtype)
    values[:8, :2] = [0, 2]  # training constants can vary in the test cohort
    values[8:, 4] = np.nan  # applying moments must not refit an empty test SNP
    path = tmp_path / "trained.bed"
    to_bed(path, values, count_A1=True)
    reader = Bed(path, count_A1=True)
    _, trained = reader[:8, :].read(dtype=dtype, xp=np).standardize(standardizer, return_trained=True)
    test_reader = reader[8:, :][::-1, :]
    standardized = test_reader.read(dtype=dtype, xp=np).standardize(trained).val
    expected = standardized @ standardized.T
    # Nine well-scaled SNP contributions and three block additions round in dtype.
    tolerance = 64 * np.finfo(dtype).eps
    for xp in (np, cupy):
        kernel = test_reader.read_kernel(trained, block_size=4, dtype=dtype, order=order, xp=xp)
        assert isinstance(kernel.val, xp.ndarray) and kernel.val.dtype == dtype
        assert kernel.val.flags["C_CONTIGUOUS" if order == "C" else "F_CONTIGUOUS"]
        actual = kernel.val if xp is np else cupy.asnumpy(kernel.val)
        np.testing.assert_allclose(actual, expected, rtol=tolerance, atol=tolerance * reader.sid_count)


def test_explicit_module_and_existing_array(cupy):
    reader = Bed(sample_file("small.bed"), count_A1=True)
    with patch.dict("os.environ", {"ARRAY_MODULE": "cupy"}):
        assert isinstance(reader.read(xp=np).val, np.ndarray)
        assert isinstance(reader.read().val, cupy.ndarray)
        data = reader.read(xp=np)
        data.standardize()
        assert isinstance(data.val, np.ndarray)
        assert isinstance(reader.read_kernel(Unit(), xp=np).val, np.ndarray)
    with patch.dict("os.environ", {"ARRAY_MODULE": "numpy"}):
        data = reader.read(xp=cupy)
        data.standardize()
        assert isinstance(data.val, cupy.ndarray)
        assert isinstance(data[::2, :].read(xp=cupy).val, cupy.ndarray)


@pytest.mark.parametrize("dtype", [np.float32, np.float64])
@pytest.mark.parametrize("order", ["F", "C"])
def test_foreign_dosages_and_trained(cupy, dtype, order):
    rng = np.random.default_rng(1)
    val = rng.uniform(0.0, 2.0, (41, 6)).astype(dtype)
    val[::5, 0] = np.nan
    val[:, 1] = 1.0
    iid = np.array([[str(index), str(index)] for index in range(len(val))])
    sid = np.array([str(index) for index in range(val.shape[1])])
    for standardizer in [Unit(), Beta(1, 25)]:
        # Beta is only defined for allele frequencies in [0, 1].
        source = val if isinstance(standardizer, Unit) else val[:, :2]
        data_np = SnpData(iid=iid, sid=sid[:source.shape[1]], val=np.array(source, order=order), xp=np)
        data_cp = SnpData(iid=iid, sid=data_np.sid, val=cupy.array(source, order=order), xp=cupy)
        expected, trained_np = data_np.standardize(standardizer, return_trained=True)
        actual, trained_cp = data_cp.standardize(standardizer, return_trained=True)
        tolerance = 64 * np.finfo(dtype).eps
        np.testing.assert_allclose(cupy.asnumpy(actual.val), expected.val, rtol=tolerance, atol=tolerance)
        fresh_np = SnpData(iid=iid, sid=data_np.sid, val=np.array(source, order=order), xp=np).standardize(trained_np)
        fresh_cp = SnpData(iid=iid, sid=data_np.sid, val=cupy.array(source, order=order), xp=cupy).standardize(trained_cp)
        np.testing.assert_allclose(cupy.asnumpy(fresh_cp.val), fresh_np.val, rtol=tolerance, atol=tolerance)


@pytest.mark.parametrize("dtype", [np.float32, np.float64])
def test_large_offset_against_centered_reference(cupy, dtype):
    # The legacy Rust CPU standardizer uses E[x*x]-E[x]**2 and loses
    # precision here, so use an independent centered float64 reference.
    rng = np.random.default_rng(12)
    values = (1e4 + rng.random((100, 4))).astype(dtype)
    mean = values.astype(np.float64).mean(axis=0).astype(dtype)
    deviation = values.astype(np.float64).std(axis=0).astype(dtype)
    expected = (values - mean) * (1 / deviation)
    data = SnpData(iid=[[str(index), str(index)] for index in range(100)],
                   sid=["a", "b", "c", "d"], val=cupy.array(values, order="F"), xp=cupy)
    actual = data.standardize(Unit()).val
    # Parallel reductions can round the large mean by a few ulps. Subtracting
    # it and dividing by a small deviation amplifies that error by this ratio.
    absolute_tolerance = 8 * np.finfo(dtype).eps * np.max(np.abs(values)) / np.min(deviation)
    np.testing.assert_allclose(cupy.asnumpy(actual), expected,
                               rtol=8 * np.finfo(dtype).eps, atol=absolute_tolerance)


@pytest.mark.parametrize("order", ["C", "F"])
def test_standardize_workspace_is_not_matrix_sized(cupy, order):
    rng = np.random.default_rng(8)
    values = rng.integers(0, 3, size=(1000, 5000), dtype=np.int8)
    pool = cupy.cuda.MemoryPool()
    with cupy.cuda.using_allocator(pool.malloc):
        data = SnpData(
            iid=[[str(index), str(index)] for index in range(values.shape[0])],
            sid=[str(index) for index in range(values.shape[1])],
            val=cupy.array(values, dtype=np.float32, order=order), xp=cupy)
        cupy.cuda.get_current_stream().synchronize()
        baseline = pool.total_bytes()
        data.standardize(Unit())
        cupy.cuda.get_current_stream().synchronize()
        # The 19 MiB matrix must not need a matrix-sized mask/temporary.
        # Four MiB allows reduction scratch space and several per-SNP vectors.
        assert pool.total_bytes() - baseline < 4 * 2**20
        assert bool(cupy.isfinite(data.val).all())
