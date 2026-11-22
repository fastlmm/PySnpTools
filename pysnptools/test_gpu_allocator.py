"""Use CuPy allocations without taking over caller-provided CuPy memory."""

from unittest.mock import patch

import numpy as np
import pytest

from pysnptools.kernelreader import KernelData
from pysnptools.snpreader import SnpData
from pysnptools.standardizer import Beta, Identity, Unit


@pytest.fixture
def cupy():
    return pytest.importorskip("cupy")


def data(values, xp):
    return SnpData(iid=[["f", str(sample)] for sample in range(len(values))],
                   sid=[str(marker) for marker in range(values.shape[1])], val=values, xp=xp)


@pytest.mark.parametrize("dtype", [np.float32, np.float64])
@pytest.mark.parametrize("order", ["C", "F"])
@pytest.mark.parametrize("standardizer", [Unit(), Beta(2, 3), Beta(1, 25), Beta(0.5, 0.5)])
def test_host_constructor_cupy_standardization(cupy, dtype, order, standardizer):

    values = np.array([[.25, 1, .5, 0, 2], [.75, 1, np.nan, 0, 2],
                       [1.75, 1, 1.5, 0, 2], [1.25, 1, .75, 0, 2]], dtype=dtype, order=order)
    if isinstance(standardizer, Unit):
        values[:, 0] -= 4
    original = values.copy()
    target = data(values, cupy)
    pointer = target.val.data.ptr
    observed_mean = np.nanmean(values.astype(np.float64), axis=0)
    if isinstance(standardizer, Unit):
        deviations = np.nanstd(values.astype(np.float64), axis=0)
        deviations[deviations == 0] = np.inf
        factors = 1 / deviations
    else:
        frequency = np.minimum(observed_mean / 2, 1 - observed_mean / 2)
        variable = np.nanstd(values.astype(np.float64), axis=0) > 0
        factors = np.zeros(len(observed_mean))
        frequency = frequency[variable]
        # Independent closed-form densities, evaluated only for variable SNPs.
        if (standardizer.a, standardizer.b) == (2, 3):
            factors[variable] = 12 * frequency * (1 - frequency) ** 2
        elif (standardizer.a, standardizer.b) == (1, 25):
            factors[variable] = 25 * (1 - frequency) ** 24
        else:
            factors[variable] = 1 / (np.pi * np.sqrt(frequency * (1 - frequency)))
    expected = (values - observed_mean) * factors
    expected[np.isnan(values)] = 0
    target.standardize(standardizer)
    assert target.val.data.ptr == pointer and target.val.dtype == dtype
    assert target.val.flags["C_CONTIGUOUS" if order == "C" else "F_CONTIGUOUS"]
    tolerance = 32 * np.finfo(dtype).eps  # bounded inputs and centered reductions
    if isinstance(standardizer, Beta) and standardizer.a == 1:
        # Retain the established Beta(1,25) closed-form comparison bound.
        tolerance = max(tolerance, 1e-12)
    np.testing.assert_allclose(cupy.asnumpy(target.val), expected, rtol=tolerance, atol=tolerance)
    np.testing.assert_array_equal(values, original)


def test_caller_arrays_and_views_keep_storage(cupy):

    for source in (cupy.array([[0., 1.], [1., 1.], [2., 1.]]),
                   cupy.empty((3, 4), order="F")[:, :2], cupy.empty((3, 2)).T):
        source[...] = cupy.arange(source.size).reshape(source.shape)
        target = data(source, cupy)
        assert target.val is source
        target.standardize()


def test_slice_standardizes_only_selected_columns(cupy):

    values = cupy.empty((3, 4), order="F")
    values[...] = cupy.array([[0, 1, 7, 8], [1, 1, 7, 8], [2, 1, 7, 8]])
    data(values[:, :2], cupy).standardize()
    np.testing.assert_array_equal(cupy.asnumpy(values[:, 2:]), [[7, 8]] * 3)


@pytest.mark.parametrize("order", ["C", "F"])
def test_snp_reads_are_cupy_and_kernels_use_cupy(cupy, order):
    # Read copies and kernels use the selected CuPy allocator.

    source = data(cupy.array([[0., 1.], [1., 1.], [2., 1.]]), cupy)
    copied = source.read(dtype="float64", order=order, xp=cupy)
    assert copied.val.data.ptr != source.val.data.ptr
    for block_size in (None, 1):
        kernel = source.read_kernel(Identity(), block_size=block_size, dtype="float64", order=order, xp=cupy)
        assert isinstance(kernel.val, cupy.ndarray)
        np.testing.assert_array_equal(cupy.asnumpy(kernel.val), [[1, 1, 1], [1, 2, 3], [1, 3, 5]])
        assert kernel.val.flags["C_CONTIGUOUS" if order == "C" else "F_CONTIGUOUS"]
    host = KernelData(iid=source.iid, val=np.eye(3), xp=np)
    on_gpu = host.read(xp=cupy).val
    assert isinstance(on_gpu, cupy.ndarray)


@pytest.mark.parametrize("dtype", [np.float32, np.float64])
@pytest.mark.parametrize("order", ["C", "F"])
@pytest.mark.parametrize("host_source", [False, True])
def test_subset_reads_use_cupy_standardization(cupy, dtype, order, host_source):

    values = np.array([[.25, 1, 5, np.nan], [.75, 1, 6, .5],
                       [np.nan, 1, 7, 1.5], [1.75, 1, 8, .75],
                       [1.25, 1, 9, 1.25]], dtype=dtype, order=order)
    source = data(values if host_source else cupy.array(values, order=order), cupy)
    rows, columns = [4, 1, 4, 0], [3, 0, 1]
    selected = source[rows, columns].read(dtype=dtype, order=order, xp=cupy)
    expected = values[np.ix_(rows, columns)]
    np.testing.assert_array_equal(cupy.asnumpy(selected.val), expected)
    np.testing.assert_array_equal(selected.iid, source.iid[rows])
    np.testing.assert_array_equal(selected.sid, source.sid[columns])
    assert selected.val.flags["C_CONTIGUOUS" if order == "C" else "F_CONTIGUOUS"]

    # Independent centered float64 reference, including a missing value and
    # constant column. Allow 32 eps for the bounded fractional input reductions.
    centered = expected.astype(np.float64) - np.nanmean(expected.astype(np.float64), axis=0)
    deviation = np.sqrt(np.nanmean(centered ** 2, axis=0))
    deviation[deviation == 0] = np.inf
    expected = centered / deviation
    expected[np.isnan(expected)] = 0
    selected.standardize(Unit())
    tolerance = 32 * np.finfo(dtype).eps
    np.testing.assert_allclose(cupy.asnumpy(selected.val), expected, rtol=tolerance, atol=tolerance)
    np.testing.assert_array_equal(cupy.asnumpy(source.val), values)


@pytest.mark.parametrize("rows,columns", [([], [1, 0]), ([1, 0], []), ([], [])])
def test_empty_subsets_have_cupy_storage(cupy, rows, columns):

    source = data(cupy.arange(6, dtype=cupy.float32).reshape(3, 2), cupy)
    selected = source[rows, columns].read(dtype="float32", xp=cupy)
    assert selected.val.shape == (len(rows), len(columns))


def test_three_dimensional_subset_remains_cupy(cupy):
    from pysnptools.util import sub_matrix

    values = np.arange(24, dtype=np.float32).reshape(4, 3, 2)
    result = sub_matrix(cupy.array(values), [3, 0], [2, 1], dtype="float32")
    np.testing.assert_array_equal(cupy.asnumpy(result), values[np.ix_([3, 0], [2, 1])])


def test_subset_narrowing_warns_and_preserves_values(cupy):

    values = np.array([[.25, np.nan], [1.25, 1], [1.75, 2]], dtype=np.float64)
    source = data(cupy.array(values), cupy)
    with pytest.warns(UserWarning, match="loss of information"):
        result = source[[2, 0], [1, 0]].read(dtype="float32", xp=cupy)
    np.testing.assert_array_equal(cupy.asnumpy(result.val), values[np.ix_([2, 0], [1, 0])].astype(np.float32))


def test_trained_cupy_data_uses_cupy_application(cupy):
    source = data(np.array([[0., 1.], [1., 1.], [2., 1.]]), cupy)
    _, trained = source.standardize(return_trained=True)
    target = data(np.array([[0., 1.], [1., 1.], [2., 1.]]), cupy)
    target.standardize(trained)
    np.testing.assert_array_equal(cupy.asnumpy(target.val), cupy.asnumpy(source.val))


def test_read_narrowing_warns_and_view_permissions_are_preserved(cupy):

    source = data(cupy.array([[0.25, 1.], [1.25, 1.], [1.75, 1.]]), cupy)
    shared = source.read(dtype="float64", order="A", view_ok=True, xp=cupy)
    assert shared.val.data.ptr == source.val.data.ptr
    with pytest.warns(UserWarning, match="loss of information"):
        narrowed = source.read(dtype="float32", xp=cupy)
    assert narrowed.val.dtype == np.float32
    np.testing.assert_array_equal(cupy.asnumpy(narrowed.val), cupy.asnumpy(source.val).astype(np.float32))


def test_cupy_allocation_failure_propagates(cupy, monkeypatch):
    def fail(*args, **kwargs):
        raise cupy.cuda.memory.OutOfMemoryError(24, 0)

    caller_values = cupy.zeros((3, 2))
    with monkeypatch.context() as context:
        context.setattr(cupy, "array", fail)
        with pytest.raises(cupy.cuda.memory.OutOfMemoryError):
            data(np.zeros((3, 2)), cupy)
        assert data(caller_values, cupy).val is caller_values


def test_other_ranks_and_int8_remain_supported_storage(cupy):
    from pysnptools.pstreader import PstData

    values = np.arange(12, dtype=np.float32).reshape(2, 3, 2)
    with patch.dict("os.environ", {"ARRAY_MODULE": "cupy"}):
        target = PstData(row=["0", "1"], col=["s0", "s1", "s2"], val=values)
    np.testing.assert_array_equal(cupy.asnumpy(target.val), values)
    integer = SnpData(iid=[["f", "0"], ["f", "1"]], sid=["s"],
                      val=cupy.ones((2, 1), dtype=cupy.int8), xp=cupy, _require_float32_64=False)
    with pytest.raises(ValueError, match="float32 or float64"):
        integer.standardize()


def test_empty_containers_use_floating_cupy_storage(cupy):
    source = SnpData(iid=[], sid=["s0", "s1"], val=None, xp=cupy)
    kernel = KernelData(iid=[], val=None, xp=cupy)
    for target, shape in ((source, (0, 2)), (kernel, (0, 0))):
        assert isinstance(target.val, cupy.ndarray) and target.val.dtype == np.float64
        assert target.val.shape == shape
