"""Propagate precision warnings for native and ordinary CuPy allocations/math."""

import warnings

import bed_reader
import numpy as np
import pytest

from pysnptools.snpreader import SnpData
from pysnptools.standardizer import Identity

WarningType = bed_reader.GpuFloat64Warning


def data(values, xp):
    return SnpData(iid=[["f", str(sample)] for sample in range(len(values))],
                   sid=[str(marker) for marker in range(values.shape[1])], val=values, xp=xp)


@pytest.fixture
def cupy():
    return pytest.importorskip("cupy")


def test_allocations_warn_for_both_native_and_cupy(cupy, monkeypatch):
    from bed_reader import bed_reader as native

    for native_enabled in (True, False):
        with monkeypatch.context() as patch:
            if not native_enabled:
                patch.delattr(native, "allocate_gpu", raising=False)
            with pytest.warns(WarningType, match="consumer GPUs") as caught:
                result = data(np.ones((3, 2)), cupy)
            assert len(caught) == 1 and result.val.dtype == np.float64


def test_existing_arrays_warn_on_math_and_float32_reductions_do_not(cupy):
    values = cupy.array([[0., 1.], [1., 1.], [2., 1.]])
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        result = data(values, cupy)
    assert result.val is values and not caught  # no allocation or conversion
    with pytest.warns(WarningType, match="consumer GPUs"):
        result.standardize()
    assert result.val.dtype == np.float64
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        data(np.ones((3, 2), dtype=np.float32), cupy).standardize()
        data(np.ones((3, 2)), np).standardize()
    assert not any(issubclass(item.category, WarningType) for item in caught)


def test_read_kernel_and_default_filter(cupy):
    source = data(np.array([[0., 1.], [1., 1.], [2., 1.]]), np)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("default", WarningType)
        result = source.read(dtype="float64", xp=cupy).standardize()
    notices = [item for item in caught if issubclass(item.category, WarningType)]
    assert len(notices) == 1
    assert result.val.dtype == np.float64
    with pytest.warns(WarningType):
        kernel = source.read_kernel(Identity(), dtype="float64", xp=cupy)
    assert kernel.val.dtype == np.float64
