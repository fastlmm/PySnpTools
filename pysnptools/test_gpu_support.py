"""GPU detection belongs to PySnpTools and does not depend on bed-reader."""

import sys
from types import SimpleNamespace

import pytest

from pysnptools.util.gpu import has_gpu_support


def test_no_cupy(monkeypatch):
    monkeypatch.setitem(sys.modules, "cupy", None)
    assert not has_gpu_support()


@pytest.mark.parametrize("device_count", [0, 1])
def test_device_probe(monkeypatch, device_count):
    operations = []
    runtime = SimpleNamespace(getDeviceCount=lambda: device_count, CUDARuntimeError=RuntimeError)
    stream = SimpleNamespace(synchronize=lambda: operations.append("sync"))
    cupy = SimpleNamespace(cuda=SimpleNamespace(runtime=runtime, get_current_stream=lambda: stream),
                           empty=lambda *args, **kwargs: operations.append("allocate"))
    monkeypatch.setitem(sys.modules, "cupy", cupy)
    assert has_gpu_support() == bool(device_count)
    assert operations == (["allocate", "sync"] if device_count else [])


def test_runtime_failure_and_unexpected_error(monkeypatch):
    class DeviceError(RuntimeError):
        pass

    runtime = SimpleNamespace(CUDARuntimeError=DeviceError, getDeviceCount=lambda: 1)
    def fail(*args, **kwargs):
        raise DeviceError("device unavailable")
    cupy = SimpleNamespace(cuda=SimpleNamespace(runtime=runtime), empty=fail)
    monkeypatch.setitem(sys.modules, "cupy", cupy)
    assert not has_gpu_support()
    def unexpected(*args, **kwargs):
        raise ValueError("unexpected failure")
    cupy.empty = unexpected
    with pytest.raises(ValueError, match="unexpected failure"):
        has_gpu_support()
