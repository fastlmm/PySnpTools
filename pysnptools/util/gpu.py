"""GPU availability and precision notices for the optional CuPy backend."""

from pysnptools.util._gpu_warnings import GpuFloat64Warning

__all__ = ["GpuFloat64Warning", "has_gpu_support"]


def has_gpu_support():
    """Return whether CuPy can allocate and synchronize on the current device.

    Optional CuPy import failures and CUDA runtime errors report unavailable
    support. Other errors propagate rather than hiding a broken installation.
    No bed-reader GPU build is required.

    >>> isinstance(has_gpu_support(), bool)
    True
    """
    try:
        import cupy
    except (ImportError, OSError):
        return False
    try:
        if cupy.cuda.runtime.getDeviceCount() == 0:
            return False
        cupy.empty(1, dtype="float32")
        cupy.cuda.get_current_stream().synchronize()
    except cupy.cuda.runtime.CUDARuntimeError:
        return False
    return True
