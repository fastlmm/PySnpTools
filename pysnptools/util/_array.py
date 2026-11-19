"""Allocate new GPU matrices through bed-reader when its CUDA allocator applies."""

import warnings

import numpy as np
from pysnptools.util._gpu_warnings import warn_gpu_float64


def _allocator(shape, dtype, xp):
    if xp is not np:
        from bed_reader import bed_reader as native
        from bed_reader._gpu import empty_gpu

        # CPU-only bed-reader builds lack the native CubeCL allocator.
        if (hasattr(native, "allocate_gpu")
                and len(shape) == 2 and max(*shape, shape[0] * shape[1]) <= 2**32 - 1
                and dtype in (np.dtype("int8"), np.dtype("float32"), np.dtype("float64"))
                and xp.cuda.Device().id == 0):
            return empty_gpu
    return None


def empty(shape, dtype, order, xp):
    dtype = np.dtype(dtype)
    allocator = _allocator(shape, dtype, xp)
    if allocator is not None:
        return allocator(shape, dtype=dtype, order=order)
    # CuPy-only builds, other devices, dtypes, and ranks still use GPU memory.
    warn_gpu_float64(dtype, xp)
    return xp.empty(shape, dtype=dtype, order=order)


def asarray(values, dtype, order, xp, copy=False):
    dtype = np.dtype(dtype)
    if xp is np:
        return np.array(values, dtype=dtype, order=order) if copy else np.asarray(values, dtype=dtype, order=order)
    if isinstance(values, xp.ndarray):
        if values.dtype == np.float64 and dtype == np.float32:
            warnings.warn("Converting float64 to float32 can cause loss of information", stacklevel=2)
        fits_order = (order in ("A", "K") or
                      (order == "C" and values.flags.c_contiguous) or
                      (order == "F" and values.flags.f_contiguous))
        if not copy and values.dtype == dtype and fits_order:
            return values
        source = values
    else:
        source = np.asarray(values, dtype=dtype)
    if source.ndim != 2:
        warn_gpu_float64(dtype, xp)
        return xp.array(source, dtype=dtype, order=order)
    if order in ("A", "K"):
        order = "F" if abs(source.strides[0]) < abs(source.strides[1]) else "C"
    if _allocator(source.shape, dtype, xp) is None:
        warn_gpu_float64(dtype, xp)
        return xp.array(source, dtype=dtype, order=order)
    result = empty(source.shape, dtype, order, xp)
    if isinstance(source, np.ndarray):
        host = np.asarray(source, dtype=dtype, order=order)
        result.set(host)
        # The temporary host array can die on return; finish the copy first.
        xp.cuda.get_current_stream().synchronize()
    else:
        xp.copyto(result, source, casting="unsafe")
    return result


def zeros(shape, dtype, order, xp):
    if xp is np:
        return np.zeros(shape, dtype=dtype, order=order)
    result = empty(shape, dtype, order, xp)
    result.fill(0)
    return result


def native_owner(values):
    # Ordinary CuPy arrays retain the CuPy path.
    # CubeCL's standardizer also indexes four reduction words per dimension.
    if max(values.shape) > (2**32 - 1) // 4:
        return None
    from bed_reader._gpu import _native_owner

    return _native_owner(values)
