"""Allocate and copy arrays with the selected NumPy or CuPy backend."""

import warnings

import numpy as np
from pysnptools.util._gpu_warnings import warn_gpu_float64


def empty(shape, dtype, order, xp):
    dtype = np.dtype(dtype)
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
    warn_gpu_float64(dtype, xp)
    return (xp.array(source, dtype=dtype, order=order) if copy
            else xp.asarray(source, dtype=dtype, order=order))


def zeros(shape, dtype, order, xp):
    if xp is np:
        return np.zeros(shape, dtype=dtype, order=order)
    result = empty(shape, dtype, order, xp)
    result.fill(0)
    return result
