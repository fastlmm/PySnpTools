"""Precision defaults for PySnpTools reads, kernels, and allocations."""
import os

import numpy as np


def array_dtype(dtype=None):
    """Resolve an explicit dtype, ARRAY_DTYPE, then float64.

    ARRAY_DTYPE accepts float32 or float64. Explicit values override it;
    individual APIs still validate supported types. Existing floating arrays
    and standardizers keep their precision. NumPy/CuPy's own defaults are
    unaffected. See :doc:`/gpu` for a complete workflow example.

    >>> from unittest.mock import patch
    >>> from pysnptools.util import array_dtype
    >>> with patch.dict("os.environ", {"ARRAY_DTYPE": "float32"}):
    ...     print(array_dtype(), array_dtype("float64"))
    float32 float64
    """
    # Same semantics as bed-reader's ARRAY_DTYPE policy, kept here so PySnpTools
    # does not depend on a private bed-reader module.
    if dtype is not None:
        return np.dtype(dtype)
    configured = os.environ.get("ARRAY_DTYPE")
    if configured is None:
        return np.dtype("float64")
    if configured not in ("float32", "float64"):
        raise ValueError(f"ARRAY_DTYPE must be 'float32' or 'float64', not {configured!r}")
    return np.dtype(configured)
