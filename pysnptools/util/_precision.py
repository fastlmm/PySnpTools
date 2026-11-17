"""Precision defaults shared with bed-reader without breaking older CPU installs."""
import os

import bed_reader
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
    resolver = getattr(bed_reader, "array_dtype", None)
    if resolver is not None:
        return resolver(dtype)
    # Supported older bed-reader releases lack this new optional policy. Keep
    # their legacy defaults/explicit arguments working; never ignore a requested
    # environment policy that the installed dependency cannot implement.
    if dtype is None and "ARRAY_DTYPE" in os.environ:
        raise RuntimeError("ARRAY_DTYPE requires bed-reader with array_dtype support")
    return np.dtype("float64" if dtype is None else dtype)
