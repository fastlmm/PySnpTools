"""Precision defaults shared with bed-reader."""
from bed_reader._array_module import array_dtype as _array_dtype


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
    return _array_dtype(dtype)
