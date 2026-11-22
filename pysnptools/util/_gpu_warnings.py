"""GPU precision warnings shared by PySnpTools and FaST-LMM."""

import sys
import warnings

import numpy as np


class GpuFloat64Warning(UserWarning):
    """Float64 GPU storage or computation may have a substantial performance cost.

    Float64 remains supported and is never silently converted. See :doc:`gpu`
    for precision choices and a float32 example.
    """


def warn_gpu_float64(dtype, xp):
    if xp is np or np.dtype(dtype) != np.dtype("float64"):
        return
    # Attribute nested reader/workflow warnings to the user's call. Python's
    # normal warning filter then shows this message once per user call site.
    frame = sys._getframe(1)
    stacklevel = 2
    while frame is not None:
        modules = frame.f_globals.get("__name__", "").split(".")
        # pandas can consume our lazy map/reduce results; it is an intermediate
        # frame, not the user's precision choice or allocation call site.
        if modules[0] not in {"bed_reader", "pysnptools", "fastlmm", "pandas"} or any(
            part in {"test", "tests"} or part.startswith("test_") for part in modules[1:]
        ):
            break
        stacklevel += 1
        frame = frame.f_back
    del frame
    warnings.warn(
        "float64 arithmetic can be much slower than float32 on consumer GPUs. Float64 GPU arrays "
        "also use twice the storage. Consider float32 where its precision is "
        "sufficient; the requested dtype is unchanged.",
        GpuFloat64Warning,
        stacklevel=stacklevel,
    )
