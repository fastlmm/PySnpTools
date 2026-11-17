"""GPU precision notices, including compatibility with older bed-reader builds."""

import warnings

import numpy as np


def warn_gpu_float64(dtype, xp):
    if xp is np or np.dtype(dtype) != np.dtype("float64"):
        return
    import bed_reader

    warn = getattr(bed_reader, "_warn_gpu_float64", None)
    if warn is not None:
        warn(dtype)
    else:
        # CuPy-only use remains supported with older published bed-reader builds.
        warnings.warn(
            "float64 arithmetic can be much slower than float32 on consumer GPUs. Float64 GPU arrays "
            "also use twice the storage. Consider float32 where its precision is sufficient; "
            "the requested dtype is unchanged.",
            getattr(bed_reader, "GpuFloat64Warning", UserWarning),
            stacklevel=2,
        )
