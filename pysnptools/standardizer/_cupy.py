"""Lazily compiled CuPy kernels; importing PySnpTools does not require CuPy."""
from functools import lru_cache


@lru_cache(maxsize=1)
def _sum_squares_kernel():
    import cupy as cp

    # Normalization must propagate NaNs, unlike missing-aware SNP variance.
    return cp.ReductionKernel(
        "T value", "float64 total", "(double)value * (double)value",
        "a + b", "total = a", "0", "normalization_sum_squares")


@lru_cache(maxsize=1)
def kernels():
    import cupy as cp

    count = cp.ReductionKernel(
        "T value", "int64 count", "value == value ? 1LL : 0LL",
        "a + b", "count = a", "0", "observed_count")
    total = cp.ReductionKernel(
        "T value", "float64 total", "value == value ? (double)value : 0.0",
        "a + b", "total = a", "0", "observed_sum")
    variance = cp.ReductionKernel(
        "T value, float64 mean", "float64 total",
        "value == value ? ((double)value - mean) * ((double)value - mean) : 0.0",
        "a + b", "total = a", "0", "centered_sum")
    apply = cp.ElementwiseKernel(
        "T mean, T factor", "T value",
        "value = value == value && factor != (T)0 ? (value - mean) * factor : (T)0",
        "apply_standardize")
    return count, total, variance, apply
