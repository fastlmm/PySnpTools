"""Scoped thread limits shared by local multiprocessing runners."""
from contextlib import contextmanager
from functools import wraps
from inspect import signature
from operator import index
from unittest.mock import patch
import warnings

from numpy import bool_
from threadpoolctl import threadpool_limits


def environment(num_threads):
    """Return startup settings for BLAS, OpenMP, and bed-reader workers."""
    if num_threads is None:
        return {}
    if isinstance(num_threads, (bool, bool_)):
        raise TypeError("num_threads must be a positive integer or None")
    num_threads = index(num_threads)
    if num_threads <= 0:
        raise ValueError("num_threads must be a positive integer or None")
    return {
        name: str(num_threads)
        for name in (
            "MKL_NUM_THREADS",
            "OPENBLAS_NUM_THREADS",
            "OMP_NUM_THREADS",
            "PST_NUM_THREADS",
            "BLIS_NUM_THREADS",
        )
    }


@contextmanager
def limits(num_threads):
    """Limit loaded BLAS libraries and restore the caller's settings on exit."""
    if num_threads is None:
        yield
        return
    settings = environment(num_threads)
    with patch.dict("os.environ", settings), threadpool_limits(
        limits=index(num_threads), user_api="blas"
    ):
        yield


def renamed_num_threads(function):
    """Accept the deprecated keyword without changing the canonical signature."""
    parameters = signature(function)

    @wraps(function)
    def wrapped(*args, **kwargs):
        if "mkl_num_threads" in kwargs:
            legacy = kwargs.pop("mkl_num_threads")
            # bind_partial also detects an explicitly supplied positional None.
            # Checking the value alone cannot distinguish None from omission.
            if "num_threads" in parameters.bind_partial(*args, **kwargs).arguments:
                raise TypeError("Specify only num_threads; mkl_num_threads is its deprecated alias")
            warnings.warn(
                "mkl_num_threads is deprecated; use num_threads instead",
                DeprecationWarning,
                stacklevel=2,
            )
            kwargs["num_threads"] = legacy
        return function(*args, **kwargs)

    return wrapped


def _legacy_get(runner):
    warnings.warn("mkl_num_threads is deprecated; use num_threads instead", DeprecationWarning, stacklevel=2)
    return runner.num_threads


def _legacy_set(runner, num_threads):
    warnings.warn("mkl_num_threads is deprecated; use num_threads instead", DeprecationWarning, stacklevel=2)
    runner.num_threads = num_threads


# Keep the previously writable public attribute synchronized with its new name.
legacy_num_threads = property(_legacy_get, _legacy_set, doc="Deprecated alias of num_threads.")



def _current_threads(runner):
    return runner.__dict__["mkl_num_threads"]


def _set_current_threads(runner, num_threads):
    runner.__dict__["mkl_num_threads"] = num_threads


# Keep the historical pickle key while exposing the canonical attribute.
# Default object pickling also preserves slots added by downstream subclasses.
canonical_num_threads = property(_current_threads, _set_current_threads, doc="CPU thread limit, or None.")
