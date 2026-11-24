"""Refuse implicit NumPy/CuPy transfers of existing arrays.

``xp=None`` always means ARRAY_MODULE, then NumPy. When an existing NumPy or
CuPy array would have to move to the other module to honor that default, the
caller must say so with an explicit ``xp=``; PySnpTools raises ValueError
rather than silently copying between the CPU and the GPU. Lists, ``val=None``,
and values decoded from files have no device and simply follow the default.

Constructors and ``val`` setters check their argument directly (:func:`resolve`,
:func:`check_assignment`). Reads check in :meth:`PstData._read`, where every
in-memory source is read, but its legacy signature carries no ``xp``. So each
public read entry point (:func:`read_entry`) records the request in a context
variable, and the outermost entry decides: internal reads made on its behalf,
which pass a resolved module down, keep the caller's choice.
"""

import functools
import inspect
import os
import sys
from contextvars import ContextVar

import numpy as np
import pysnptools.util as pstutil

# (target module, whether the caller passed xp explicitly), or None outside any read.
_read_request = ContextVar("pysnptools_read_request", default=None)


def module_of(value):
    """Return numpy or cupy for a NumPy or CuPy array, else None (lists, scalars, None)."""
    if isinstance(value, np.ndarray):
        return np
    cupy = sys.modules.get("cupy")  # a CuPy array exists only if CuPy was imported
    if cupy is not None and isinstance(value, cupy.ndarray):
        return cupy
    return None


def _mismatch_error(source, target, what):
    setting = os.environ.get("ARRAY_MODULE")
    setting = "ARRAY_MODULE unset" if setting is None else f"ARRAY_MODULE={setting}"
    return ValueError(
        f"{what} is a {source.__name__} array, but xp=None selects {target.__name__} "
        f"({setting}). Pass xp='{target.__name__}' to copy it, or "
        f"xp='{source.__name__}' to keep it on its current device.")


def resolve(xp, value, what):
    """Resolve ``xp`` for an existing value; see the module docstring."""
    target = pstutil.array_module(xp)
    if xp is None:
        source = module_of(value)
        if source is not None and source is not target:
            raise _mismatch_error(source, target, what)
    return target


def check_assignment(value, xp, what):
    """Require an array assigned to ``val`` to use the object's module, ``xp``.

    A setter has no ``xp`` argument, so it never moves an array between modules."""
    source = module_of(value)
    if source is not None and source is not xp:
        raise ValueError(
            f"{what} is a {xp.__name__} array; cannot assign a {source.__name__} array. "
            f"Convert it first, or construct a new object with xp='{source.__name__}'.")


def read_entry(function):
    """Record a public read's ``xp`` request for :func:`check_source`."""
    signature = inspect.signature(function)

    @functools.wraps(function)
    def wrapper(*args, **kwargs):
        if _read_request.get() is not None:  # nested: the outermost request decides
            return function(*args, **kwargs)
        xp = signature.bind(*args, **kwargs).arguments.get("xp")
        token = _read_request.set((pstutil.array_module(xp), xp is not None))
        try:
            return function(*args, **kwargs)
        finally:
            _read_request.reset(token)

    return wrapper


def check_source(val, what):
    """Raise if an implicit (``xp=None``) read would move in-memory ``val`` to the other module."""
    request = _read_request.get()
    if request is None:
        return
    target, explicit = request
    if explicit:
        return
    source = module_of(val)
    if source is not None and source is not target:
        raise _mismatch_error(source, target, what)
