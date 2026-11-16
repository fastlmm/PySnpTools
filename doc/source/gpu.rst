GPU arrays, subsetting, and standardization
===========================================

The GPU changes described here are under development on this branch. They are
not a promise that the latest published PySnpTools contains these changes.
CuPy is optional and must be installed separately for your GPU platform.
For the concurrent explicit reader API and native GPU BED path, see
:doc:`gpu-api-preview`.

Construct arrays and select a backend
-------------------------------------

``SnpData`` and ``KernelData`` accept ``xp`` as a module or module name.
An explicit choice overrides ``ARRAY_MODULE``; without either, NumPy is used.
The numerical values live in the selected module; IDs and position metadata
remain NumPy arrays on the CPU. Moving an existing float32 or float64 array
between NumPy and CuPy preserves its dtype. Lists and integer arrays normally
become float64. A same-module floating-point input can share storage, so copy it
first when later in-place changes must not affect the caller.

This complete example needs CuPy and a working GPU, but no BED file or native
GPU bed-reader build::

    import cupy as cp
    import numpy as np
    from pysnptools.snpreader import SnpData
    from pysnptools.util import asnumpy

    assert cp.cuda.runtime.getDeviceCount() > 0
    values = np.array([[0, 2], [1, 1], [2, 0]], dtype=np.float32)
    data = SnpData(
        iid=[["family", str(sample)] for sample in range(3)],
        sid=["snp0", "snp1"], val=values, xp=cp
    )
    assert isinstance(data.val, cp.ndarray) and data.val.dtype == np.float32
    host = SnpData(iid=data.iid, sid=data.sid, pos=data.pos,
                   val=data.val, xp=np)
    np.testing.assert_array_equal(host.val, values)
    assert host.val.dtype == np.float32
    np.testing.assert_array_equal(asnumpy(data.val), values)

``asnumpy`` transfers a CuPy array to the CPU and returns a NumPy array unchanged.
``array_module("cupy")`` can log a warning and select NumPy if the CuPy import
fails. Importing CuPy explicitly, as above, avoids accidentally mistaking this
fallback for a successful GPU run. Driver/runtime failures are separate errors.

Subset without moving values to the CPU
---------------------------------------

``sub_matrix`` dispatches from the actual input array, independently of
``ARRAY_MODULE``. It supports contiguous float32/float64 arrays with two or three
dimensions and keeps the result on the input device. Row and column lists select
a Cartesian product, preserving order and duplicates. On the CuPy path, negative
or out-of-range indices raise ``IndexError`` instead of CuPy's wrapping behavior::

    from pysnptools.util import sub_matrix

    selected = sub_matrix(data.val, [2, 0, 2], [1, 0],
                          dtype=data.val.dtype, order="F")
    assert isinstance(selected, cp.ndarray)
    assert selected.dtype == np.float32 and selected.flags.f_contiguous
    np.testing.assert_array_equal(cp.asnumpy(selected), values[[2, 0, 2]][:, [1, 0]])

The default result dtype is float64; pass ``dtype`` explicitly to keep float32.
Narrowing float64 to float32 emits a precision warning. These low-level index
rules differ from reader slicing, which accepts slices and Boolean masks.

On this committed branch, reader ``read`` methods do not yet accept ``xp``.
They consult ``ARRAY_MODULE`` when constructing the returned data object, even
if the original object was constructed with an explicit ``xp``. To keep an
in-memory subset on the GPU, set the environment for the read and request the
dtype explicitly. The context below restores the prior environment::

    from unittest.mock import patch

    with patch.dict("os.environ", {"ARRAY_MODULE": "cupy"}):
        subset = data[[2, 0], :].read(dtype="float32")
    assert isinstance(subset.val, cp.ndarray) and subset.val.dtype == np.float32
    np.testing.assert_array_equal(cp.asnumpy(subset.val), values[[2, 0]])

``ARRAY_MODULE`` is process-wide. Avoid changing it concurrently with other
workflows. An in-memory CuPy array does not imply that all readers, writers, or
utility functions accept CuPy inputs.

Standardize and reuse training statistics
-----------------------------------------

``Unit`` and ``Beta`` standardize numerical values in place. On this branch,
keep ``ARRAY_MODULE`` consistent with the data device when calling Unit: it
consults the environment and can fail when asked to convert CuPy values to
NumPy implicitly. Unit centers each SNP and divides by its observed standard deviation;
missing values become zero. A SNP with no observations raises ``ValueError``.
A constant training SNP has an infinite stored standard deviation so a trained
standardizer zeros that SNP in later data too. Beta centers and weights SNPs
using the beta density at the minor allele frequency; it is not Unit scaling.
Its SciPy beta-density calculation transfers per-SNP values through the CPU.

Using the ``data`` from the first example::

    from pysnptools.standardizer import Unit, Beta

    with patch.dict("os.environ", {"ARRAY_MODULE": "cupy"}):
        standardized, trained = data.standardize(Unit(), return_trained=True)
        assert standardized is data
        assert isinstance(data.val, cp.ndarray) and data.val.dtype == np.float32
        fresh = SnpData(iid=data.iid, sid=data.sid,
                        val=cp.asarray(values), xp=cp)
        fresh.standardize(trained)
        np.testing.assert_allclose(cp.asnumpy(fresh.val), cp.asnumpy(data.val),
                                   rtol=1e-6, atol=1e-6)
        weighted = SnpData(iid=data.iid, sid=data.sid,
                           val=cp.asarray(values), xp=cp).standardize(Beta(1, 25))
        assert isinstance(weighted.val, cp.ndarray)

The tolerance here is for this small, well-scaled float32 example, not a general
numerical guarantee. Applying trained statistics requires matching SNP IDs and
order; reuse statistics fitted on training data when transforming test data.
The training call has already changed ``data.val`` in place.

Read files and construct a kernel
---------------------------------

On this branch, selecting CuPy for a PySnpTools BED read transfers values after
the existing file read; it does not itself request native GPU BED decoding.
FaST-LMM's GPU BED integration is a separate path. For an existing BED dataset::

    from pysnptools.snpreader import Bed

    with patch.dict("os.environ", {"ARRAY_MODULE": "cupy"}):
        bed = Bed("cohort.bed", count_A1=True)
        block = bed[:100, :200].read(dtype="float32").standardize(Unit())
        kernel = bed.read_kernel(Unit(), block_size=1000, dtype="float32")
    assert isinstance(block.val, cp.ndarray)
    assert isinstance(kernel.val, cp.ndarray)

The BED file must have accompanying BIM and FAM files. Keep ``count_A1`` explicit.
``read_kernel(Unit())`` forms the product of standardized SNPs and their
transpose; use kernel standardization separately if a unit average diagonal is
needed. ``block_size`` is a suggested SNP block size and does not cap total GPU
memory. The result is still a dense samples-by-samples matrix. In the committed
implementation, the small-SNP/low-rank path can read all SNPs despite a requested
block size, and kernel accumulation can use float64 before the requested output
cast. Do not infer an all-float32 computation or a memory guarantee from the
returned dtype.
