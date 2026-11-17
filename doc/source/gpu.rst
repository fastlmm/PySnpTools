GPU arrays, readers, and standardization
========================================

This guide describes the current development implementation, including the
reader GPU API. Matching development versions of PySnpTools and bed-reader are
required; installing the published packages does not necessarily provide these
features. NumPy runs on the CPU. CuPy support here is experimental CUDA support;
Python HIP support is not implemented.

Backend, dtype, and order
-------------------------

An explicit ``xp=np``, ``xp=cp``, ``xp="numpy"``, or ``xp="cupy"`` overrides
``ARRAY_MODULE``. Only ``xp=None`` consults that environment variable; absent
both, NumPy is used. Missing CuPy raises ``ImportError``. An unknown module name
raises ``ValueError``; a successful CuPy import does not establish driver or
hardware availability. There is no automatic CPU retry after a GPU error.

.. list-table:: Reader and constructor options
   :header-rows: 1
   :widths: 24 22 54

   * - API or option
     - Default
     - Behavior
   * - ``SnpData(..., xp=...)``, ``KernelData(..., xp=...)``
     - ``xp=None``
     - Selects the module for values; metadata remains NumPy. Float32/float64
       inputs keep their dtype across transfers. Lists and integers normally
       become float64. Same-module floating arrays can share storage.
   * - ``SnpReader.read``, ``KernelReader.read``
     - ``xp=None``, ``dtype=np.float64``, ``order="F"``
     - Return ``SnpData`` or ``KernelData``. Every read chooses its backend
       independently of the constructor and its requested output dtype.
   * - ``SnpReader.read_kernel``
     - ``standardizer=None``, ``block_size=None``, ``order="A"``
     - A standardizer is required despite the ``None`` signature default;
       omitting it raises ``AssertionError``. Returns ``KernelData``.
       Output dtype defaults to float64; ``xp=None`` uses the environment.
   * - ``order``
     - As above
     - ``F`` varies the first index fastest; ``C`` varies the last index
       fastest. ``A`` allows the reader to choose. BED maps ``A`` to ``F``;
       blocked kernel accumulation maps ``A`` to ``C``.
   * - ``view_ok``
     - ``False``
     - Allows sharing when practical for in-memory reads; it does not guarantee
       a view. BED reads and kernel products allocate output regardless.
   * - ``num_threads``
     - ``None``
     - CPU thread count, also controlled by ``PST_NUM_THREADS``, ``NUM_THREADS``,
       then ``MKL_NUM_THREADS``. Ignored by native GPU BED decoding and CuPy math.
   * - ``force_python_only`` on reads
     - ``False``
     - Selects a CPU reader's Python path where supported. Native GPU BED reads
       reject ``True`` with ``ValueError``; it does not select a CPU backend.

Ordinary numerical reads and standardization support float32 and float64.
BED also supports ``dtype="int8", _require_float32_64=False`` (missing = -127);
int8 data must be converted to floating point before standardization. Float16
is not a supported numerical read dtype. In-memory constructors preserve array
dtypes, not an arbitrary requested layout; use ``read(order=..., dtype=..., xp=...)``
for explicit output requirements. Dense kernels have shape samples by samples.

The following complete example checks both directions of environment overrides,
actual array types, dtype, order, and a duplicate-preserving subset. It needs
CuPy and a working GPU, but no native BED GPU build or data download::

    import cupy as cp
    import numpy as np
    from unittest.mock import patch
    from pysnptools.snpreader import SnpData
    from pysnptools.util import asnumpy, sub_matrix

    assert cp.cuda.runtime.getDeviceCount() > 0
    cp.zeros(1).sum().item()
    values = np.array([[0, 2], [1, 1], [2, 0]], dtype=np.float32)
    data = SnpData(iid=[["family", str(sample)] for sample in range(3)],
                   sid=["snp0", "snp1"], val=values.copy(), xp=np)
    with patch.dict("os.environ", {"ARRAY_MODULE": "numpy"}):
        gpu = data.read(dtype="float32", order="C", xp=cp)
        subset = gpu[[2, 0, 2], :].read(dtype="float32", order="F", xp=cp)
        subset.standardize()  # uses CuPy despite ARRAY_MODULE=numpy
    assert isinstance(gpu.val, cp.ndarray) and gpu.val.flags.c_contiguous
    assert isinstance(subset.val, cp.ndarray) and subset.val.flags.f_contiguous
    assert subset.val.dtype == np.float32
    with patch.dict("os.environ", {"ARRAY_MODULE": "cupy"}):
        host = gpu.read(dtype="float32", order="F", xp=np)
        host.standardize()  # stays on NumPy
    assert isinstance(host.val, np.ndarray) and host.val.flags.f_contiguous
    assert host.val.dtype == np.float32
    selected = sub_matrix(gpu.val, [2, 0, 2], [1, 0], dtype=gpu.val.dtype)
    assert isinstance(selected, cp.ndarray) and selected.dtype == np.float32
    np.testing.assert_array_equal(asnumpy(selected), values[[2, 0, 2]][:, [1, 0]])

``sub_matrix`` chooses from the actual input array, independently of the
environment. It accepts contiguous 2-D/3-D float32/float64 arrays and
one-dimensional integer index sequences, preserving order and duplicates.
It rejects negative, out-of-range, floating-point, and Boolean indices with
``IndexError``. Reader slicing also accepts slices and Boolean masks.
``sub_matrix`` defaults to float64: pass ``dtype=values.dtype`` to keep float32.
Narrowing float64 to float32 warns of precision loss. ``asnumpy`` copies a CuPy
array to the host and returns a NumPy array unchanged.

Standardization and trained statistics
--------------------------------------

Unit and Beta modify values in place on their existing device; they do not
consult ``ARRAY_MODULE``. Values must be contiguous floating-point arrays.
``return_trained=True`` returns ``(data, trained_standardizer)``; otherwise the
return is the same data object. Training through Unit/Beta returns NumPy
``trained.stats`` with shape ``(sid_count, 2)``: mean and population standard
deviation (``ddof=0``), in the input floating dtype. Even GPU training transfers
these small per-SNP statistics to the CPU.

.. list-table:: Standardizer behavior
   :header-rows: 1
   :widths: 23 77

   * - Standardizer
     - Behavior
   * - ``Unit()``
     - Subtract the observed mean and divide by the observed population
       standard deviation. Impute NaNs with zero afterward. The final column
       variance is generally below one when observations are missing; a
       constant column becomes zero.
   * - ``Beta(a, b)``
     - Center and multiply by the beta density at the minor allele frequency
       ``min(mean / 2, 1 - mean / 2)``; do not divide by standard deviation.
       Native and CuPy paths validate means in [0, 2] and beta parameters.
       CuPy reduces on the GPU but evaluates the per-SNP density with SciPy
       on the CPU, then transfers weights back.
   * - ``UnitTrained``
     - Reuse fitted moments. Can select/reorder known SNP IDs automatically;
       an unknown ID raises ``KeyError``. Unique training IDs avoid ambiguity.
   * - ``BetaTrained``
     - Reuse fitted moments and beta parameters. Requires exactly the same
       SNP IDs in the same order; a mismatch raises ``AssertionError``.
       Kernel construction checks the full SNP IDs/order, then reuses the
       corresponding fitted moments for each block without refitting them.

Fitting a wholly missing SNP raises ``ValueError``. Applying already-trained
moments does not refit or require an observation: an all-NaN test column becomes
zero. A constant training SNP stores infinite standard deviation and remains
zero in test data even if its new values vary. Unit accepts general finite
values, including negative values; Beta is for allele counts/dosages with valid
allele-frequency means. CPU and GPU rounding need not be identical, especially
for large offsets and small variance.

For Unit/Beta standardizer methods, ``block_size`` is deprecated, ignored, and
warns when supplied directly. ``SnpData.standardize(block_size=...)`` currently
ignores it without forwarding a warning. ``force_python_only=True`` selects the
NumPy implementation on CPU but is ignored for CuPy arrays; it never transfers
them to the CPU. ``num_threads`` applies to native CPU standardization and is
ignored for CuPy. Constant training columns receive zero Beta weight before
multiplication, including when the beta density is infinite at allele-frequency
endpoints. This applies to both native and Python CPU execution.

This independent example exercises training, imputation, constant columns, and
applying training moments to different values, on both devices::

    import cupy as cp
    import numpy as np
    from pysnptools.snpreader import SnpData
    from pysnptools.standardizer import Unit, Beta
    from pysnptools.util import asnumpy

    iid = [["family", str(sample)] for sample in range(3)]
    sid = ["varying", "constant", "missing"]
    values = np.array([[0, 1, 0], [1, 1, np.nan], [2, 1, 2]], dtype=np.float32)
    for xp in (np, cp):
        for standardizer in (Unit(), Beta(1, 25)):
            data = SnpData(iid=iid, sid=sid, val=xp.array(values, order="F"), xp=xp)
            result, trained = data.standardize(standardizer, return_trained=True)
            assert result is data and isinstance(data.val, xp.ndarray)
            assert data.val.dtype == np.float32 and data.val.flags.f_contiguous
            assert isinstance(trained.stats, np.ndarray)
            assert np.isinf(trained.stats[1, 1])
            np.testing.assert_array_equal(asnumpy(data.val)[:, 1], 0)
            assert float(data.val[1, 2]) == 0
            reordered = SnpData(iid=iid, sid=[sid[2], sid[0]],
                                val=xp.array(values[:, [2, 0]], order="C"), xp=xp)
            if isinstance(standardizer, Unit):
                reordered.standardize(trained)  # UnitTrained remaps known IDs
                np.testing.assert_allclose(asnumpy(reordered.val),
                                           asnumpy(data.val)[:, [2, 0]], rtol=1e-6, atol=1e-6)
            else:
                try:
                    reordered.standardize(trained)
                except AssertionError:
                    pass  # BetaTrained requires exact IDs/order
                else:
                    raise AssertionError("BetaTrained must reject reordered IDs")
            missing_test = SnpData(iid=iid, sid=sid,
                                   val=xp.full((3, 3), xp.nan, dtype=xp.float32), xp=xp)
            missing_test.standardize(trained)  # applying moments does not refit
            np.testing.assert_array_equal(asnumpy(missing_test.val), 0)
            changed = values.copy()
            changed[:, 1] = [0, 1, 2]
            fresh = SnpData(iid=iid, sid=sid, val=xp.array(changed, order="F"), xp=xp)
            fresh.standardize(trained)
            tolerance = 32 * np.finfo(np.float32).eps  # small, well-scaled example
            np.testing.assert_allclose(asnumpy(fresh.val), asnumpy(data.val),
                                       rtol=tolerance, atol=tolerance)
            kernel_source = SnpData(iid=iid, sid=sid,
                                    val=xp.array(changed, order="F"), xp=xp)
            unblocked = kernel_source.read_kernel(trained, dtype="float32", xp=xp)
            blocked = kernel_source.read_kernel(trained, block_size=2,
                                                dtype="float32", xp=xp)
            assert isinstance(blocked.val, xp.ndarray) and blocked.val.dtype == np.float32
            expected_kernel = asnumpy(unblocked.val)
            np.testing.assert_allclose(asnumpy(blocked.val), expected_kernel,
                                       rtol=tolerance,
                                       atol=tolerance * np.max(np.abs(expected_kernel)))
            empty = SnpData(iid=iid, sid=["empty"],
                            val=xp.full((3, 1), xp.nan, dtype=xp.float32), xp=xp)
            try:
                empty.standardize(standardizer)
            except ValueError as error:
                assert "No individual observed" in str(error)
            else:
                raise AssertionError("wholly missing training SNP must fail")

Native GPU BED reads and kernel construction
--------------------------------------------

GPU BED decoding requires the matching bed-reader Python extension built with
``cubecl-cuda``, CuPy, a working CUDA driver, and a supported NVIDIA GPU.
``bed_reader.has_gpu_support()`` reports the compiled capability, not runtime
health. Prebuilt CPU wheels are not sufficient; there is no PySnpTools GPU extra.
Metadata stays on the CPU. Native GPU reads expose Rust-owned device memory
through CuPy, whose views keep that owner alive.

Missing native GPU support raises ``RuntimeError``. Native reads also reject
``force_python_only=True`` and unsupported dtypes/orders with ``ValueError``.
Other file formats can read on CPU and then transfer values; requesting CuPy
does not promise native GPU I/O for every reader. Custom standardizers must
support the chosen array module.

This complete example creates its own BED/BIM/FAM files. It checks selected
values, int8 missing representation, and blocked Unit/Beta kernels against
independent NumPy products. Kernel construction uses an internal fused BED
read/standardize operation; callers use the public read methods::

    import cupy as cp
    import numpy as np
    from pathlib import Path
    from tempfile import TemporaryDirectory
    from bed_reader import has_gpu_support, to_bed
    from pysnptools.snpreader import Bed
    from pysnptools.standardizer import Unit, Beta

    if not has_gpu_support():
        raise RuntimeError("this example requires a bed-reader build with the cubecl-cuda Cargo feature")
    rng = np.random.default_rng(19)
    values = rng.integers(0, 3, size=(12, 9)).astype(np.float64)
    values[0, 0] = np.nan
    with TemporaryDirectory() as directory:
        path = Path(directory) / "cohort.bed"
        to_bed(path, values, count_A1=True)
        bed = Bed(path, count_A1=True)
        data = bed[[4, 1, 4], :].read(dtype="float32", order="F", xp=cp)
        assert isinstance(data.val, cp.ndarray) and data.val.flags.f_contiguous
        assert data.val.dtype == np.float32
        np.testing.assert_array_equal(cp.asnumpy(data.val), values[[4, 1, 4]])
        integer = bed.read(dtype="int8", _require_float32_64=False, xp=cp)
        assert isinstance(integer.val, cp.ndarray) and integer.val.dtype == np.int8
        assert int(integer.val[0, 0]) == -127
        for standardizer in (Unit(), Beta(1, 25)):
            cpu = bed.read(dtype="float32", xp=np).standardize(standardizer)
            kernel = bed.read_kernel(standardizer, block_size=3, dtype="float32",
                                     order="C", xp=cp)
            assert isinstance(kernel.val, cp.ndarray)
            assert kernel.val.dtype == np.float32 and kernel.val.flags.c_contiguous
            # Three blocks of well-scaled float32 values accumulate roundoff.
            tolerance = 64 * np.finfo(np.float32).eps
            np.testing.assert_allclose(cp.asnumpy(kernel.val), cpu.val @ cpu.val.T,
                                       rtol=tolerance, atol=tolerance * bed.sid_count)
            host_kernel = kernel.read(dtype="float32", order="F", xp=np)
            assert isinstance(host_kernel.val, np.ndarray)
            assert host_kernel.val.dtype == np.float32 and host_kernel.val.flags.f_contiguous

``read_kernel(Unit())`` returns the product of standardized SNPs and their
transpose. It does not normalize the average diagonal to one; apply
``kernel.standardize()`` separately for that. The requested positive integer
``block_size`` limits each SNP block, including when SNPs are fewer than samples;
``None`` loads all SNPs. It is not a total memory budget. The dense kernel,
CuPy/CubeCL allocator pools, multiplication temporaries, and downstream solver
workspaces also use GPU memory. A lazy ``SnpKernel`` subset can still construct
the full kernel before selecting its rectangle. Native GPU BED buffers are
limited to ``2**32 - 1`` elements, with further native standardization index
limits; these are indexing limits, not claims that the device has that memory.

PySnpTools forms GPU products and kernel accumulations in the requested floating
dtype, while standardization can use float64 reductions. An output dtype does
not describe every intermediate. A read's ``view_ok`` cannot make a matrix
product into a view. Allocation and CUDA errors propagate rather than selecting
a different backend or reducing the problem silently.
