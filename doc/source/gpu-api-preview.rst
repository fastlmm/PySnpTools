Explicit GPU reader API preview
===============================

This page describes the coordinated GPU API changes currently being developed
alongside this documentation. They are present in the original development
checkout but are not part of this documentation branch's committed code base.
Use :doc:`gpu` for the committed API. Installing the published package is not
sufficient to run this preview; the matching PySnpTools and GPU-enabled
bed-reader implementations are required.

Select the device at each read
------------------------------

``SnpReader.read``, ``SnpReader.read_kernel``, and ``KernelReader.read`` gain an
``xp`` argument accepting NumPy/CuPy modules or their names. An explicit value
overrides ``ARRAY_MODULE``. Otherwise the environment selects the module, with
NumPy as the final default. Specify both ``xp`` and ``dtype`` at read boundaries
so the result does not depend on a surrounding workflow's environment.

This complete in-memory example uses no files. It checks explicit backend
selection against an opposing environment, selects a subset on the GPU, and
standardizes on the array's existing device::

    import cupy as cp
    import numpy as np
    from unittest.mock import patch
    from pysnptools.snpreader import SnpData
    from pysnptools.standardizer import Unit

    assert cp.cuda.runtime.getDeviceCount() > 0
    values = np.array([[0, 2], [1, 1], [2, 0]], dtype=np.float32)
    data = SnpData(iid=[["family", str(sample)] for sample in range(3)],
                   sid=["snp0", "snp1"], val=values, xp=np)
    with patch.dict("os.environ", {"ARRAY_MODULE": "numpy"}):
        gpu = data.read(dtype="float32", xp=cp)
        subset = gpu[[2, 0], :].read(dtype="float32", xp=cp)
        subset.standardize(Unit())
    assert isinstance(subset.val, cp.ndarray) and subset.val.dtype == np.float32
    np.testing.assert_allclose(cp.asnumpy(subset.val), [[1, -1], [-1, 1]],
                               rtol=1e-6, atol=1e-6)
    with patch.dict("os.environ", {"ARRAY_MODULE": "cupy"}):
        host = gpu.read(dtype="float32", xp=np)
        host.standardize(Unit())
    assert isinstance(host.val, np.ndarray)

The standardization change is significant: an existing array determines the
device, regardless of ``ARRAY_MODULE``. It resolves the Unit device mismatch
described in the committed guide. Subsetting still needs a read to materialize
values, and that read has its own device and dtype choice. A default float64
read does not preserve float32 merely because the input is float32.

Read BED data and form a kernel
-------------------------------

With the matching CUDA-enabled bed-reader build and CuPy, BED reads request
GPU decoding directly. Kernel construction can read and standardize BED blocks
on the GPU before CuPy multiplication. A separate public ``read().standardize()``
sequence remains supported; callers do not need private fused-read methods.
For an existing BED dataset with accompanying BIM/FAM files::

    from pysnptools.snpreader import Bed
    from pysnptools.kernelreader import SnpKernel

    bed = Bed("cohort.bed", count_A1=True)
    block = bed[:100, :200].read(dtype="float32", xp=cp)
    block.standardize(Unit())
    kernel = bed.read_kernel(Unit(), block_size=7, dtype="float32", xp=cp)
    kernel_again = SnpKernel(bed, Unit(), block_size=7).read(
        dtype="float32", xp=cp
    )
    assert isinstance(kernel.val, cp.ndarray) and kernel.val.dtype == np.float32
    np.testing.assert_allclose(cp.asnumpy(kernel.val), cp.asnumpy(kernel_again.val),
                               rtol=1e-5, atol=1e-5)

The requested SNP block size is respected even when there are fewer SNPs than
samples, unlike the shortcut described in the committed guide. It limits the
SNP blocks, not the dense kernel, eigensolver workspace, or allocator pools.
Only BED decoding is native GPU I/O here; other formats can read on the CPU and
transfer values to the selected device.

CuPy Unit/Beta standardization also supports general floating-point arrays,
including dosages, with per-SNP workspaces. Missing values become zero after
standardization; fitting a wholly missing SNP raises. Applying trained
statistics preserves the training SNP order and constant-column behavior.
The new CuPy Beta path evaluates its weights on the GPU, replacing the committed
implementation's SciPy round trip. Float32 reductions and general floating-point
statistics need not be bit-identical to CPU results.

A BED GPU request without a GPU-enabled native build raises instead of quietly
reading the BED values on the CPU. This differs from the committed FaST-LMM
native-reader fallback. Module-name resolution can still fall back to NumPy
when importing CuPy fails; import CuPy explicitly and pass the module when
checking GPU execution. Driver errors and insufficient memory remain errors.
