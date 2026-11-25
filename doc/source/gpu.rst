GPU arrays, readers, and standardization
========================================

PySnpTools reads genotypes on the CPU and transfers the requested data to
CuPy for GPU standardization and kernel construction. NumPy is the default.

**~10× faster with GPU and float32.** On a synthetic dataset with 14,925
individuals and 100,000 SNPs, FaST-LMM's single-kernel ``single_snp`` ran about
10 times as fast with float32 on an RTX 3060 as with the default CPU float64,
and about 5 times as fast as CPU float32. The CPU runs used 16 threads on a
Ryzen 9 7950X.

PySnpTools supplied standardized genotypes in CuPy arrays. Most of the time
savings came from FaST-LMM running its own computations on the GPU, including
the eigendecomposition and association tests.
See the `FaST-LMM documentation <https://fastlmm.github.io/FaST-LMM/>`_.

We recommend float32 for GPU work. Consumer GPUs like the RTX 3060 process
float32 arithmetic much faster than float64, and float32 arrays use half the
memory. Compare against float64 on representative data if precision matters.

Installing GPU support
----------------------

Install PySnpTools with the ``cupy`` extra::

    pip install "pysnptools[cupy]==1.1.0b2"

GPU support requires NVIDIA CUDA on x86-64 Windows or Linux, including WSL2.
Python 3.15 is currently CPU-only. See
`CuPy's installation guide <https://docs.cupy.dev/en/stable/install.html>`_
for driver and CUDA setup. A standard CPU bed-reader installation is sufficient.

Read, standardize, and build a kernel
-------------------------------------

Use your PLINK BED file with its matching BIM and FAM files. ``snps.val`` and
``kernel.val`` are CuPy arrays. Sample and SNP metadata stay on the CPU.

.. code-block:: python

    from pysnptools.snpreader import Bed
    from pysnptools.standardizer import Unit

    bed = Bed("cohort.bed", count_A1=True)
    snps = bed.read(xp="cupy", dtype="float32")
    snps.standardize(Unit())
    kernel = bed.read_kernel(Unit(), xp="cupy", dtype="float32")

``read_kernel(Unit())`` reads and standardizes the SNPs itself. If you only
need the kernel, call it directly on the reader. Set ``block_size`` to limit
how many SNPs it processes at once. The dense sample-by-sample kernel and
CuPy workspaces still need GPU memory.

The :download:`GPU workflow notebook <../ipynb/gpu_workflow.ipynb>` covers
dosages, trained standardizers, blocked kernels, and PCA.

Backend and precision
---------------------

* Explicit ``xp=`` and ``dtype=`` override the process-wide ``ARRAY_MODULE``
  and ``ARRAY_DTYPE`` settings. File reads and kernels default to NumPy and
  float64 when neither is set. Choose ``"cupy"`` or ``"numpy"`` for the backend
  and ``"float32"`` or ``"float64"`` for precision.
* To move existing arrays across devices, pass ``xp=`` explicitly. For example,
  ``kernel.read(xp="numpy", dtype="float32")`` copies the kernel back to the CPU.
  Standardization works in place on the array's current device.
* Use ``pysnptools.util.gpu.has_gpu_support()`` to check for a working CUDA device.
  Requesting CuPy when it is missing raises ``ImportError``. GPU errors propagate
  to the caller.

Reference details
-----------------

Devices and copies
    With ``xp=None``, PySnpTools never copies existing values between the CPU and
    the GPU. A constructor (``SnpData``, ``KernelData``, ``DistData``, ``PstData``)
    given an array from the other module, or a read whose in-memory source is on
    the other module, raises ``ValueError``; pass ``xp=`` to copy. Lists,
    ``val=None``, and values read from files follow ``ARRAY_MODULE``, then NumPy.
    Assigning to ``.val`` requires an array from the object's own module. Writers
    accept CuPy data and copy it to NumPy. Metadata always stays on the CPU.

Precision
    ``ARRAY_DTYPE`` (``float32`` or ``float64``) sets the default for reads,
    kernels, GPU allocation, ``sub_matrix``, and list or integer construction;
    an explicit ``dtype`` wins. Existing floating arrays and trained standardizers
    keep their dtype. GPU float64 reads, allocations, and standardization emit
    ``pysnptools.util.gpu.GpuFloat64Warning`` (a ``UserWarning``) without
    changing the dtype; filter it if float64 is intended. Narrowing float64 to
    float32 warns about precision loss. Float32 storage does not remove every
    float64 operation: reductions can use float64.

Standardizers
    ``Unit`` and ``Beta`` work in place on the array's current device and need
    contiguous float32 or float64 values. ``Unit`` divides by the observed
    population standard deviation (``ddof=0``), then sets NaNs to zero; constant
    columns become zero. ``Beta(a, b)`` centers and multiplies by the beta density
    at the minor allele frequency; on CuPy it reduces on the GPU but evaluates the
    per-SNP density with SciPy on the CPU. Fitting a wholly missing SNP raises
    ``ValueError``; ``Beta`` also rejects means outside [0, 2]. CPU and GPU
    results can differ in the last bits. ``force_python_only`` and
    ``num_threads`` affect only CPU standardization; CuPy arrays stay on the GPU.

Trained statistics
    ``return_trained=True`` returns NumPy statistics (mean and population standard
    deviation per SNP), even after GPU training. ``UnitTrained`` selects and
    reorders known SNP IDs (an unknown ID raises ``KeyError``; reordering with
    repeated training IDs raises ``ValueError``). ``BetaTrained`` requires the
    same SNP IDs in the same order. Applying trained statistics never refits: a
    constant training SNP stays zero, and an all-NaN test column becomes zero.

Kernels and memory
    ``read_kernel`` requires a standardizer and returns standardized SNPs times
    their transpose, in the requested dtype. A positive integer ``block_size``
    limits each SNP block; it is not a memory budget. The dense samples-by-samples
    kernel, CuPy's memory pool, multiplication temporaries, and FaST-LMM's solver
    workspaces also use GPU memory. ``DiagKtoN`` needs a finite, positive trace
    and raises ``ValueError`` otherwise. Allocation and CUDA errors propagate;
    there is no automatic CPU retry.

Readers
    Every file reader decodes on the CPU (bed-reader for BED), then copies the
    requested subset to CuPy, so a custom reader needs only its ordinary CPU
    ``_read`` method. ``force_python_only`` and ``num_threads`` apply to that CPU
    decoding. BED also reads ``int8`` with ``_require_float32_64=False``
    (missing is -127), which must be converted to floating point before
    standardization.
