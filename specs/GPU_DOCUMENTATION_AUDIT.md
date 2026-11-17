<!-- TODO0 consider deleting this spec once the work below is implemented and released. -->

# GPU documentation review

Reviewed November 16, 2026, from `8c68afc` (immediately before GPU work)
through `366d6bb`, plus the current working changes. Nothing was published.

## Scope and corrections

Read the README, changelog, documentation landing page, GPU guide, former API
preview, and changed docstrings in SNP/kernel readers, data containers,
standardizers, and array utilities. Checked their claims against implementation
and examples. The current inventory includes 34 changed class/method docstrings.

- Kept the README's GPU introduction to one feature bullet linking the guide
  and notebook.
- Corrected reader summaries that promised NumPy-only results or a fresh values
  allocation even when `view_ok=True` permits sharing.
- Explained independent per-read backend/dtype selection, float32 preservation,
  CPU metadata, Cargo feature requirements, native CUDA device 0, ignored CPU
  options, missing/constant SNPs, and trained-statistics precision and alignment.
- Documented `SnpKernel.read_snps(xp=...)`, its ignored `view_ok` argument, and
  GPU behavior, with a self-contained doctest.
- Corrected normalization prose and a missing NumPy import. Distinguished fitting
  `DiagKtoN` from reusing its scalar on different data.
- Added direct notebook downloads to the generated landing page and guide.

## Confirmed bugs fixed

1. **Invalid kernel block sizes:** a negative size could silently return an
   all-zero kernel. Validate at `SnpKernel` construction, before reads and
   allocations. Reject nonpositive, noninteger, and Boolean sizes with
   `ValueError`; keep `None` and positive Python/NumPy integers. The regression
   checks the public eager/lazy APIs and an analytic nonzero kernel.
2. **Trained statistics across dtypes:** a misplaced Boolean grouping skipped
   dtype/ownership checks for F-order statistics. Applying float32 statistics
   to float64 native CPU targets (and the reverse) failed with a misleading
   ndarray type error. Both layouts now require a compatible working statistics
   array. Original fitted moments remain unchanged. The new regression failed
   in four native F-order cases before the fix; it checks both standardizers,
   both conversion directions, C/F layout, native/Python execution, missing
   values, and a training constant against explicit formulas.

Neither fix changes the numerical algorithm for valid, already-compatible input.
Widening the target dtype does not recover precision lost during training.

## Normalization follow-up fixed

<!-- TODO0 define and enforce consistent invalid normalization behavior before release (implemented; may no longer apply). -->

Before the follow-up, `DiagKtoN` had inconsistent behavior for a zero-trace kernel. With a 2-by-2 zero
`KernelData`, NumPy execution produces NaNs with divide-by-zero/invalid warnings;
CuPy execution raises `ZeroDivisionError`. This was reproduced on the installed
development wheel. The guide/example uses finite positive-trace data, and the
type documentation states that precondition.

The implementation now validates the existing scalar reduction before division
and mutation. Both devices raise `ValueError` for zero, negative, NaN, or infinite
denominators and for empty inputs, preserving values on rejection. A separate GPU
sum-of-squares reduction propagates NaNs instead of skipping them as the reused
missing-aware variance reduction did. It keeps float64 GPU accumulation and
uses no additional matrix scan or device-to-host transfer. Valid CPU inputs keep
their established scalar dtype and rounding. A nonfinite or underflowed scaling
factor is also rejected. Applying an already-fitted scalar to zero test data
remains supported.

Focused CPU/GPU regressions passed: both dtypes, C/F order, SNP/kernel data,
both return modes, invalid/empty inputs, analytic valid results, and trained
reuse on zero test data (176 subtests), including against the installed wheel.
The full follow-up suite passed 754 tests in 297 seconds; the 42 focused GPU
cases, complete installed-wheel notebook, and strict Sphinx build also passed.
These results and the changed behavior are reported directly in the conversation.

Legacy API concerns remain visible in the guide: `SnpData.standardize(block_size)`
and `SnpKernel.read_snps(view_ok)` ignore those arguments, while direct Unit/Beta
calls warn about `block_size`. A future deprecation pass should make those warnings
consistent without silently changing sharing or standardization behavior.

## Notebook

`doc/ipynb/gpu_workflow.ipynb` contains 20 cells and a rendered synthetic PCA plot:

1. Check the native build and actual CUDA execution; report package/device versions.
2. Write local BED/BIM/FAM data under the user's cache, with missing and constant SNPs.
3. Decode on the GPU; check exact hard calls, duplicate subsets, dtypes, C/F layout,
   environment overrides, CPU metadata, and device-preserving `sub_matrix`.
4. Check fractional Unit/Beta inputs against independent centered-variance and
   polynomial beta-density formulas in both floating dtypes and memory orders.
5. Fit on training rows, apply without refitting, remap UnitTrained SNP IDs, and
   reuse GPU-trained float32 moments on a CPU float64 target.
6. Build held-out BetaTrained kernels in blocks and compare against an independent
   float64 product. Build and normalize a full-cohort Unit kernel for GPU PCA.
7. Plot only two transferred coordinates; verify eigenpair residuals on the GPU.
8. Demonstrate warm wall-clock timing with device synchronization before and after
   each operation, without claiming a tiny example is faster than the CPU.

The complete notebook passed in a fresh kernel outside the source checkout using
the installed wheel, CPython 3.13, CuPy 14.2, and an NVIDIA RTX 3060. Hard calls use
exact comparisons. Floating checks state roundoff-based tolerances and kernel
absolute-product scales. Temporary BED files are removed at the end.

## Validation

- Wheel and source distribution built successfully; the isolated environment was
  refreshed with the development wheel. The notebook and all three independent
  guide examples passed against that installed wheel.
- The 42 focused GPU API/general-standardization regression cases passed.
- Both newly added CPU regressions passed, including 24 subtests. Final complete
  source canonical validation passed 752 tests in 285 seconds, with HDF5 and BGEN
  dependencies available and no skipped tests.
- Strict Sphinx HTML generation passed. Local links/fragments and notebook
  downloads resolve. Rendered reader and `read_snps` signatures expose `xp=None`
  as the last parameter. The PCA output was visually inspected.
- The canonical suite needs legacy source fixtures absent from the wheel; an
  attempted installed-wheel canonical invocation failed during fixture setup.
  Full canonical coverage therefore uses the source checkout. This does not
  replace the separate installed-wheel notebook and GPU tests.

Logs and generated HTML are under `~/.cache/python-docs-audit/`, including
`pysnptools-doc-review-{final-canonical,gpu,sphinx,examples,notebook}.log` and
`PySnpTools/review-html/`. The reusable notebook runner is
`/mnt/c/deldir/python-docs-audit/run_pysnptools_notebook.py`.

Not qualified here: the complete Python/OS matrix, HIP, other GPUs, large or
near-capacity performance, or acceptance against future published GPU packages.
