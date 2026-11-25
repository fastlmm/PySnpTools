# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.1.0] - in beta (changes since 1.0.1)

### Highlights

* **GPU support with CuPy.** Every reader accepts `xp="cupy"`: it decodes on the CPU
  (bed-reader is CPU-only), copies the data to the GPU, and standardization and kernel
  construction run there. Install the `cupy` extra (`pip install "pysnptools[cupy]"`).
  FaST-LMM's single-kernel `single_snp` uses this to run about 10 times as fast on a GPU
  with float32 as on the CPU with float64. See the GPU guide.
* **Backend and precision defaults.** The `ARRAY_MODULE` and `ARRAY_DTYPE` environment
  variables set the defaults for reads, kernels, `sub_matrix`, and new containers; an
  explicit `xp` or `dtype` wins. `DistData`, `PstData`, `DistReader.read`, and
  `PstReader.read` gain an `xp` parameter. Float32 and float64 are preserved when
  moving values between NumPy and CuPy.

### Behavior changes

* Requesting CuPy without it installed (`xp="cupy"`, `array_module("cupy")`, or
  `ARRAY_MODULE=cupy`) raises `ImportError`. 1.0.1 logged a warning and used NumPy.
* With `xp=None`, existing values are never copied between the CPU and the GPU.
  Constructing a container from an array of the other module, reading such in-memory
  values, or assigning one to `.val` raises `ValueError`; pass `xp=` to copy. Lists,
  `val=None`, and values read from files follow `ARRAY_MODULE`.
* `Unit().standardize` no longer moves NumPy values to CuPy when `ARRAY_MODULE=cupy`;
  every standardizer works on the values' existing device. Trained statistics are
  always NumPy arrays.
* Rename the local runners' thread limit from `mkl_num_threads` to `num_threads`
  (the old name remains as a deprecated alias). The limit now also covers OpenBLAS,
  BLIS, OpenMP, spawned processes, and reducers. `LocalMultiThread` runs each task in a
  copy of the caller's context.
* `read_kernel` honors an explicit `block_size` even when there are fewer SNPs than
  individuals, so such kernels can differ in the last bits of rounding. Invalid block
  sizes, including negative ones that could silently give an all-zero kernel, raise
  `ValueError`. `SnpData.standardize(block_size=...)` warns that `block_size` is
  deprecated, like the standardizers.
* `sub_matrix` rejects floating-point indices, Boolean masks, invalid shapes, and
  out-of-range values with `IndexError` instead of truncating, wrapping, or panicking.
  It also accepts CuPy arrays, keeping them on the GPU.
* GPU float64 work emits `pysnptools.util.gpu.GpuFloat64Warning`, because float64 can
  be much slower than float32 on consumer GPUs. `pysnptools.util.gpu.has_gpu_support()`
  reports whether CuPy can use the current device.

### Fixes

* Constant fractional SNPs (for example, dosages that are all 0.1) now standardize to
  zero with Unit and Beta on every path and train an infinite deviation. Before, a
  rounded mean left a tiny variance, so they standardized to +/-1 on the Python path
  and kept small residues on the native path. Constant 0/1/2 genotypes were unaffected.
* Blocked kernels with `DiagKtoN` as the SNP standardizer fitted a separate factor to
  each block; they are now scaled once and match the unblocked kernel. Blocked kernels
  with `BetaTrained` now validate the SNP order and reuse the matching fitted moments.
* `DiagKtoN` rejects zero, negative, and nonfinite normalization denominators with
  `ValueError` before changing values.
* Applying Unit or Beta trained statistics works for Fortran-order targets.
* Python Beta standardization of constant SNPs with singular endpoint densities (for
  example, `Beta(0.5, 0.5)`) returns zeros without invalid arithmetic, and rejects
  allele-count means outside [0, 2], matching the native path.
* `UnitTrained` rejects reordering by a repeated training SNP ID, and rejects zero or
  otherwise invalid manually supplied deviations before changing values. Fitting a
  wholly missing SNP raises on every path.
* Every writer accepts CuPy-backed data and copies it to NumPy. Containers keep their
  backend across pickling.
* Pheno, Dense, Dat, and Ped text parsing stays on the CPU until a read chooses its
  backend.

### Other

* Require bed-reader 1.3 (during the beta, 1.3.0b3 or later) and threadpoolctl.
* `sub_matrix` validates plain-integer list, tuple, and range indices without a
  per-element Python loop.
* Wheels and source distributions exclude generated `.pytest_cache` and `__pycache__`
  folders.
* New GPU guide and GPU workflow notebook.

## [1.0.1] - 2026-11-14

* Make `Hashdown` downloads (used by `example_file`) more robust: retry
  network errors and HTTP 5xx responses up to three times with backoff, time
  out stalled connections, and write each file under a temporary name so an
  interrupted download never leaves a partial file in the cache.

* Correct 69 stale sample-file hashes to match the existing immutable source
  commit, restoring example downloads without changing the sample data. Add
  regression checks against generated data and distributed metadata.

* Require h5py 3.15.0 or later for the HDF5 extra on Python 3.14 and
  later, while retaining the existing requirement on earlier Python versions.

## [1.0.0] - 2026-11-13

* Require Python 3.11.5 or later and refresh stable runtime and development
  dependencies, including security fixes for Pillow and Soup Sieve.
* Make HDF5 support optional: existing HDF5 users must install
  `pysnptools[hdf5]`. BGEN remains available through `pysnptools[bgen]`.
* Load optional backends only when used, with actionable installation errors
  and preservation of failures inside an installed backend.
* Prepare Python 3.15 core support; optional formats remain qualified on
  Python 3.11 through 3.14 pending upstream binary support.
* Require bed-reader 1.2.0 or later, the first release with Python 3.15 and
  Windows ARM64 wheels.
* Remove the use of `six` from `LocalMultiThread`, so PySnpTools no longer
  needs it (0.5.15 imported it without declaring it).

## [0.5.15] - 2026-08-06

* Add support for Python 3.14 while retaining Python 3.10 through 3.13.
* Use a Python-version-specific NumPy requirement for Python 3.14 compatibility.
* Repair the tutorial sample-data downloads with immutable references and add
  `example_file_synth`.
* Migrate the build backend to `uv_build` and test wheel and source artifacts
  in clean Python 3.10 and 3.14 environments.
* Modernize dependency groups, package metadata, and the cross-platform CI
  matrix.
* Remove `wheel` from runtime dependencies.

## [0.5.14] - 2024-11-2

Add support for Python 3.13.

## [0.5.12] - 2024-5-28

### Changed 0.5.12

* Make support for BGEN format optional.

## [0.5.11] - 2023-11-8

### Changed 0.5.11

* Improved code formatting and linting.
* Tested with bed-reader 1.0.0.
* Removed dependency 'six'.
* Fixed minor bug <https://github.com/fastlmm/PySnpTools/issues/5>
* Updated sample `*.npz` and `*.memmap` files to Python 3 format.
