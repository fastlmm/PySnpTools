# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

* Fix blocked kernel construction with `BetaTrained`: validate the complete
  SNP ID order, reuse the matching fitted moments in each block, and preserve
  the original trained standardizer when returning training information.

* Fix Python CPU Beta standardization of constant SNPs with singular endpoint
  densities (for example, `Beta(0.5, 0.5)`): return zeros without invalid
  multiplications during fitting or application of trained statistics.

* Validate `sub_matrix` indices before conversion on both NumPy and CuPy.
  Reject floating-point indices, Boolean masks, invalid shapes, and out-of-range
  values with `IndexError` instead of truncating, wrapping, or reaching a native panic.
* Preserve float32/float64 when moving numerical values between NumPy and CuPy.
* Support CuPy arrays in `sub_matrix`, preserving the input device, requested
  layout, and dtype; reject negative or out-of-range CuPy indices and warn when
  narrowing float64 to float32.
* Support Beta standardization on CuPy and reject wholly missing SNPs when
  fitting Python/CuPy Unit or Beta standardizers, matching the native path.
* Add a GPU guide with examples and current backend and reader limitations.

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
