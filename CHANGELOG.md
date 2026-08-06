# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
