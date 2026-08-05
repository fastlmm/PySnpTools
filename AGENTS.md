# Coding Notes for Agents

This file contains repository-wide guidance for PySnpTools, a mature numerical
Python package used directly by FaST-LMM. Preserve compatibility, data
integrity, and numerical correctness unless a task explicitly changes those
goals.

## General Policies

- When work is interrupted or reaches a stopping point, report the current
  status, what remains, and the recommended next step. If the next step is
  within the current task and safe to perform, perform it instead of merely
  recommending it.
- Inspect existing code, tests, documentation, project configuration, and
  package-data rules before introducing a new pattern. Prefer focused changes
  that follow nearby conventions over unrelated cleanup or broad rewrites.
- Do not silently skip required Python versions, operating systems,
  architectures, dependency configurations, tests, doctests, or artifact
  checks. If a required environment or tool is unavailable, fail clearly and
  report what is missing.
- Avoid silent clamping, coercion, truncation, precision loss, or fallback
  behavior. Validate inputs and fail clearly unless the public API explicitly
  documents another behavior.
- Preserve behavior on every supported Python version, not only the interpreter
  used for development.

## Python 3.14 and FaST-LMM Coordination

For Python 3.14 support work, read the coordination and prerequisite-release
requirements in
[`../FaST-LMM/specs/PYTHON_3_14_SUPPORT_SPEC.md`](../FaST-LMM/specs/PYTHON_3_14_SUPPORT_SPEC.md)
before making changes.

PySnpTools must publish a Python 3.14-compatible release before FaST-LMM can
complete its Python 3.14 release. Support must cover Python 3.10 through 3.14
and every maintained platform and architecture. Keep the FaST-LMM specification
aligned with implementation decisions that materially change dependency
requirements, the validation matrix, artifact coverage, or the release plan.

Final acceptance must install and test published `bed-reader` and PySnpTools
artifacts rather than relying on unpublished sibling checkouts. Do not make
unnecessary bed-reader changes as part of the FaST-LMM support work; use its
published Python 3.14-compatible release unless a demonstrated blocker requires
otherwise.

## Genetic-Data and Numerical Correctness

- Preserve sample and marker ordering, IID/SID alignment, allele counting,
  missing-value behavior, position metadata, array shapes, dtypes, and memory
  ordering across readers, writers, standardizers, and kernels.
- Treat changes to BED, BGEN, HDF5, NPZ, distributed data, caching, map-reduce,
  standardization, and indexing behavior as behavior changes rather than
  mechanical refactors.
- Do not update expected numerical or textual output merely to make a test pass.
  First confirm that the new result is correct and does not hide a precision,
  ordering, alignment, or compatibility regression.
- Use deterministic seeds in tests involving randomness. When exact equality is
  inappropriate, use an explicit, justified tolerance rather than an
  unnecessarily broad one.
- Preserve representative small read/write round trips and end-to-end tests for
  each supported data format and optional feature.

## Error Handling

- Preserve useful diagnostics and original exceptions when translating errors.
  Use exception chaining (`raise ... from error`) when adding context.
- Catch only exceptions that can be handled meaningfully. Do not replace useful
  failures with generic messages, sentinel values, empty data, or broad
  catch-all behavior.
- Do not suppress warnings solely to quiet tests or CI. Fix the cause or
  document a narrow, justified suppression.
- File, URL, cache, and distributed-execution errors should identify the
  relevant operation and resource without leaking credentials or other secrets.

## Dependencies and Packaging

- Treat NumPy, SciPy, pandas, h5py, bed-reader, cloudpickle, and optional BGEN
  dependency upgrades as behavior migrations, not just resolution or import
  fixes. Review changed APIs and defaults, then test affected behavior.
- Keep dependency-upgrade changes focused when practical. Avoid mixing them
  with broad code or documentation churn unless required by the migration.
- Declare direct runtime imports as direct project dependencies with honest
  lower bounds. Use Python-version markers when a newer interpreter needs a
  newer dependency without unnecessarily raising requirements elsewhere.
- Keep feature dependencies such as BGEN support in published optional extras;
  keep development-only tools in development dependency groups rather than
  runtime metadata.
- Keep `pyproject.toml`, the lockfile, classifiers, package discovery,
  package-data configuration, CI, and documentation consistent.
- Build and test both the wheel and source distribution in isolated
  environments. Verify metadata, licenses, authorship files, hashdown manifests,
  test fixtures required at runtime, and other intended package data.
- Test installed artifacts without relying on the repository being on
  `PYTHONPATH`. A passing test against the source checkout is not sufficient
  evidence that the package is correct.

## Test Data, Hashdown Manifests, and Generated Files

- Treat hashdown JSON files as integrity-sensitive manifests. When their source
  URLs, commit pins, file lists, or hashes change, verify every affected
  download and explain why the change is required.
- Keep remote test-data references immutable where practical. Do not replace a
  commit-pinned URL with a moving branch merely for convenience.
- Do not casually edit binary genetic datasets, generated documentation,
  downloaded fixtures, or cache contents. Identify the source of truth first.
- Treat built documentation under generated output directories as output, not
  source. Edit the corresponding Sphinx source or configuration, regenerate,
  and inspect the resulting diff.
- If generated output must be patched urgently, make the matching source or
  generator change in the same work so regeneration does not revert it.

## API and Code Design

- Preserve documented public APIs and serialized or file-format compatibility.
  Treat changes to public signatures, defaults, accepted input shapes or types,
  result schemas, import paths, and exceptions as compatibility-sensitive.
- Prefer one clear canonical API path. Keep aliases or compatibility shims only
  when downstream users require them, and document the canonical path.
- Keep implementation details private. Do not expose helpers publicly merely to
  make internal code or tests convenient.
- Place functionality with the abstraction it belongs to. When moving an
  abstraction, account for related helpers, tests, examples, and documentation
  rather than leaving partial duplicate implementations.
- Prefer concise names when module or class context already supplies context.
  Avoid obscure abbreviations while retaining conventional genetics and linear
  algebra notation where it is established and unambiguous nearby.

## Tests and Validation

- Add or update tests for every behavior change and regression fix. A regression
  test should fail for the original defect and pass for the corrected behavior.
- Run the narrowest relevant tests while iterating, then run the repository's
  complete required validation appropriate to the change before handing work
  back.
- Compare the suites assembled by legacy `tests/test.py` with pytest collection
  before changing test orchestration. Do not retire either path until the
  replacement demonstrably covers the required unit tests, integration tests,
  and doctests.
- Test both minimal dependencies and all optional dependencies. BGEN tests must
  be clearly distinguished so missing native support cannot silently reduce
  required coverage.
- For packaging and dependency work, validate supported boundary interpreters,
  including Python 3.10 and 3.14 during the current modernization, and test the
  actual wheel and source distribution.
- Do not disable, deselect, or weaken tests, doctests, warnings, lint rules, or
  CI jobs merely to obtain a passing result. Any exception must be narrow,
  documented, and reported.

## Comments and Documentation

- Preserve useful TODOs, diagnostic notes, data-format explanations, and
  debugging comparisons while their underlying issue remains unresolved.
- Keep public documentation synchronized with signatures, defaults, exceptions,
  supported versions, optional dependencies, installation steps, and examples.
- Prefer one complete, executable example for a public workflow and link to it
  from related documentation instead of maintaining multiple drifting copies.
- Use American English. In Markdown, place blank lines around headings, lists,
  and fenced code blocks, and keep list markers consistent.

## Development and Release Safety

- Use project-local environments and the repository's documented toolchain. Do
  not silently install system-wide packages or alter unrelated global
  configuration.
- Do not publish a real release to PyPI or another package index. Agents may
  prepare version changes, release notes, artifacts, and commands, but a person
  must perform or explicitly authorize publication.
- Ordinary CI and release qualification must not rely on unpublished sibling
  checkouts, unreviewed prereleases, or accidental local state.
- Before handing work back, summarize validation performed, anything not run,
  and any data-integrity, numerical, compatibility, packaging, or release risk
  that remains.
