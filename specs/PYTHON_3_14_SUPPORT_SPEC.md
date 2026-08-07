# PySnpTools Python 3.14 Support Specification

<!-- Short-term implementation plan. Remove or archive after the release. -->

## Status

Proposed.

## Objective

Release PySnpTools with supported, non-experimental Python 3.14 compatibility
while retaining Python 3.10, 3.11, 3.12, and 3.13 support. The release must
install from its published wheel and source distribution, preserve genetic-data
and numerical behavior, and be suitable for use as a released FaST-LMM
dependency.

This is the repository-local implementation plan. Cross-repository sequencing
remains: publish compatible PySnpTools and fastlmmclib releases before final
FaST-LMM qualification, then test FaST-LMM against those published packages.

## Current State

- `pyproject.toml` requires Python 3.10 or newer and advertises Python 3.10
  through 3.13.
- CI tests Python 3.10 through 3.13 on Linux, Windows, Intel macOS, and Apple
  Silicon macOS, but does not test Python 3.14.
- CI uses `astral-sh/setup-uv@v3`, separately installs Python, manually
  activates `.venv`, and runs lint and builds redundantly in every matrix job.
- CI runs the legacy `tests/test.py` suite with minimal dependencies, then
  only `pysnptools/distreader/test_bgen2.py` after installing all extras. It
  has not established that normal pytest collection, doctests, and the legacy
  suite provide equivalent coverage.
- `uv.lock` is ignored, so ordinary CI does not use a committed frozen
  resolution.
- The package uses setuptools. NumPy is unnecessarily listed as an isolated
  build requirement even though no PySnpTools extension build has been
  identified.
- `wheel` is declared as a runtime dependency even though it is a packaging
  tool and no runtime import has been identified.
- Development tools are split between a published `dev` extra and legacy
  `[tool.uv].dev-dependencies`.
- License metadata uses the deprecated table form.
- The workflow uploads only an sdist and has no isolated artifact tests or
  Trusted Publishing release job.
- The two hashdown manifests are runtime package data and point at immutable
  remote objects by content hash.
- The tutorial request described in
  [issue #10](https://github.com/fastlmm/PySnpTools/issues/10) is broken:
  `pysnptools.hashdown.json` lists `doc/ipynb/all.bed`, `.bim`, and
  `.fam` at a pinned source location where they are unavailable.

## Supported Matrix

Required support is:

| Platform | Python versions |
| --- | --- |
| Linux x86-64 | 3.10, 3.11, 3.12, 3.13, 3.14 |
| Windows x86-64 | 3.10, 3.11, 3.12, 3.13, 3.14 |
| macOS Intel | 3.10, 3.11, 3.12, 3.13, 3.14 |
| macOS Apple Silicon | 3.10, 3.11, 3.12, 3.13, 3.14 |

Use explicit, currently supported GitHub-hosted runner labels. Recheck runner
availability immediately before implementation. Python 3.14 jobs are required
and must not use `continue-on-error`.

## Tutorial Data and Hashdown Integrity

Resolve issue #10 using this implementation:

1. Use the five recovered root-level files in
   [`fastlmm/bed-sample-files`](https://github.com/fastlmm/bed-sample-files)
   as the source of truth.
2. Pin `synth.hashdown.json` to immutable commit
   `b17af28faf1ea6152fc84010dbf2f3fa38a84fe1`; do not use the moving
   `main` branch.
3. Keep one helper module and expose `example_file_synth(...)` beside the
   existing `example_file(...)` and `example_file_bgen(...)` functions.
4. Update the tutorial and paper notebooks to request the root-level filenames
   through `example_file_synth(...)`.
5. Remove the five unavailable `doc/ipynb/...` entries from
   `pysnptools.hashdown.json`.
6. Verify the exact served bytes against the historical hashes. The recovered
   Windows text files must use their original LF representation.
7. Add focused retrieval tests for BED, BIM, FAM, phenotype, and covariate
   files. Keep network retrieval in an appropriate integration or release-gate
   job if ordinary offline unit-test reliability requires separation.
8. Validate every entry in all three hashdown manifests before release, or
   document a deliberately narrower validation set and why excluded data
   cannot be checked.
9. Execute both `doc/ipynb/tutorial.ipynb` and `doc/ipynb/paper.ipynb`
   from start to finish in a clean environment with an empty
   `PYSNPTOOLS_CACHE_HOME`. Confirm that every cell succeeds, the synthetic
   files are downloaded from the pinned sample-data revision, and the computed
   results remain correct. Regenerate and inspect stored notebook outputs only
   after the execution succeeds.

Close issue #10 only after the fixed immutable references and regression test
are present.

## Dependencies

### NumPy

Declare the coordinated Python-version-specific runtime requirement:

```toml
"numpy>=1.22.0; python_version < '3.14'",
"numpy>=2.3.5; python_version >= '3.14'",
```

Do not raise the minimum for Python 3.10 through 3.13 solely for Python 3.14.
Test NumPy 2.x behavior involving dtypes, memory order, missing values,
standardization, indexing, and serialized NPZ data.

### Remaining dependencies

Verify that stable releases resolve, install, and pass affected tests across
the supported matrix for:

- SciPy
- pandas
- h5py
- psutil
- cloudpickle
- more-itertools
- `bed-reader[samples]`
- optional `cbgen` and `bgen-reader`
- all test, lint, documentation, and build tools

Determine and test honest direct-dependency lower bounds. Use Python-version
markers where Python 3.14 needs a newer release without unnecessarily raising
older interpreters' bounds. Required CI and release qualification must use
stable published packages, not VCS dependencies or prereleases.

### Direct dependency refresh checklist

Review every current runtime requirement during this release and record the
tested decision in the implementation change:

| Dependency | Required decision |
| --- | --- |
| NumPy | Keep the explicit pre-3.14 and 3.14-or-newer markers above. |
| SciPy | Test the current stable line, including the 1.18 line identified during planning, and raise or mark the lower bound only where compatibility evidence requires it. |
| pandas | Test the current stable line, including the 3.0 line identified during planning, and preserve the older lower bound where it remains truthful. |
| psutil | Test both the declared minimum and current stable release; retain or raise the minimum based on results. |
| h5py | Verify stable Python 3.14 wheels and test HDF5 reads and writes at the declared boundary and current stable release. |
| cloudpickle | Test the declared minimum and current stable release, including multiprocessing and map-reduce serialization paths. |
| `bed-reader[samples]` | Require a published Python 3.14-compatible release and test PySnpTools against that artifact. |
| more-itertools | Test the declared minimum and current stable release; do not raise the minimum without a code or support requirement. |
| wheel | Remove it from runtime dependencies unless an actual runtime import is demonstrated. |

Do not mechanically replace every lower bound with the newest version shown by
an editor. The normal stable-resolution job proves compatibility with current
releases; the lower-bound job proves that published minimums remain honest.
Use a Python-version marker when only Python 3.14 requires the newer line.

Remove `wheel` from runtime dependencies unless a runtime use is demonstrated.
Remove NumPy from `[build-system].requires` when changing the backend; restore
it only if an isolated build proves it is required and document the reason.

### Optional BGEN support

Determine whether stable `cbgen` and `bgen-reader` artifacts support Python
3.14 on every claimed platform. If they do, run the BGEN suite across the
supported combinations. If they do not, make the limitation explicit with
accurate markers, CI, documentation, and release notes. Do not silently skip
BGEN tests or use unpublished artifacts to claim support.

## Project and Build Tooling

Use `uv` as the standard interface for interpreter installation, dependency
resolution, environments, commands, builds, and publishing.

- Pin a reviewed `uv` version.
- Pin GitHub Actions to reviewed full commit SHAs with release comments.
- Select the matrix interpreter through the current `setup-uv` action.
- Use `uv run`; do not manually activate virtual environments.
- Remove `uv.lock` from `.gitignore`, generate it for the complete supported
  matrix, and commit it.
- Use `uv sync --frozen` and `uv run --frozen` in ordinary CI.
- Add required lower-bound testing with `--resolution lowest-direct` on the
  applicable Python boundary versions.
- Add a scheduled highest-stable dependency solve to detect ecosystem drift
  without rewriting the committed lockfile in ordinary CI.
- Keep prerelease ecosystem testing separate and non-blocking if it is retained.

Migrate the pure-Python package from setuptools to `uv_build` unless package
inspection finds an unsupported requirement. Configure the existing flat
`pysnptools` package layout explicitly. Remove the parallel setuptools path
after isolated builds prove the migration.

During the backend migration:

1. Preserve the intended package modules.
2. Include `LICENSE.md`, `AUTHORS.txt`, applicable RST content,
   `pysnptools.hashdown.json`, `bgen.hashdown.json`, and any runtime test
   helper intentionally exposed by the installed package.
3. Exclude repository tests, generated documentation under `doc/build`,
   notebooks, caches, and development-only files unless they are intentionally
   required in the sdist.
4. Compare wheel and sdist manifests with version 0.5.14 artifacts and explain
   every intentional addition or removal.
5. Build with no unpublished source overrides.
6. Build a wheel from the generated sdist and verify that it is equivalent to
   the wheel built directly from the checkout.

Replace legacy `[tool.uv].dev-dependencies` with PEP 735 dependency groups.
Separate test, lint, and documentation tools. Retain a published `dev` extra
only if downstream use justifies it; feature extras such as `bgen` remain
published optional dependencies. Verify that normal `uv sync` and `uv run`
commands no longer emit the `tool.uv.dev-dependencies` deprecation warning.

## Metadata

Update package metadata to:

- Add the Python 3.14 classifier and retain Python 3.10 through 3.13.
- Keep `requires-python = ">=3.10"`.
- Apply tested dependency bounds and markers.
- Replace the deprecated license table with PEP 639 metadata:

  ```toml
  license = "Apache-2.0"
  license-files = ["LICENSE.md"]
  ```

- Use HTTPS project URLs and conventional labels such as `Homepage`,
  `Documentation`, `Issues`, and `Source`.

Inspect the built metadata rather than validating only `pyproject.toml`.

## Test-Suite Work

Keep the established unittest orchestration for this release:

```bash
uv run python tests/test.py
```

The legacy runner passed 725 tests on Python 3.14.5. Root pytest collection
found only 199 items, and explicitly collecting `pysnptools/test.py` failed
because pytest attempted to construct the parameterized
`NaNCNCTestCases` class without its required arguments. The legacy suite
constructs those cases correctly. This is a test-harness incompatibility, not a
product failure.

Do not migrate the legacy suite to pytest as part of Python 3.14 support. Keep
`tests/test.py` as the canonical core command and run the optional BGEN tests
separately with:

```bash
uv run pytest pysnptools/distreader/test_bgen2.py
```

Confirm that the legacy suite continues to execute its intended doctest cases,
and keep the BGEN job visibly separate so optional-dependency coverage cannot
silently disappear. A future test-infrastructure project may redesign the
parameterized unittest cases and reconsider unified pytest discovery.

Required coverage must include representative:

- BED read/write round trips and sample/variant subsetting.
- Standardization for float32 and float64.
- IID/SID alignment, missing values, and array order.
- HDF5 and NPZ readers and writers.
- Kernel and distance readers.
- Map-reduce and cache behavior that is portable in hosted CI.
- Optional BGEN behavior when supported.
- Hashdown retrieval and integrity validation.

Do not update numerical expectations merely to make Python 3.14 pass. Confirm
the cause and correctness of every changed tolerance or expected value.

## Continuous Integration

Split the current all-in-one matrix into focused jobs:

1. Run one pinned Ruff check.
2. Run the complete test matrix on Python 3.10 through 3.14 across Linux,
   Windows, Intel macOS, and Apple Silicon macOS using
   `uv run python tests/test.py`.
3. Run lower-bound tests on applicable boundary interpreters.
4. Run optional BGEN tests distinctly so missing optional dependencies cannot
   silently reduce core coverage, using
   `uv run pytest pysnptools/distreader/test_bgen2.py`.
5. Build wheel and sdist once.
6. Install and smoke-test the exact artifacts in clean Python 3.10 and 3.14
   environments without the repository on `PYTHONPATH`.
7. Run the tutorial-data retrieval integration test.
8. Run a scheduled highest-stable dependency test.

Set minimal permissions, timeouts, concurrency cancellation, frozen resolution,
and action SHA pins. Upload both wheel and sdist with unambiguous artifact
names.

The installed-artifact smoke test must import major public modules, load a tiny
dataset, perform an inexpensive read or standardization, and verify the
installed hashdown manifests are present. Test both the directly built wheel
and a wheel rebuilt from the sdist.

## Release

Add a Trusted Publishing workflow:

1. Trigger from a matching immutable version tag.
2. Build wheel and sdist once.
3. Test those exact files.
4. Transfer the tested files to a separate protected publishing job without
   rebuilding.
5. Grant `id-token: write` only to that job.
6. Publish with PyPI Trusted Publishing and retain attestations.
7. Never publish from pull-request workflows.

The human maintainer chooses the version and authorizes publication. The tag
must point to the exact source commit used for the PyPI artifacts. Add release
notes covering Python 3.14, dependency-bound changes, the tutorial-data repair,
the BGEN support status, and contributor workflow changes.

After publication, create clean environments that have no sibling checkout on
`PYTHONPATH`, install the published PySnpTools package on Python 3.10 and
3.14, and run the artifact smoke tests. FaST-LMM final qualification must use
this published version.

## Acceptance Criteria

The work is complete when:

- Metadata advertises Python 3.10 through 3.14.
- Python 3.14 resolves NumPy 2.3.5 or newer while older Python versions retain
  their intended lower bound.
- Every direct runtime dependency has a justified, tested lower bound, and
  packaging-only tools are not runtime dependencies.
- The direct dependency refresh checklist has a recorded result for both the
  declared boundary and current stable release of each dependency.
- Stable required dependencies resolve and install on every supported
  platform and Python version.
- BGEN support on Python 3.14 is either fully tested from stable artifacts or
  accurately excluded and documented.
- Issue #10 is fixed with immutable, valid tutorial-data references and a
  passing retrieval/integrity regression test.
- The committed lockfile supports the complete matrix; frozen, lower-bound,
  and scheduled freshness jobs pass.
- The complete unit, integration, and doctest coverage is known and passes.
- Normal `uv sync` and `uv run` commands emit no project-configuration
  deprecation warnings.
- The tutorial and paper notebooks execute from start to finish in a clean
  environment with an empty sample-data cache.
- The wheel and sdist contain the intended modules, metadata, licenses, and
  hashdown manifests.
- Direct-checkout and sdist-derived wheels have equivalent intended contents.
- Clean Python 3.10 and 3.14 environments install and exercise the exact
  artifacts successfully.
- Required CI passes on Linux, Windows, Intel macOS, and Apple Silicon macOS.
- Trusted Publishing releases the already-tested artifacts.
- The published version has a matching Git tag at its exact source commit.
- Release notes document support, dependency, data, optional-feature, and
  workflow changes.

## Out of Scope

- Dropping Python 3.10 through 3.13.
- Changing genetic-data formats or established numerical behavior.
- Rewriting PySnpTools or its dependencies in Rust.
- Modifying `bed-reader` without a demonstrated PySnpTools release blocker.
- FaST-LMM or fastlmmclib implementation work.
- Unrelated API redesign, performance work, or broad refactoring.
