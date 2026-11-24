<!-- TODO0 consider deleting this spec once the work below is implemented and released. -->

# PySnpTools 1.1.0b1: handoff

Written November 23, 2026, for the agent taking over. Read `AGENTS.md` here
and in FaST-LMM first; Carl's personal rules are in `~/.codex/AGENTS.md`.
Never publish, approve a PyPI deployment, push a release tag, or merge to
`main` yourself: prepare everything and give Carl the exact commands.

## Where things stand

- Branch `cupy-sub-matrix`, open as PR #25 against `main`.
- Last pushed commit: `d48b746`. CI run 45313830400 for it was in progress at
  handoff; the two runs before it (on `bba4f6b` and `89db1f1`) passed 39 of 39.
- Committed locally, **not yet pushed**: `80fac2d` (CI runs the pytest-only
  test modules) and this spec. CI uses `cancel-in-progress: true`, so push
  only after the current run finishes, or it cancels that run. Then monitor
  the new run in the background and report each job (see Carl's CI rule in
  `~/.codex/AGENTS.md`).
- Version `1.1.0b1`; requires `bed-reader[samples]>=1.3.0b3` (published on
  PyPI and crates.io). `uv.lock` is refreshed against the published b3.
- Full local validation passed against the published b3 (legacy runner, full
  pytest, GPU tests, BGEN, strict docs, wheel and sdist installed in clean
  venvs on Python 3.11 and 3.14 with CuPy and 3.15 without). Later commits
  passed the legacy runner, full pytest, GPU tests, and docs locally, but
  the full artifact validation must be **rerun on the final head** (see
  below).

## Carl's review against `main`

Carl reviews `git diff origin/main...cupy-sub-matrix` in eight batches. For
each batch: the agent reviews first, fixes findings (with regression tests),
then gives Carl a short file list with clickable line links and what to look
for. Carl answers in small steps ("a bit at a time"); keep replies short.
Progress is tracked in FaST-LMM `specs/REVIEW_INVENTORY.md` (section "Carl's
PySnpTools review against `main`").

| Batch | Status |
| --- | --- |
| 1. GPU and precision core | Approved (`69112e4`; delta at `9ae0a1d`) |
| 2. Standardizers | Approved (`077f35b`) |
| 3. Array utilities and containers | Approved (`30ebc30`) |
| 4. SNP readers | Approved (`bba4f6b`) |
| 5. Kernels, distributed readers, writers | Approved (`89db1f1`) |
| 6. Thread control | Approved (`d48b746`) |
| 7. Tests | Agent review done; CI gap fixed in `80fac2d`. **Waiting on Carl.** |
| 8. Docs and packaging | Open: `CHANGELOG.md`, `README.md`, `doc/`, `pyproject.toml`, `uv.lock`, `.github/`, `specs/` |

After each approval, update the batch row in FaST-LMM
`specs/REVIEW_INVENTORY.md` (commit there too).

### Batch 8 checklist for the agent review

- `CHANGELOG.md` 1.1.0b1 entry: complete and accurate against the diff
  (strict `xp=None` rule, writers copy to NumPy, pickling, Beta and
  `UnitTrained` fixes, explicit `block_size` rounding change, `block_size`
  warning, constant fractional SNPs, b3 requirement, withdrawn 1.1.0a1 GPU
  features). The date line says 2026-11-23; update it to the release date.
- `doc/source/gpu.rst`, `index.rst`, `gpu-api-preview.rst`, the notebook
  `doc/ipynb/gpu_workflow.ipynb`, and `README.md`: no stale bed-reader GPU,
  1.3.0a1, or "silently copies" statements; examples obey the strict rule.
- `pyproject.toml`, `uv.lock`, `.github/workflows/{ci,release}.yml`.
- `specs/`: `GPU_DOCUMENTATION_AUDIT.md` and
  `SAMPLE_MANIFEST_INTEGRITY_SPEC.md` carry "consider deleting" TODO0
  markers; ask Carl whether to delete them before release.

## Decisions Carl made (do not reopen)

- `warn_gpu_float64` stays private (`pysnptools.util._gpu_warnings`).
- Strict `xp=None`: `xp=None` means `ARRAY_MODULE`, then NumPy. An existing
  NumPy or CuPy array is never copied to the other module implicitly:
  constructors (`SnpData`, `KernelData`, `DistData`, `PstData`), `.val`
  assignment, and reads or kernels of in-memory values raise `ValueError`.
  Lists, `val=None`, and file-decoded values (including text readers and
  memory maps) follow `ARRAY_MODULE`. Core: `pysnptools/util/_device.py`.
- Writers accept CuPy-backed data and copy it to NumPy explicitly
  (`_device.numpy_for_write`).
- `block_size` for standardization stays accepted with a
  `DeprecationWarning` (it has been a no-op since 2015); not removed.
  `read_kernel`'s `block_size` is a different, live feature.

## Remaining steps, in order

1. Push `80fac2d` and this spec after CI run 45313830400 finishes; monitor CI.
2. Get Carl's batch 7 approval; do the batch 8 agent review, fix, present.
3. Rerun the full validation on the final head:
   `bash /mnt/c/deldir/pysnptools-1.1.0b1/validate.sh source` then
   `... validate.sh artifacts` (logs in `/mnt/c/deldir/pysnptools-1.1.0b1/logs/`).
   Confirm the logged commit is the final head and the tree is clean.
4. Release gate (`AGENTS.md`, `RELEASING.md` "Publish"): re-read FaST-LMM
   `specs/todo.md` (shared), this `specs/` folder, open issues and PRs
   (`gh issue list`, `gh pr list`), and `git grep -niE "todo0|cmk0"`. List
   every item for Carl as included or deferred, with the reason. At handoff:
   issue #26 is fixed on this branch (PR #27 is superseded: recommend closing
   it, and closing #26 when #25 merges); PR #24 (Dependabot `uv_build`
   0.12.24 to 0.13.7) was deferred as a separate build-backend change.
5. Give Carl the publish commands (RELEASING.md requires a clean commit on
   `main`): merge PR #25 (Carl's action), confirm `main` CI, confirm
   1.1.0b1 is not on PyPI, then
   `git tag -a v1.1.0b1 -m "PySnpTools 1.1.0b1"` and
   `git push origin v1.1.0b1`; Carl approves the `pypi` environment and
   creates the GitHub release from the changelog.
6. After Carl publishes, remind him to consider yanking PySnpTools 1.1.0a1
   (its unbounded `bed-reader>=1.3.0a1` and `bed-reader[cupy]` extra no
   longer work with the CPU-only bed-reader).

## FaST-LMM side (Part 2 of the original request)

Branch `gpu-bed-reader`, pushed at `1a4c88d`. Already done there:

- `a958625`: the CPU-only predictor methods read covariates and phenotypes
  with `xp=np` (found by running FaST-LMM under `ARRAY_MODULE=cupy` against
  the strict rule).
- `11833d2`: stop passing the ignored `block_size` to `SnpData.standardize`.
- FaST-LMM's full suite (230 tests) passed against this PySnpTools with
  `ARRAY_MODULE` unset; under `ARRAY_MODULE=cupy` the remaining failures are
  FaST-LMM test code that builds NumPy containers or assigns NumPy to GPU
  `.val` (intended under the strict rule; several also failed before).

After PySnpTools 1.1.0b1 is on PyPI, Part 2 (from Carl's original request):

1. Raise bounds to `bed-reader>=1.3.0b3` and `pysnptools>=1.1.0b1` in both
   `dependencies` and the `cupy` extra; relock; bump 1.1.0a1 to 1.1.0b1.
2. `RELEASE_NOTES.md`: include the CPU speedups, not only the GPU. Numbers in
   `specs/SYNTHETIC_BENCHMARK.md` (November 22, 14,925 x 100,000): GPU
   float32 33.2 s; CPU float32 158.9 s; CPU float64 395.8 s (both CPU at 16
   threads). CPU float32 is 2.5x faster than float64 because the kernel is
   built in float32 and decompositions use SciPy float32 LAPACK (NumPy
   silently upcast float32 to float64). 16 threads beat 32 on this 16-core
   machine. Top-10 SNP ranking identical across configurations; max
   |log10 P| difference vs CPU float64 is 1.4e-3 (CPU f32) and 3.4e-3 (GPU
   f32). Mention the strict `xp=None` behavior inherited from PySnpTools.
3. `specs/todo.md`: move finished bed-reader items to Done (1.3.0b2, 1.3.0b3,
   notes on issues #44 and #45), remove the stale `b535ef3` entry and the
   old 12.7x / 384.8 s headline, update "alpha" wording to beta. No TODO
   markers in `todo.md`; never touch its opening `todo000` line.
4. Final checks against the published bed-reader and PySnpTools wheels from
   PyPI (not sibling checkouts): clean-venv wheel and sdist installs, the
   full suite (run `tests/test.py` from the `tests/` folder, as CI does), GPU
   tests. Then start CI and monitor it.

## Pitfalls learned in this work

- Ask before any CPU- or GPU-heavy benchmark; normal test suites are fine.
- Scripts live in `/mnt/c/deldir/pysnptools-1.1.0b1/` (`validate.sh`,
  `b3-local.sh`, `cupy-env-compare.sh`, `fastlmm-cupy-env-compare.sh`,
  `readers_under_cupy.py`); venvs in `~/.cache/pysnptools-1.1.0b1/` and
  `~/.cache/fastlmm-device-check/`.
- Export `QCTOOLPATH=$HOME/.local/bin/qctool` or two BGEN doctests fail.
- `pysnptools/test.py` belongs to the legacy runner (`tests/test.py`); with
  `pytest -k` it hits a collection error, so run single tests by node ID.
- Use `env -u ARRAY_MODULE` or `unset`, never `ARRAY_MODULE=` (an empty value
  raises `ValueError`).
- Test runs can delete or create files in the repositories (FaST-LMM
  `tmp/feature_selection/`, PySnpTools `distreader/localcache1/` and
  `peertopeer1/`). Check `git status` before committing; never `git add -A`.
- Some files use CRLF (FaST-LMM `fastlmm/inference/*.py`, PySnpTools
  `snpreader/dense.py`); edit them byte-wise to preserve line endings.
- No AI attribution in commits, PRs, or release notes; suggest a 1-2 line
  commit message for each commit; write "PySnpTools" in full.
