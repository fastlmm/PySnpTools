# Releasing PySnpTools

This is the standing release process for PySnpTools. Short-term release plans
and compatibility investigations belong in `specs/`; keep this file limited to
the process that should apply to every release.

The maintainer chooses the version and explicitly approves publication. A
release must be built from a clean commit on `main`, and its `vX.Y.Z` tag must
point to that exact commit.

## Release order

PySnpTools must be published before a FaST-LMM release that depends on its new
behavior or compatibility. If a release requires a newer `bed-reader`, publish
and verify that dependency first.

Final qualification must use published dependencies, not sibling source
checkouts, Git dependencies, or prereleases.

## One-time repository setup

- Maintain the tag-triggered `.github/workflows/release.yml` workflow.
- Configure a PyPI Trusted Publisher for the `fastlmm/PySnpTools` repository,
  that workflow filename, and a protected GitHub `pypi` environment.
- Require maintainer approval for the `pypi` environment and restrict it to
  version tags.
- Give the publish job only the permissions it needs, including
  `id-token: write`. Other jobs need only `contents: read`.
- Pin third-party GitHub Actions to reviewed full commit SHAs and pin `uv` to a
  reviewed version.
- Keep a scheduled CI job that resolves current stable dependencies without
  rewriting or relying on the committed lockfile. Ordinary CI should continue
  to use the frozen lockfile for reproducibility.
- Pass the distributions built and tested by CI to the publish job as an
  artifact. Never rebuild them in the publish job.

The release workflow must publish with PyPI Trusted Publishing and retain the
generated attestations. Do not store PyPI passwords or API tokens in the
repository or GitHub Actions.

## Prepare the release

1. Start from a clean release branch based on current `main`.
2. Review open issues, pull requests, and dependency advisories for anything
   that affects the release.
3. Set the version in `pyproject.toml` and update `CHANGELOG.md` with the release
   date, compatibility changes, deprecations, and user-visible fixes.
4. Confirm that dependency bounds and Python-version markers describe versions
   actually tested in CI. When adding support for a new Python feature release,
   audit every direct dependency against both its declared lower bound and its
   current stable release. Record whether each bound was retained, raised, or
   given a Python-version marker, and test the resulting minimum-dependency
   environment rather than relying only on the locked current environment.
5. Regenerate and commit `uv.lock`, then verify it:

   ```console
   uv lock --check
   ```

6. Run the core and optional suites from the locked environment:

   ```console
   uv sync --frozen --all-extras
   uv run --frozen --no-sync python tests/test.py
   uv run --frozen --no-sync pytest pysnptools/distreader/test_bgen2.py
   ```

7. Execute the tutorial and paper notebooks from start to finish in a clean
   notebook environment. Inspect numerical results and committed output changes
   rather than accepting changed output mechanically.
8. Validate every entry in the package's hashdown manifests, including the
   packaged tutorial sample files.
9. Build the source distribution and wheel without local source overrides:

   ```console
   uv build --no-sources
   ```

10. Inspect both artifact manifests for required metadata, licenses, hashdown
    files, and package data. Confirm that development files, caches, notebooks,
    and generated output are absent unless intentionally distributed.
11. Install and smoke-test the exact wheel and source distribution in clean
    environments on the oldest and newest supported Python versions, outside
    the source checkout and without the repository on `PYTHONPATH`.
12. Wait for all required CI jobs to pass on every supported operating system
    and Python version. Resolve warnings that indicate a compatibility,
    packaging, or security problem.

## Publish

1. Merge the reviewed release change into `main` and confirm that `main` is
   clean and up to date.
2. Confirm that the version is not already present on PyPI and that the tag does
   not already exist.
3. Create and push an annotated tag at the release commit:

   ```console
   git tag -a vX.Y.Z -m "PySnpTools X.Y.Z"
   git push origin vX.Y.Z
   ```

4. Verify that the tag workflow builds and tests the exact artifacts intended
   for publication.
5. Review the workflow summary and approve its protected `pypi` environment.
6. Confirm that PyPI received both the wheel and source distribution and that
   their attestations are present.
7. Create the matching GitHub release from the changelog entry.

## Verify the published release

- Install `pysnptools==X.Y.Z` from PyPI into a fresh environment with no source
  checkout on `PYTHONPATH` and run an import and representative read/write
  smoke test.
- Re-run the clean artifact checks against the published distribution on the
  oldest and newest supported Python versions.
- Run the relevant InstallTest scenarios.
- Check the rendered documentation and tutorial links.
- Close or update the issues and pull requests resolved by the release, linking
  to the published version or its CI evidence.
- Continue with dependent releases only after these checks pass.

## Failed releases

Do not move or replace a published version tag, and do not upload different
files under an existing version. If a published release is unusable, yank it on
PyPI with a concise reason and publish a corrected version. Preserve the failed
release's tag and evidence for traceability.

TestPyPI may be used through a separate Trusted Publisher when a publishing
workflow itself needs qualification. It does not replace clean local artifact
tests or final testing against the real PyPI dependency ecosystem.
