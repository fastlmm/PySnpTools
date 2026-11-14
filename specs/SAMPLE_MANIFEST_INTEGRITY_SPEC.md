# Sample manifest integrity audit

<!-- TODO0 consider deleting this spec once the work below is implemented and released. -->

## Findings on November 14, 2026

All 244 entries in the three packaged manifests were downloaded from their
immutable URLs and hashed (155,593,376 bytes). The BGEN and synthetic manifests
passed. Before this repair, the main PySnpTools manifest had 69 mismatches: 56 BIM files, 11 FAM
files, one Python example, and one NPZ file. The pre-repair manifest was byte-identical to
the one in the `v1.0.0` release; the retry change did not introduce these defects.

The main manifest pins PySnpTools commit
`ed14e050b2b75e7f4ddb73d512fbe928bbdb2b85`. Downloaded mismatches inspected
against that commit match the Git blobs, so these are stale manifest entries
rather than evidence of corrupted HTTP downloads.

For example, `pysnptools/examples/tutorial.py` expects MD5
`0d3bf3c2a2627cce0281ad0520201bea`, matching commit
`ef4f34e5cef881dabbc1f44ff6031a45b0e9ff09`. The pinned file has MD5
`5777c8ec12233d0468f47706689ebd29`; the source change is an equivalent chained
comparison. Genetic fixtures need more care: one historical BIM file contains
coordinate `122520001.0`, while the pinned file contains `131790273`.
These are not all whitespace or serialization-only changes.

All 69 old hashes match the original pinned snapshot,
`cf248cbf762516540470d693532590a77c76fba2`, used in May 2020. By July 2020,
44 of these entries disagreed with the newer source pin; by November 2023,
all 69 did. The pin advanced during fixture and Python 3 format updates while
these hashes retained older values. An initial path-history search missed
the original matching reader-name archive; checking the original pin found it.
Converting LF and CRLF does not explain any of the 69 mismatches on its own.

## Repair constraints

Further validation supports correcting the catalog against the already-pinned
files rather than restoring older fixture bytes:

- All 44 distributed BIM/FAM pieces match the sample IDs, marker IDs, order,
  and coordinates in the independently hashed `metadata.npz` exactly.
- The other 23 affected generated-data files match the current tracked
  fixtures exercised by the canonical suite exactly. Historical FAM changes
  include the deliberate `Fix gen iid` commit, not just whitespace. Historical
  BIM files also predate the addition of explicit allele columns.
- All three tracked versions of the reader-name NPZ archive contain exactly
  the same ordered 44 filenames; the archive bytes differ.
- The Python example's only relevant historical source change is the
  equivalent chained assertion described above.

These findings support correcting the catalog for the already-pinned dataset.
They do not justify changing fixture contents.

- Do not replace hashes merely because a download produced different bytes.
- Preserve sample/marker ordering, coordinates, phenotypes, and genotype data.
  Review the fixture set as a coherent dataset rather than independent files.
- Prefer restoring intended, independently verified immutable fixture content.
  A deliberate dataset migration needs its own numerical review and tests.
- Recheck every manifest entry, including currently unused samples, after any
  source pin or hash change. Preserve failing integrity diagnostics.
- Complete this repair before qualifying 1.0.1, as requested by the maintainer.

## Resolution

Corrected the 69 stale hashes without changing the URL, commit pin, file list,
or fixture contents. All 229 downloaded main-catalog files were compared byte
for byte with their pinned Git blobs. All 244 downloaded entries across the
three catalogs now pass MD5 verification.

Two new download regression tests failed with the original catalog and pass
with the corrected hashes. One compares a generated dataset's full genotype
matrix, IDs, and coordinates with the deterministic generator. The other
checks a distributed piece against the independently hashed metadata.
The coordinate comparison follows `Bed`'s documented representation of PLINK
zero positions as NaN; IDs and nonmissing coordinates are compared exactly.

## Evidence

Local audit details and downloaded files are under
`~/.cache/fastlmm-release-2026/pysnptools-1.0.1-manifests/`:
`audit.json` records all 244 results and 69 mismatches;
`historical-matches.json` records the 68 historical source commits.
These outputs are outside the repository and describe the original failures.
The historical source URLs and genetic fixtures were preserved by the repair.
