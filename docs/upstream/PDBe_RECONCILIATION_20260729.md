# PDBe upstream reconciliation — 2026-07-29

## Scope

This reconciliation advances the reviewed PDBe upstream boundary from
`dcc5e034c9aa64ab1d9405aeb97db81ac2acaad9` through
`722bf67` on `PDBeurope/AFDB-Integration-Kit` `main`.

## Commit disposition

| Upstream commit | Disposition | AFCDB treatment |
| --- | --- | --- |
| `314e714` Allow absent gene metadata | Ported | Gene metadata is optional in model-summary and collection schemas; when present, all three metadata schemas require a non-empty, non-whitespace string. Embedded validation schemas were updated for parity. |
| `98c49c7` Normalize ModelCIF provenance metadata | Ported selectively | Added canonical AlphaFold/AlphaFold-Multimer, ipSAE, DSSP/PyDSSP provenance normalization and wired it into generation, exporters, and CLI paths. AFCDB-specific pipeline structure was retained. |
| `065d4d3` Refresh ModelCIF e2e provenance examples | Not ported | The commit refreshes generated example outputs. Generated examples are outside the AFCDB runtime reconciliation boundary and the referenced ColabFold example tree is not present in this branch. Equivalent runtime behavior is covered by focused tests. |
| `722bf67` Merge branch `vad-modelcif-provenance-review` | No unique delta | Merge-only commit; its constituent changes are accounted for above. |

## Verification

- `tests/test_validator.py`: absent gene accepted; null, empty, and whitespace placeholders rejected.
- `tests/test_modelcif_generate.py`: monomer and complex provenance normalization, strict AlphaFold version handling, and generator integration.
- Full test suite: `324 passed, 1 skipped`.

The reconciliation remains selective: upstream generated artifacts and repository
layout changes are not imported unless they affect AFCDB runtime behavior.
