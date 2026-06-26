# Curated Real ColabFold Examples

This fixture set contains a small curated subset of real ColabFold outputs for
future Step 6 conversion and manifest tests.

Contents:

- `monomers/`: 3 single-chain ColabFold zip-derived examples.
- `homodimers/`: 3 two-chain examples copied from the local homodimer source
  directory.
- `heterodimers/`: 3 two-chain examples copied selectively from the large
  heterodimer tar archive.
- `manifest.json`: machine-readable validation summary for every copied
  example, including source path, copied filenames, sizes, pLDDT length, PAE
  dimension, chain spans derived from the copied PDB, and chain UniProt
  accessions.

Naming:

- Fixture directories and copied fixture files use a single AF-style
  `model_entity_id`.
- Monomer examples use reserved test fixture IDs:
  `AF-0000000300000001` through `AF-0000000300000003`.
- Homodimer examples keep their source AF IDs because those are already
  single-model IDs in the corrected merged ColabFold manifest.
- Heterodimer source files were named with two component AF IDs joined by an
  underscore. The fixture copies use reserved single test IDs:
  `AF-0000000300000101` through `AF-0000000300000103`. The original composite
  source ID and component reassigned AF IDs are preserved in `manifest.json`.

Fixture index:

| Fixture ID | Type | Source ID | Chains / accessions | Intended coverage |
| --- | --- | --- | --- | --- |
| `AF-0000000300000001` | Monomer | `ACATN_HUMAN_19de7` | A: `O00400` | Single-chain raw ColabFold conversion with long pLDDT/PAE arrays. |
| `AF-0000000300000002` | Monomer | `C76C2_ARATH_6db51` | A: `O64637` | Single-chain conversion with a different AlphaFold2 model rank/source. |
| `AF-0000000300000003` | Monomer | `CDK9_CAEEL_5ca86` | A: `Q9TVL3` | Single-chain conversion and monomer input CSV provenance. |
| `AF-0000000065760001` | Homodimer | `AF-0000000065760001` | A/B: `Q6GZX4` | Two-chain homomer spans with raw `meta_v1` confidence JSON. |
| `AF-0000000066074510` | Homodimer | `AF-0000000066074510` | A/B: `Q46806` | Higher-quality two-chain homomer fixture with raw `meta_v1` confidence JSON. |
| `AF-0000000065760002` | Homodimer | `AF-0000000065760002` | A/B: `Q6GZX3` | Larger homomer confidence and chain-span coverage. |
| `AF-0000000065760003` | Homodimer | `AF-0000000065760003` | A/B: `Q197F8` | Largest curated homomer fixture in this set. |
| `AF-0000000300000101` | Heterodimer | `AF_0000000066426974_AF_0000000066426875` | A: `A0ABS2QMZ4`; B: `A0ABS2QMF5` | Real heteromer with normalized single fixture ID and preserved component IDs. |
| `AF-0000000300000102` | Heterodimer | `AF_0000000066576660_AF_0000000066577761` | A: `A0ABS8Y874`; B: `A0ABS8RL98` | Real heteromer chain-accession mapping coverage. |
| `AF-0000000300000103` | Heterodimer | `AF_0000000066908168_AF_0000000066909825` | A: `A0AAV6YE70`; B: `A0AAV6YAU6` | Real heteromer conversion coverage with short chains. |

Selection rules used here:

- Keep only the minimal files needed for later conversion tests.
- Prefer `rank_001` PDB plus the matching raw ColabFold score JSON.
- For monomers, also keep the tiny input `.csv` file because it preserves the
  original job/example identifier and sequence.
- Do not copy MSA files, PNGs, environment folders, tarballs, logs, or extra
  ranked models.
- Normalize copied filenames to:
  - `<model_entity_id>-model_v1.pdb`
  - `<model_entity_id>-scores_v1.json`
  - `<model_entity_id>-input.csv` for monomer input CSVs

Validation performed before copying:

- `plddt`, `pae`, and `max_pae` are present in the selected JSON.
- PDB residue count equals `len(plddt)`.
- PAE is square and its dimension equals `len(plddt)`.
- Chain IDs and flattened residue spans are recorded in `manifest.json`.

How tests should use this set:

- Read `manifest.json` for expected metadata instead of inferring chain
  accessions or source provenance from filenames.
- Use these fixtures for integration-style conversion and manifest checks.
  Keep synthetic tests for small edge cases that are easier to isolate.
- Treat these as curated regression fixtures, not as a refreshed example corpus.
  A broad metadata refresh should be a separate repo-shape or data-preparation
  task.

Notes:

- Homodimer source JSON files are named `*-meta_v1.json`, but they contain the
  required raw ColabFold-like fields (`plddt`, `pae`, `max_pae`).
- Chain UniProt accessions are recorded from the corrected merged ColabFold
  manifest and related corrected mapping outputs under
  `/mnt/disks/toolkit-data/mapping_file_heterodimer/afid_overlap_remediation_bundle_2026-02-24/`.
- Heterodimer examples were extracted only temporarily under `/tmp`; only the
  selected files were copied into the repository.
