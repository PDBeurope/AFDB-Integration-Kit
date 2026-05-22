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

Notes:

- Homodimer source JSON files are named `*-meta_v1.json`, but they contain the
  required raw ColabFold-like fields (`plddt`, `pae`, `max_pae`).
- Chain UniProt accessions are recorded from the corrected merged ColabFold
  manifest and related corrected mapping outputs under
  `/mnt/disks/toolkit-data/mapping_file_heterodimer/afid_overlap_remediation_bundle_2026-02-24/`.
- Heterodimer examples were extracted only temporarily under `/tmp`; only the
  selected files were copied into the repository.
