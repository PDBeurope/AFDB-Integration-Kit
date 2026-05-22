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
  dimension, and chain spans derived from the copied PDB.

Selection rules used here:

- Keep only the minimal files needed for later conversion tests.
- Prefer `rank_001` PDB plus the matching raw ColabFold score JSON.
- For monomers, also keep the tiny input `.csv` file because it preserves the
  original job/example identifier and sequence.
- Do not copy MSA files, PNGs, environment folders, tarballs, logs, or extra
  ranked models.

Validation performed before copying:

- `plddt`, `pae`, and `max_pae` are present in the selected JSON.
- PDB residue count equals `len(plddt)`.
- PAE is square and its dimension equals `len(plddt)`.
- Chain IDs and flattened residue spans are recorded in `manifest.json`.

Notes:

- Homodimer source JSON files are named `*-meta_v1.json`, but they contain the
  required raw ColabFold-like fields (`plddt`, `pae`, `max_pae`).
- Heterodimer examples were extracted only temporarily under `/tmp`; only the
  selected files were copied into the repository.
