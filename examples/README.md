# Examples Overview

This directory now has two roles:

- [`colabfold_e2e/`](./colabfold_e2e/): the primary Step 9 runnable reference
  for raw ColabFold-like outputs -> AFDB integration artifacts.
- legacy root example files such as `AF-0000000000000001-*`: older hand-curated
  sample assets kept for backward compatibility and ad hoc manual checks.

## Recommended Reference

Use [`examples/colabfold_e2e/`](./colabfold_e2e/) when you want the current
copy-pasteable workflow. It includes:

- normalized `*-meta_v1.json` + `*-model_v1.pdb` inputs for three curated
  monomer fixtures,
- converted confidence and PAE JSONs,
- merged manifests,
- per-model and per-chain metadata JSONs plus batch files,
- ModelCIF input JSONs,
- generated ModelCIF files,
- DSSP-enriched mmCIF files,
- enriched PDB files,
- BCIF outputs with backend/fallback notes,
- exact regeneration commands in `config/commands.txt`.

Regenerate it from the repo root with:

```bash
.venv/bin/python scripts/generate_colabfold_e2e_example.py \
  --duckdb /mnt/disks/toolkit-data/uniprot_extract_2025_04_merged_5way/db/uniprot_2025_04_merged_5way.duckdb \
  --output-dir examples/colabfold_e2e
```

See [`examples/colabfold_e2e/README.md`](./colabfold_e2e/README.md) for the
selected fixtures, caveats, and generated file layout.
