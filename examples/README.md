# Examples Overview

This directory contains the committed end-to-end reference examples for raw
ColabFold-like outputs -> AFDB integration artifacts.

## Recommended References

The two committed runnable references are:

- [`examples/colabfold_monomer_e2e/`](./colabfold_monomer_e2e/): three curated
  monomer fixtures
- [`examples/colabfold_complex_e2e/`](./colabfold_complex_e2e/): one homodimer
  and one heterodimer fixture

Both references include:

- normalized `*-meta_v1.json` + `*-model_v1.pdb` inputs
- converted confidence and PAE JSONs
- merged manifests
- per-model and per-chain metadata JSONs plus batch files
- ModelCIF input JSONs
- generated ModelCIF files
- DSSP-enriched mmCIF files
- enriched PDB files
- BCIF outputs with backend/fallback notes
- exact regeneration commands in `config/commands.txt`

The committed DuckDB subset
[`examples/uniprot_example_subset.duckdb`](./uniprot_example_subset.duckdb)
contains only the six UniProt accessions needed by the monomer and complex
references, so the examples no longer depend on a VM-local UniProt database.

Regenerate the monomer reference from the repo root with:

```bash
.venv/bin/python scripts/generate_colabfold_e2e_example.py \
  --duckdb examples/uniprot_example_subset.duckdb \
  --output-dir examples/colabfold_monomer_e2e
```

Regenerate the complex reference with:

```bash
.venv/bin/python scripts/generate_colabfold_e2e_example.py \
  --duckdb examples/uniprot_example_subset.duckdb \
  --output-dir examples/colabfold_complex_e2e \
  --example-id AF-0000000065760001 \
  --example-id AF-0000000300000101
```

See [`examples/colabfold_monomer_e2e/README.md`](./colabfold_monomer_e2e/README.md)
and [`examples/colabfold_complex_e2e/README.md`](./colabfold_complex_e2e/README.md)
for the selected fixtures, caveats, and generated file layout.
