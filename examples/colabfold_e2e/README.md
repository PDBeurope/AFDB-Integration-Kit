# ColabFold End-to-End Reference

This directory is a small runnable reference that starts from curated
ColabFold-like fixture outputs and produces the AFDB integration artifacts the
toolkit is expected to emit.

Selected fixtures:

- `AF-0000000300000001` -> `O00400`
- `AF-0000000300000002` -> `O64637`
- `AF-0000000300000003` -> `Q9TVL3`

Input source:

- Raw fixture source files live under
  [`tests/fixtures/colabfold_real_examples/`](../../tests/fixtures/colabfold_real_examples/README.md).
- This example copies the selected fixture PDB and score JSON files into
  [`input/`](./input/) using normalized AFDB-style names:
  `*-model_v1.pdb` and `*-meta_v1.json`.

## Regenerate

Run the repo-owned helper from the repository root:

```bash
.venv/bin/python scripts/generate_colabfold_e2e_example.py \
  --duckdb /mnt/disks/toolkit-data/uniprot_extract_2025_04_merged_5way/db/uniprot_2025_04_merged_5way.duckdb \
  --output-dir examples/colabfold_e2e
```

That helper executes the same script sequence used to populate this directory
and records the exact subprocess commands in
[`config/commands.txt`](./config/commands.txt).

## Flow

The generated files follow the old Nextflow end-to-end order:

1. Normalize raw ColabFold fixture inputs into [`input/`](./input/).
2. Validate PDB/score consistency via [`validation_results.tsv`](./validation_results.tsv).
3. Convert ColabFold scores to AFDB confidence + PAE JSONs in [`scores/`](./scores/).
4. Write per-model manifests in [`chain_manifests/`](./chain_manifests/) and [`model_manifests/`](./model_manifests/), then merge them into [`merged_manifests/`](./merged_manifests/).
5. Export per-model metadata JSONs in [`model_jsons/`](./model_jsons/) and batch them into [`model_batches/`](./model_batches/).
6. Export per-chain metadata JSONs in [`chain_jsons/`](./chain_jsons/) and batch them into [`chain_batches/`](./chain_batches/).
7. Export ModelCIF generator input JSONs in [`modelcif_input/`](./modelcif_input/).
8. Generate ModelCIF files in [`modelcif/`](./modelcif/).
9. Run DSSP and write enriched mmCIF files in [`dssp/`](./dssp/).
10. Generate enriched PDB files in [`modelpdb/`](./modelpdb/).
11. Attempt BCIF generation and write results to [`bcif/`](./bcif/).

## Tooling Caveats

- Reproducing the metadata JSONs and ModelCIF inputs requires the DuckDB at
  `/mnt/disks/toolkit-data/uniprot_extract_2025_04_merged_5way/db/uniprot_2025_04_merged_5way.duckdb`.
- `mkdssp` was available locally and was used for the DSSP step.
- Mol* `cif2bcif` was not on `PATH`, so the example uses the explicit
  `biotite` backend instead.
- The local `biotite` backend could not convert the DSSP-enriched CIF files.
  Each `.bcif` file in [`bcif/`](./bcif/) was therefore generated from the
  pre-DSSP ModelCIF file instead. The per-model status is recorded in
  [`run_summary.json`](./run_summary.json).
- `pydssp` is not installed in this environment, so it was not used.
- iPSAE, clash, and interface analysis are intentionally not part of this
  Step 9 reference because the pivot target is the raw ColabFold -> final
  AFDB artifact path.
