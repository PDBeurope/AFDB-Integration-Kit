# ColabFold Complex End-to-End Reference

This directory is a small runnable reference that starts from curated
ColabFold-like complex fixture outputs and produces the AFDB integration
artifacts the toolkit is expected to emit.

Selected fixtures:

- Homodimer: `AF-0000000066074510` -> `Q46806` / `Q46806`
- Heterodimer: `AF-0000000300000101` -> `A0ABS2QMZ4` / `A0ABS2QMF5`

Input source:

- Raw fixture source files live under
  [`tests/fixtures/colabfold_real_examples/`](../../tests/fixtures/colabfold_real_examples/README.md).
- This example copies the selected fixture PDB and score JSON files into
  [`input/`](./input/) using normalized AFDB-style names:
  `*-model_v1.pdb` and `*-meta_v1.json`.

For the monomer reference example, see
[`examples/colabfold_monomer_e2e/`](../colabfold_monomer_e2e/README.md).

## Regenerate

Run the repo-owned helper from the repository root:

```bash
.venv/bin/python scripts/generate_colabfold_e2e_example.py \
  --duckdb examples/uniprot_example_subset.duckdb \
  --output-dir examples/colabfold_complex_e2e \
  --example-id AF-0000000066074510 \
  --example-id AF-0000000300000101 \
  --python-bin .venv/bin/python \
  --dssp-algorithm mkdssp \
  --bcif-backend biotite
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
7. Enrich the complex model and chain metadata with computed iPSAE/interface
   summary metrics from [`ipsae/ipsae_summary.csv`](./ipsae/ipsae_summary.csv).
8. Export ModelCIF generator input JSONs in [`modelcif_input/`](./modelcif_input/).
9. Generate ModelCIF files in [`modelcif/`](./modelcif/).
10. Run DSSP and write enriched mmCIF files in [`dssp/`](./dssp/).
11. Generate enriched PDB files in [`modelpdb/`](./modelpdb/).
12. Attempt BCIF generation and write results to [`bcif/`](./bcif/).

## Complex Metrics Source

- [`uniprot/templates/colabfold_example_modelcif_metadata.json`](../../uniprot/templates/colabfold_example_modelcif_metadata.json)
  is static scaffolding for ModelCIF metadata categories, provider text, and
  software parameter definitions.
- That template must not contain fake computed iPSAE values.
- The computed `complexPredictionAccuracy_*` values in the committed complex
  JSON outputs come from the iPSAE summary CSV and are injected into:
  - [`model_jsons/`](./model_jsons/) and [`model_batches/`](./model_batches/)
  - [`chain_jsons/`](./chain_jsons/) and [`chain_batches/`](./chain_batches/)
  - ModelCIF global QA metrics during ModelCIF generation
- The lightweight complex e2e example currently carries only the iPSAE-derived
  metric set. The separate clash/interface-analysis fields
  (`complexPredictionAccuracy_N_clash_backbone`,
  `complexPredictionAccuracy_N_clash_heavyAtom`, `numberOfInteractions`) are
  not part of this committed example unless a later stage adds them explicitly.

## Tooling Caveats

- Reproducing the metadata JSONs and ModelCIF inputs requires the committed
  DuckDB subset
  [`examples/uniprot_example_subset.duckdb`](../uniprot_example_subset.duckdb),
  which contains only the UniProt accessions used by the monomer and complex
  examples.
- `mkdssp` was available locally and was used for the DSSP step.
- Mol* `cif2bcif` was not on `PATH`, so the example uses the explicit
  `biotite` backend instead.
- The committed `.bcif` files were generated directly from the
  DSSP-enriched CIF files in [`dssp/`](./dssp/). The per-model status is
  recorded in [`run_summary.json`](./run_summary.json).
- `pydssp` is not installed in this environment, so it was not used.
- The heterodimer fixture metadata is exported as fragment chains using local
  chain ranges from the curated fixture set. This example does not reconstruct
  authoritative UniProt residue offsets beyond what is committed in the
  fixture corpus.
- Those local fragment ranges are preserved through the generated ModelCIF
  entity/reference alignment metadata and the heterodimer complex now gets an
  AFDB-style fallback complex name (`Complex of .../...`) in the JSON and PDB
  outputs.
- The committed complex metadata schemas require the computed iPSAE blocks on
  both per-record files and batch files, so removing those fields from a
  complex `model_jsons/*.json`, `model_batches/*.json`, `chain_jsons/*.json`,
  or `chain_batches/*.json` file will fail validation.
- Legacy PDB `DBREF` records remain a format caveat for long UniProt
  accessions. The mmCIF/ModelCIF and JSON artifacts are the authoritative
  metadata outputs for the heterodimer example.
