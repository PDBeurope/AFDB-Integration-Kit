# Examples Overview

This directory contains the committed end-to-end reference examples for raw
ColabFold-like outputs -> AFDB integration artifacts.

## Recommended References

The two committed runnable references are:

- [`examples/colabfold_monomer_e2e/`](./colabfold_monomer_e2e/): three curated
  monomer fixtures
- [`examples/colabfold_complex_e2e/`](./colabfold_complex_e2e/): one homodimer
  and one heterodimer fixture

There is also a fragment-metadata reference:

- [`examples/fragment_metadata_e2e/`](./fragment_metadata_e2e/): validates
  custom fragment names across seven real collaborator monomer/dimer
  topologies. Its committed runner stops before coordinates; a companion
  deterministic synthetic-coordinate runner continues through final ModelCIF.

The two coordinate references both include:

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

For the complex reference, the committed `complexPredictionAccuracy_*` fields
in the JSON metadata are computed outputs from the iPSAE enrichment stage, not
template literals. The ModelCIF template under
[`uniprot/templates/colabfold_example_modelcif_metadata.json`](../uniprot/templates/colabfold_example_modelcif_metadata.json)
only supplies static ModelCIF scaffolding plus software parameter definitions;
the computed complex metrics are injected from the enriched complex model JSONs
before ModelCIF export.

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
  --example-id AF-0000000066074510 \
  --example-id AF-0000000300000101
```

See [`examples/colabfold_monomer_e2e/README.md`](./colabfold_monomer_e2e/README.md)
and [`examples/colabfold_complex_e2e/README.md`](./colabfold_complex_e2e/README.md)
for the selected fixtures, caveats, and generated file layout.

Regenerate the fragment-metadata pre-structure reference with:

```bash
.venv/bin/python scripts/generate_fragment_metadata_e2e_example.py
```

It recreates `examples/fragment_metadata_e2e/generated/` and validates the
DuckDB-to-manifest-to-metadata path only. See
[`examples/fragment_metadata_e2e/README.md`](./fragment_metadata_e2e/README.md)
for its coverage.

When the local donor pairs are available under
`/mnt/disks/toolkit-data/viruses/sample_data`, reproduce the synthetic
coordinate continuation with:

```bash
.venv/bin/python scripts/generate_fragment_coordinate_e2e_example.py
```

It writes to
`/mnt/disks/toolkit-data/viruses/fragment_metadata_synthetic_e2e/generated/`
and validates seven PDB/meta pairs, converted scores, individual and batched
metadata, ModelCIF input, and seven final ModelCIF files. The assets are
deterministic software fixtures and are explicitly unsuitable for scientific
interpretation.

## Validate Committed E2E Outputs

Run these checks from the repo root after installing the project environment.
They validate one representative monomer and one representative complex from
the committed reference trees.

The fragment-metadata reference has no committed coordinate artifacts. Run
the appropriate pre-structure or external synthetic-coordinate command above.

### Metadata JSONs

Use `run-schema-validation` with the dedicated example-output schemas for
committed `model_jsons/*.json`, `chain_jsons/*.json`, and
`config/provider.json` files, and use the same validator on the committed batch
files:

```bash
.venv/bin/python main.py run-schema-validation \
  -i examples/colabfold_monomer_e2e/model_jsons/AF-0000000300000001.json \
  -t model-summary
.venv/bin/python main.py run-schema-validation \
  -i examples/colabfold_monomer_e2e/model_batches/AF-metadata-1-of-1.json \
  -t model-summary
.venv/bin/python main.py run-schema-validation \
  -i examples/colabfold_monomer_e2e/chain_jsons/AF-0000000300000001.json \
  -t collection-doc
.venv/bin/python main.py run-schema-validation \
  -i examples/colabfold_monomer_e2e/chain_batches/AF-chain-metadata-1-of-1.json \
  -t collection-doc
.venv/bin/python main.py run-schema-validation \
  -i examples/colabfold_monomer_e2e/config/provider.json \
  -t provider

.venv/bin/python main.py run-schema-validation \
  -i examples/colabfold_complex_e2e/model_jsons/AF-0000000300000101.json \
  -t model-summary
.venv/bin/python main.py run-schema-validation \
  -i examples/colabfold_complex_e2e/model_batches/AF-metadata-1-of-1.json \
  -t model-summary
.venv/bin/python main.py run-schema-validation \
  -i examples/colabfold_complex_e2e/chain_jsons/AF-0000000300000101.json \
  -t collection-doc
.venv/bin/python main.py run-schema-validation \
  -i examples/colabfold_complex_e2e/chain_batches/AF-chain-metadata-1-of-1.json \
  -t collection-doc
.venv/bin/python main.py run-schema-validation \
  -i examples/colabfold_complex_e2e/config/provider.json \
  -t provider
```

For these committed e2e references, use `model-summary` for `model_jsons/*.json`
and `model_batches/*.json`, `collection-doc` for `chain_jsons/*.json` and
`chain_batches/*.json`, and `provider` for `config/provider.json`. The
canonical `model` schema remains for full model metadata entries only.

The `validate-metadata-file` command uses the same shared metadata schema
validator as `run-schema-validation` and requires the same explicit `--type`
value, for example:

```bash
.venv/bin/python main.py validate-metadata-file \
  --file examples/colabfold_monomer_e2e/model_jsons/AF-0000000300000001.json \
  --type model-summary
```

### Score JSONs

Validate the confidence JSON, the PAE JSON, and the confidence/PAE length
relationship:

```bash
.venv/bin/python main.py validate-plddt-file \
  --file examples/colabfold_monomer_e2e/scores/AF-0000000300000001-confidence_v1.json
.venv/bin/python main.py validate-pae-file \
  --file examples/colabfold_monomer_e2e/scores/AF-0000000300000001-predicted_aligned_error_v1.json
.venv/bin/python main.py validate-relationships-pair \
  --plddt-file examples/colabfold_monomer_e2e/scores/AF-0000000300000001-confidence_v1.json \
  --pae-file examples/colabfold_monomer_e2e/scores/AF-0000000300000001-predicted_aligned_error_v1.json

.venv/bin/python main.py validate-plddt-file \
  --file examples/colabfold_complex_e2e/scores/AF-0000000300000101-confidence_v1.json
.venv/bin/python main.py validate-pae-file \
  --file examples/colabfold_complex_e2e/scores/AF-0000000300000101-predicted_aligned_error_v1.json
.venv/bin/python main.py validate-relationships-pair \
  --plddt-file examples/colabfold_complex_e2e/scores/AF-0000000300000101-confidence_v1.json \
  --pae-file examples/colabfold_complex_e2e/scores/AF-0000000300000101-predicted_aligned_error_v1.json
```

### ModelCIF Dictionary Validation

Validate representative ModelCIF files with `gemmi` and the ModelCIF
dictionary:

```bash
gemmi validate -p -d mmcif_ma.dic \
  examples/colabfold_monomer_e2e/modelcif/AF-0000000300000001-model_v1.cif
gemmi validate -p -d mmcif_ma.dic \
  examples/colabfold_complex_e2e/modelcif/AF-0000000300000101-model_v1.cif
```

### Manual Coordinate-File Viewer Checks

Open representative PDB, ModelCIF, and BCIF files in the Mol* web viewer:

1. Go to https://molstar.org/viewer/.
2. Drag and drop the file into the browser window, or use **Open Files** in
   the left panel.
3. Confirm the structure opens correctly, no error messages are shown in the
   viewer, and the structure looks structurally correct by eye.

Representative files:

```text
examples/colabfold_monomer_e2e/modelpdb/AF-0000000300000001-model_v1.pdb
examples/colabfold_monomer_e2e/modelcif/AF-0000000300000001-model_v1.cif
examples/colabfold_monomer_e2e/bcif/AF-0000000300000001-model_v1.bcif

examples/colabfold_complex_e2e/modelpdb/AF-0000000300000101-model_v1.pdb
examples/colabfold_complex_e2e/modelcif/AF-0000000300000101-model_v1.cif
examples/colabfold_complex_e2e/bcif/AF-0000000300000101-model_v1.bcif
```

Optionally open the same files in ChimeraX or another preferred structure
viewer such as PyMOL. The expected result is a clean import with no
parser/import errors and a structure that looks correct by eye.
