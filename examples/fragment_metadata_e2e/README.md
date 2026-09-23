# Fragment Metadata E2E Reference

This is a small, runnable pre-structure reference for the custom fragment-name
workflow. It proves that a name assigned to a UniProt residue range reaches the
canonical chain manifest, model metadata, chain metadata, and ModelCIF input
for the topology cases we expect to support. It deliberately stops before any
coordinate-dependent work. A companion synthetic-coordinate run now carries
the same seven cases through final ModelCIF generation.

## Coverage

| Model ID | Topology | Source |
| --- | --- | --- |
| `AF-0000000211971324` | Full-length monomer | Real collaborator AF-ID/topology; coordinates unavailable |
| `AF-0000000212005744` | Fragment monomer | Real collaborator AF-ID/topology; coordinates unavailable |
| `AF-0000000212013504` | Full-length homodimer | Real collaborator AF-ID/topology; coordinates unavailable |
| `AF-0000000212013519` | Full-length heterodimer | Real collaborator AF-ID/topology; coordinates unavailable |
| `AF-0000000212039399` | Fragment homodimer | Real collaborator AF-ID/topology; coordinates unavailable |
| `AF-0000000212039401` | Fragment heterodimer | Real collaborator AF-ID/topology; coordinates unavailable |
| `AF-0000000212039400` | Mixed fragment/full-length heterodimer | Real collaborator AF-ID/topology; coordinates unavailable |

`AF-0000000212039400` is the important mixed case: chain A is the `P27409`
fragment `1..46`, while chain B is the full-length `P28711` protein.

The seven real rows are copied exactly into
[`config/source_collaborator_wide.csv`](./config/source_collaborator_wide.csv)
from the collaborator-style wide manifest. The single 12-row
[`config/canonical_input_manifest.csv`](./config/canonical_input_manifest.csv)
is the A/B manifest used by the downstream checks. Collaborator `chain_a_id` and
`chain_b_id` values are occupancy/source identifiers, not authoritative protein
names; this workflow never derives names from them.

The names in [`config/fragment_metadata.json`](./config/fragment_metadata.json)
and all UniProt names/sequences represented by
[`config/mock_uniprot_seed.json`](./config/mock_uniprot_seed.json) are invented
development data. They demonstrate mechanics only and are not biological
annotations to reuse in production.

## Run

From the repository root, run:

```bash
.venv/bin/python scripts/generate_fragment_metadata_e2e_example.py
```

The script recreates only `examples/fragment_metadata_e2e/generated/`, builds a
small mock DuckDB, loads fragment names, normalizes the real wide rows,
enriches the 12-row canonical manifest, runs the existing metadata and
ModelCIF-input exporters, combines metadata batches, and validates the result.

The generated layout is:

```text
generated/
  mock_uniprot.duckdb
  config/                 # derived dataset/provider configuration and model IDs
  logs/stages.log          # exact commands and command output
  reports/                 # loader and manifest-enrichment reports
  manifests/               # wide-derived and canonical enriched manifests
  model_metadata/          # one JSON record per model
  chain_metadata/          # one JSON collection per model
  modelcif_input/          # one ModelCIF generator input per model
  model_batches/           # combined model records
  chain_batches/           # combined chain records
  run_summary.json
```

Successful completion means there are seven models and 12 chain rows, the two
custom fragment names have been joined by accession plus range, repeated
components share an entity, distinct fragments of one accession do not, and
every consumer carries the expected name. Useful inspection points are:

```bash
rg 'Development fragment' \
  examples/fragment_metadata_e2e/generated/manifests/canonical_enriched.csv
rg 'Development fragment' \
  examples/fragment_metadata_e2e/generated/model_metadata/
rg 'Development fragment' \
  examples/fragment_metadata_e2e/generated/chain_metadata/
rg 'Development fragment' \
  examples/fragment_metadata_e2e/generated/modelcif_input/
```

## Synthetic-coordinate continuation

Real collaborator coordinate assets were not available for this development
cycle. The companion runner therefore creates deterministic, internally
consistent software fixtures from two explicitly named local donor pairs and
runs all seven cases through final ModelCIF:

```bash
.venv/bin/python scripts/generate_fragment_coordinate_e2e_example.py
```

The default output is outside the repository at:

```text
/mnt/disks/toolkit-data/viruses/fragment_metadata_synthetic_e2e/generated/
```

The runner safely recreates only the final `generated/` directory. Its
`source_inputs/` contains an external snapshot of the four small source
fixtures, while `input/` contains seven PDB/meta pairs. Later directories
contain converted confidence/PAE JSONs, merged manifests, individual and
batched metadata, ModelCIF inputs, and seven final ModelCIF files. Exact
commands are in `logs/commands.log`, and `run_summary.json` is written only
after all checks pass.

This focused run deliberately excludes iPSAE. Its complex model and chain
metadata batches are therefore pre-iPSAE integration fixtures, not final
production-schema-ready deliverables: the production schemas require the
iPSAE metric bundle for complexes. No placeholder metrics are invented. A
production run must calculate and enrich those metrics before final schema
validation and release.

Coordinates and scores are deterministic synthetic test data. Residues are
normalized to alanine, repeated donor segments are translated, and pLDDT/PAE
values are remapped consistently. The geometry and scores are unsuitable for
science. Source paths, SHA256 checksums, target lengths/ranges, and crop/tile
details are recorded in `reports/synthetic_provenance.json`.

The established full E2E starts from raw `<AF-ID>-model_v1.pdb` and raw
`<AF-ID>-meta_v1.json`; the latter must contain pLDDT, PAE, and max-PAE data
(and may include ipTM). It does not take mmCIF as input: mmCIF/ModelCIF is
generated downstream.

## GCP handoff

To extend the seven real cases through the full coordinate workflow, provide:

- the bucket URI or prefix holding the files;
- the precise object naming and directory layout for each AF-ID;
- the required authentication method and preferred tool: `gcloud storage` or
  `gsutil`;
- whether the score files are raw `*-meta_v1.json` files or a different
  confidence/PAE layout requiring conversion.

All seven AF-IDs can be replaced with authentic assets once that information
is available. The synthetic-coordinate runner remains useful as a deterministic
regression test, but its outputs must not be treated as biological models.

## Production integration

With a real UniProt DuckDB, use the same preparation order:

```bash
add-fragment-metadata --db uniprot.duckdb --fragments fragments.json \
  --report fragment_load_report.json
enrich-fragment-manifest --input collaborator_manifest.csv \
  --output canonical_chain_manifest.csv --db uniprot.duckdb \
  --report fragment_manifest_report.json
```

Use `--strict` on enrichment when every requested fragment must resolve; omit
it to retain unmatched rows with an empty `protein_name` and review the report.
The loader records a run ID and supports a conflict-safe
`add-fragment-metadata --restore RUN_ID` if the annotations need to be undone.
Pass the enriched canonical manifest to the existing converter and metadata and
ModelCIF exporters. The mock database in this example is only a fixture; do
not use it instead of the production UniProt cache.

## Design lessons captured here

- Fragment data lives in a separate DuckDB child table, leaving `entry`
  unchanged.
- A fragment is identified by UniProt accession and a 1-based, inclusive
  residue range.
- The collaborator-wide adapter used here maps its A/B columns to the
  canonical one-row-per-chain representation.
- Entity IDs represent components: identical full proteins or identical
  fragments share an entity; different ranges are distinct entities.
- Homo/heteromer classification uses the same component identity. Different
  ranges of one accession are heteromeric; repeated copies of one range are
  homomeric.
- A non-empty manifest `protein_name` wins; otherwise consumers retain the
  established UniProt-name fallback.
- No Nextflow dependency is involved, and manifests without `protein_name`
  retain their previous downstream behaviour.
