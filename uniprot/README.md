# UniProt Metadata Toolkit

This module hosts the utilities we use to extract UniProt records, build a
lightweight DuckDB cache, and prepare AlphaFold metadata JSON files.  All tools
stream the official UniProt flat files, so you can work with Swiss-Prot and
TrEMBL releases without unpacking them in advance.  Everything under
`uniprot/outputs` is considered source-data preparation and can be regenerated
at any time.

## Directory layout

```
uniprot/
  README.md                  # This document
  data/
    uniprot_sprot.dat.gz -> /path/to/your/uniprot_sprot.dat.gz
    uniprot_trembl.dat.gz -> /path/to/your/uniprot_trembl.dat.gz
  outputs/
    parquet/                  # Parquet subset extracted from UniProt
    db/                       # DuckDB cache built from the parquet subset
  templates/
    colabfold_example_modelcif_metadata.json    # Starter template for ModelCIF metadata exports
  scripts/
    shard_uniprot.py          # Split UniProt releases into shards for parallel runs
    extract_subset.py         # Stream UniProt releases into entry.parquet
    merge_parquet_shards.py   # Merge shard-level parquet outputs into one file
    build_duckdb.py           # Materialise entry.parquet inside DuckDB
    add_custom_metadata.py    # Apply or restore field-level entry annotations
    add_fragment_metadata.py  # Load or restore named sequence fragments
    enrich_fragment_manifest.py # Stream-normalise and enrich chain manifests
    export_model_metadata.py  # Emit one model-level JSON document (per model)
    export_chain_metadata.py  # Emit per-chain JSON documents (one per chain)
    combine_metadata.py       # Combine metadata JSON entries into batches
  (Programmatic wrappers live under `afdb_integration_kit.uniprot.*`)
```

All files produced by these utilities are ignored by git; users are expected to
generate them locally.

## Prerequisites

| Requirement | Notes |
|-------------|-------|
| Python ≥ 3.12 | Install dependencies with `pip install -r requirements.txt`. |
| tqdm | Installed via the requirements file (progress bar for subset extraction). |
| DuckDB | Pulled in as a Python dependency. |
| Nextflow & Java (optional) | Required only if you want to run the workflow wrapper. |

The `data/` directory contains symlinks to the UniProt releases.  Update them if
your data lives elsewhere.

## Optional: shard releases for parallel extraction

For large TrEMBL runs, pre-shard the release once so you can process shards in
parallel later.  Shards are stable by accession (md5 hash modulo the shard
count) and keep valid UniProt flat-file syntax. Examples below use the 2025_04
release, but older releases work the same as long as you keep the original
Swiss-Prot/TrEMBL `.dat.gz` files.

```bash
python3 uniprot/scripts/shard_uniprot.py \
  -o uniprot/outputs/shards \
  -r 2025_04 \
  --shard-count 8 \
  uniprot/data/uniprot_sprot.dat.gz
```

Outputs land under `uniprot/outputs/shards/<release>/<source>/` with filenames
that include `sprot`/`trembl` for downstream `reviewed` detection.  Repeat the
command for TrEMBL and choose a higher `--shard-count` if needed.  Run this once
per release.

## Typical workflow

1. **Extract the UniProt subset**

   ```bash
   python3 uniprot/scripts/extract_subset.py --mapping /path/to/uniprot_afid_mapping.csv -o uniprot/outputs/parquet -r 2025_04 uniprot/data/uniprot_sprot.dat.gz uniprot/data/uniprot_trembl.dat.gz
   ```

   The command streams both flat files, writes the matching entries into
   `uniprot/outputs/parquet/entry.parquet`, shows a running count of matched
   accessions, and stops as soon as every accession in the mapping is found.

   Shard-aware (run one command per shard, ideally in parallel):

   ```bash
   python3 uniprot/scripts/extract_subset.py \
     --mapping /path/to/uniprot_afid_mapping.csv \
     --shard-count 8 \
     --shard-index 0 \
     -o uniprot/outputs/parquet/shard-00 \
     -r 2025_04 \
     uniprot/outputs/shards/2025_04/sprot/sprot-shard-00.dat.gz
   ```

   Repeat for each shard index and source (sprot/trembl). Concatenate the shard
   Parquet files afterward (e.g., `duckdb` or `pyarrow` can append
   `uniprot/outputs/parquet/shard-*/entry.parquet` into a single dataset).


   CLI entrypoints (after `pip install -e .`):

   - `afdb-uniprot-shard -o uniprot/outputs/shards -r 2025_04 --shard-count 8 uniprot/data/uniprot_sprot.dat.gz`
   - `afdb-uniprot-extract --mapping /path/to/uniprot_afid_mapping.csv --shard-count 8 --shard-index 0 -o uniprot/outputs/parquet/sprot-shard-00 -r 2025_04 uniprot/outputs/shards/2025_04/sprot/sprot-shard-00.dat.gz`
   - `afdb-uniprot-merge -o uniprot/outputs/parquet/entry.parquet uniprot/outputs/parquet/sprot-shard-*/entry.parquet`
   - `afdb-uniprot-build-db --parquet-dir uniprot/outputs/parquet --db uniprot/outputs/db/uniprot_2025_04.duckdb --force`
   - `add-custom-metadata --db uniprot/outputs/db/uniprot_2025_04.duckdb --annotations custom_metadata.json`
   - `add-fragment-metadata --db uniprot/outputs/db/uniprot_2025_04.duckdb --fragments fragments.json`
   - `enrich-fragment-manifest --input collaborator_manifest.csv --output canonical_chain_manifest.csv --db uniprot/outputs/db/uniprot_2025_04.duckdb`


2. **Build the DuckDB cache**

   ```bash
   python3 uniprot/scripts/build_duckdb.py --parquet-dir uniprot/outputs/parquet --db uniprot/outputs/db/uniprot_2025_04.duckdb --force
   ```

### Add custom metadata to an existing cache

`add-custom-metadata` makes controlled, in-place changes to rows in the
existing DuckDB `entry` table. Existing exporters and workflows continue to
read that table without any configuration or interface changes.

The annotations file must be a JSON object keyed by the exact value of
`entry.primary_ac`. Each value must be an object whose keys are actual column
names in the target database's `entry` table:

```json
{
  "O00400": {
    "protein_full_names": "Preferred protein name",
    "gene_names": "SLC33A1",
    "organism_common_names": ["human"]
  },
  "Q46806": {
    "protein_short_names": ["Hydantoinase"]
  }
}
```

Run it after building the database and before exporting metadata:

```bash
add-custom-metadata \
  --db uniprot/outputs/db/uniprot_2025_04.duckdb \
  --annotations custom_metadata.json
```

The command discovers the `entry` schema from the specified database at
runtime. It does not add columns. `primary_ac` is the immutable lookup key and
cannot itself be changed. A supported column may be populated when its current
value is `NULL`, but JSON `null` is not accepted as an annotation value and
does not mean "clear this field."

Input keys use DuckDB column names, not downstream JSON names. For example,
use `protein_full_names`; an exporter-facing key such as
`uniprotDescription` is an unknown column and is skipped. The accepted JSON
types depend on the inspected DuckDB type:

- `VARCHAR`: a JSON string.
- Integer types: a JSON integer; booleans are not accepted as integers.
- `BOOLEAN`: a JSON boolean.
- `VARCHAR[]`: either a JSON list containing only strings or one JSON string.
  A single string is normalised to a one-item list before comparison and
  storage.

Other DuckDB types are reported as unsupported and skipped rather than
coerced. Unknown accessions, unknown or immutable columns, invalid values, and
unsupported types are per-item warnings: valid independent changes in the
same document still apply.

`sequence`, `length`, and `md5` form one validation bundle. If any one is
present for an accession, all three must be present. `sequence` must be a
non-empty uppercase string containing only
`ACDEFGHIKLMNPQRSTVWYBJOUXZ`; `length` must equal its residue count; and `md5`
must match the MD5 hash of the sequence, case-insensitively. A valid hash is
stored in lowercase. If the bundle is invalid, all three fields are skipped
as one unit while unrelated valid fields for that accession remain eligible.
This checks only the internal DuckDB record; it does not inspect coordinate
files, pLDDT arrays, or chain manifests.

Preview an annotation document with no database changes:

```bash
add-custom-metadata \
  --db uniprot/outputs/db/uniprot_2025_04.duckdb \
  --annotations custom_metadata.json \
  --dry-run
```

Every invocation prints a concise human-readable summary and per-item
warnings to stdout. Add `--report` for a machine-readable JSON report:

```bash
add-custom-metadata \
  --db uniprot/outputs/db/uniprot_2025_04.duckdb \
  --annotations custom_metadata.json \
  --report custom_metadata_report.json
```

The report includes the database and input paths, input SHA-256, run ID,
dry-run state, requested/found/updated/unchanged/skipped/conflicted counts,
old and new values for planned or applied changes, and explicit reasons for
every skipped or conflicted item. In the count summary, `requested` and
`found` count accessions for an annotation run; `updated` and `unchanged`
count fields. Additional `requested_fields` and `valid_fields` counts are
included in the JSON report.

For a non-dry run, all field updates and their history records are committed
in one transaction. History is stored in `custom_annotation_history` with the
run ID, application time, source hash, accession, column name, and exact old
and new values as JSON. Save the printed run ID to reverse the run later:

```bash
add-custom-metadata \
  --db uniprot/outputs/db/uniprot_2025_04.duckdb \
  --restore 6f4c4332-c504-4d92-b33c-28af2e4d9b9e
```

Restore is also transactional and supports `--dry-run` and `--report`. Each
field is restored only if its current value still equals the value written by
that run. If a later annotation or another process changed the value, the
restore reports a conflict and leaves that field untouched; other
non-conflicting fields from the run are still restored.

### Add named fragment metadata and enrich a chain manifest

Use this two-step workflow when a model contains a named portion of a UniProt
entry, such as a polyprotein-derived fragment. It is intentionally separate
from `add-custom-metadata`: the existing `entry` table is never altered.

1. After building the DuckDB cache, load authoritative fragment names with
   `add-fragment-metadata`.
2. Before ColabFold conversion, metadata export, or ModelCIF input export,
   run `enrich-fragment-manifest` to create the canonical one-row-per-chain
   manifest used by those consumers.

Do not derive an authoritative protein name from collaborator `chain_a_id` or
`chain_b_id`. They are source/occupancy identifiers only: the normalized
canonical `chain_id` is the literal `A` or `B`. Supply names in the fragment
JSON (or in an already canonical manifest); the enricher only joins names by
UniProt accession and residue range.

#### Fragment JSON and DuckDB contract

`add-fragment-metadata` accepts one JSON object. Each key is an exact UniProt
accession and each value is a list of fragment records with exactly these
fields:

```json
{
  "P27409": [
    {
      "sequence_start": 1,
      "sequence_end": 46,
      "protein_name": "Example fragment 1"
    }
  ]
}
```

Ranges are 1-based and inclusive. `sequence_start` and `sequence_end` must be
JSON integers, positive, and ordered (`sequence_start <= sequence_end`);
`protein_name` must be a non-empty string. The accession must exist in
`entry.primary_ac`, and the inclusive end must not exceed the corresponding
`entry.sequence` length. Invalid records and unknown accessions are reported
without preventing independent valid records in the same input from loading.
Duplicate `(accession, start, end)` records within one document are skipped.

The command creates this separate table lazily on a non-dry run:

```sql
CREATE TABLE fragment_metadata (
    uniprot_ac VARCHAR NOT NULL,
    sequence_start INTEGER NOT NULL,
    sequence_end INTEGER NOT NULL,
    protein_name VARCHAR NOT NULL,
    PRIMARY KEY (uniprot_ac, sequence_start, sequence_end)
);
```

Thus a range is the identity of a fragment. Re-loading the same key with a
different name updates that name; an identical name is unchanged.

Load and preview examples:

```bash
add-fragment-metadata \
  --db uniprot/outputs/db/uniprot_2025_04.duckdb \
  --fragments fragments.json

add-fragment-metadata \
  --db uniprot/outputs/db/uniprot_2025_04.duckdb \
  --fragments fragments.json \
  --dry-run \
  --report fragment_load_report.json
```

Every run prints a summary. `--report` writes JSON containing paths, the input
SHA-256, dry-run state, run ID (for a real load), requested/found/updated/
unchanged/skipped/conflicted counts, planned or applied old/new values, and
per-record skip or conflict reasons. A dry run opens the database read-only:
it writes neither fragment data nor fragment/history tables.

Non-dry loads record old and new names in `fragment_metadata_history`, together
with a run ID, application time, and source SHA-256. Restore a run with:

```bash
add-fragment-metadata \
  --db uniprot/outputs/db/uniprot_2025_04.duckdb \
  --restore 6f4c4332-c504-4d92-b33c-28af2e4d9b9e \
  --report fragment_restore_report.json
```

Restore also supports `--dry-run`. It restores updates or removes rows that
were inserted by the selected run only when the current name still equals the
name written by that run. Later changes are reported as conflicts and left in
place, while non-conflicting records are restored transactionally.

#### Normalize and enrich collaborator or canonical manifests

`enrich-fragment-manifest` reads a comma- or tab-delimited manifest and writes
a canonical CSV/TSV manifest (using the detected delimiter). Input and output
paths must differ. It auto-detects either of these inputs:

- The 19-column collaborator-wide manifest, identified by `afdb_id`,
  `chain_a_id`, and `chain_a_uniprot`. `afdb_id` becomes
  `model_entity_id`; each populated `chain_a_id`/`chain_b_id` becomes one
  output row with the literal canonical `chain_id` `A`/`B`;
  `chain_*_uniprot` becomes `uniprot_ac`;
  `chain_*_is_polyprotein?` becomes `is_fragment`; and
  `chain_*_start`/`chain_*_end` become the inclusive sequence range. Chain B
  is emitted only when `chain_b_id` is populated. For older wide manifests,
  the command also accepts `chain_a`/`chain_b` occupancy columns and the
  no-question-mark `chain_*_is_polyprotein` aliases.
- An existing canonical manifest, identified by `model_entity_id`, `chain_id`,
  and `uniprot_ac`. Its canonical values are normalized and its additional
  columns are preserved.

The output begins with these canonical columns:

```text
model_entity_id,entity_id,chain_id,uniprot_ac,is_fragment,sequence_start,sequence_end,protein_name
```

For canonical input, additional non-canonical columns are retained after these
columns. For collaborator-wide input, target-level extra columns are retained
on each normalized chain row, while source chain-specific columns (including
the collaborator occupancy, accession, range, polyprotein, and length fields)
are consumed during normalization and are not copied to the output.

`is_fragment` accepts common boolean spellings (`true`/`false`, `1`/`0`,
`yes`/`no`, and `t`/`f`). Fragment ranges are required, positive integers and
are 1-based inclusive; ranges for non-fragments may be blank. `protein_name`
is optional in canonical input. A supplied non-empty fragment name is
preserved; an empty or absent one is looked up in `fragment_metadata` using
`(uniprot_ac, sequence_start, sequence_end)`. Non-fragment rows do not receive
a fragment lookup.

Entity IDs are reassigned sequentially within each model. Repeated chains of
the same whole UniProt entry share an entity ID; repeated chains of the same
fragment share one when their `(uniprot_ac, sequence_start, sequence_end)`
matches. Different ranges of the same accession receive different entity IDs.
Keep rows for each model contiguous in the input so that the streaming entity
assignment remains scoped to that model.

For example:

```bash
enrich-fragment-manifest \
  --input collaborator_manifest.csv \
  --output canonical_chain_manifest.csv \
  --db uniprot/outputs/db/uniprot_2025_04.duckdb \
  --report fragment_manifest_report.json
```

The command streams input rows and output rows rather than retaining the
normalized manifest in memory; it only preloads the small fragment-name
lookup. Its report records the detected schema, input/output row counts,
fragment rows, enriched and preserved names, and unmatched or ambiguous
fragment lookups (with bounded samples). By default, unmatched or ambiguous
lookups are non-blocking and leave `protein_name` empty. With `--strict`, any
such lookup returns a non-zero status and leaves the destination unchanged;
the JSON report is still written when requested.

An empty or missing `protein_name` preserves legacy behavior in downstream
consumers: they use their normal UniProt fallback (`protein_full_names`, then
`protein_short_names`, then the entry name, then the accession). When a
non-empty manifest name is present, it takes precedence. Consequently,
non-fragment manifests and older canonical manifests without `protein_name`
remain compatible with the existing pipeline.


3. **Export per-model metadata records (model-level Solr schema)**

   ```bash
   python3 uniprot/scripts/export_model_metadata.py \
     --model-entity-id AF-0000000000000004 \
    --db uniprot/outputs/db/uniprot_2025_04.duckdb \
     --config /path/to/dataset_config.json \
     --mapping /path/to/uniprot_afid_mapping.csv \
     --model-manifest /path/to/uniprot_model_metadata.csv \
     --out /path/to/per_accession/models/AF-0000000000000004.json
   ```

   - `--mapping` (chain manifest): `model_entity_id,entity_id,chain_id,uniprot_ac,protein_name,sequence_start,sequence_end,is_fragment,is_isoform,entity_type,average_plddt,fraction_plddt_very_low,fraction_plddt_low,fraction_plddt_confident,fraction_plddt_very_high`
   - `--model-manifest` (model manifest, optional): `model_entity_id,ipTM,average_plddt,name,isAMdata`
   - Output: one JSON document per model.

4. **Export per-chain metadata records (chain-level Solr schema)**

   ```bash
   python3 uniprot/scripts/export_chain_metadata.py \
     --model-entity-id AF-0000000000000004 \
    --db uniprot/outputs/db/uniprot_2025_04.duckdb \
     --config /path/to/dataset_config.json \
     --mapping /path/to/uniprot_afid_mapping.csv \
     --model-manifest /path/to/uniprot_model_metadata.csv \
     --out /path/to/per_accession/chains/AF-0000000000000004.json
   ```

   - `--mapping` (chain manifest): `model_entity_id,entity_id,chain_id,uniprot_ac,protein_name,sequence_start,sequence_end,is_fragment,is_isoform,entity_type,average_plddt,fraction_plddt_very_low,fraction_plddt_low,fraction_plddt_confident,fraction_plddt_very_high`
   - Output: an array of JSON documents, one per chain for the requested model.

5. **Combine per-accession JSON files**

   ```bash
   python3 uniprot/scripts/combine_metadata.py --input-dir /path/to/per_accession --output-dir /path/to/batches --output-prefix AF-metadata
   ```

   After all single-entry JSONs exist, combine them into batches ready for
   publication.  The script defaults to 10,000 records per file—adjust with
   `--chunk-size` if you need a different batch size.  Each output is named
   `<prefix>-N-of-M.json`.


5. **Produce ModelCIF generator metadata**

   ```bash
   python3 uniprot/scripts/export_modelcif_input.py \
     --model-id AF-0000000000000004 \
     --manifest examples/config/subset_uniprot_afid_mapping.csv \
   --db uniprot/outputs/db/uniprot_2025_04.duckdb \
     --template uniprot/templates/colabfold_example_modelcif_metadata.json \
     --out examples/AF-0000000000000004_model_gen.json
   ```

   The manifest is a CSV with one row per chain:

   ```
   model_entity_id,entity_id,chain_id,uniprot_ac
   AF-0000000000000001,1,A,Q9TVL3
   AF-0000000000000002,1,A,P12345
   AF-0000000000000003,1,A,Q98765
   ```

   Entities that share the same UniProt accession reuse the same `entity_id`.
   The exporter deduplicates `_ma_target_ref_db_details` per entity while still
   emitting a `_ma_target_entity_instance` row for every chain.

   Copy `uniprot/templates/colabfold_example_modelcif_metadata.json` to your workspace if you need
   to customise provider details, data-usage statements, or software fields
   before running large batches.
   That template is intentionally static scaffolding. It defines ModelCIF
   categories and software parameter metadata such as the iPSAE cutoffs, but it
   must not be used to hard-code computed `complexPredictionAccuracy_*` values.
   For complex examples, those computed interface metrics are produced later by
   the iPSAE enrichment stage, written into the enriched model JSON metadata,
   and then copied from that JSON into ModelCIF global QA metrics.


6. **Run the Nextflow workflows**

   Nextflow requires Java (e.g. `sudo apt-get install -y openjdk-17-jre`) and a
   local Nextflow installation (`curl -s https://get.nextflow.io | bash`). We
   provide two entry points:

   _AF metadata batches (per-accession + combine):_

   ```bash
  nextflow run workflow/af_metadata.nf \
    --db uniprot/outputs/db/uniprot_2025_04.duckdb \
    --config /path/to/dataset_config.json \
    --mapping /path/to/uniprot_afid_mapping.csv \
    --model_manifest /path/to/uniprot_model_metadata.csv \
    --per_outdir /path/to/per_accession/models \
    --batch_outdir /path/to/batches \
    -w "$PWD/work"
   ```

   _ModelCIF generator metadata:_

   ```bash
  nextflow run workflow/modelcif_metadata.nf \
    --db uniprot/outputs/db/uniprot_2025_04.duckdb \
     --manifest examples/config/subset_uniprot_afid_mapping.csv \
     --template uniprot/templates/colabfold_example_modelcif_metadata.json \
     --output_dir examples/modelcif_metadata \
     -w "$PWD/work/modelcif"
   ```

   The ModelCIF workflow consumes the same manifest used by
   `export_modelcif_input.py` and writes one JSON file per model into the
   requested output directory.

## Notes
- The subset extractor caches the target accessions in memory only—no UniProt
  contents are loaded wholesale.  Memory usage therefore scales with the number
  of accessions, not the size of the release.
- `dataset_config.json` should live alongside the model artefacts you intend to
  package (e.g., one lives under `examples/` in this repo).  It holds
  dataset-level values such as `toolUsed`, `modelCreatedDate`, version lists, and
  placeholder pLDDT fractions (initially `null` so you can fill them in after
  batching).
- `uniprot_afid_mapping.csv` links model entity IDs to UniProt accessions.
  Extend or replace it with your own mapping before running the pipeline.
- Generated JSON files and DuckDB outputs are ignored by git, so you can rerun
  the pipeline without polluting the repository.

For further context on the AlphaFold integration toolkit, see the repository
root README.
