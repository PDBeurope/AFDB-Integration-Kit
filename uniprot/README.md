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
    modelcif_metadata.json    # Starter template for ModelCIF metadata exports
  scripts/
    extract_subset.py         # Stream UniProt releases into entry.parquet
    build_duckdb.py           # Materialise entry.parquet inside DuckDB
    export_model_metadata.py  # Emit one model-level JSON document (per model)
    export_chain_metadata.py  # Emit per-chain JSON documents (one per chain)
    combine_metadata.py       # Combine metadata JSON entries into batches
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

## Typical workflow

1. **Extract the UniProt subset**

   ```bash
   python3 uniprot/scripts/extract_subset.py --mapping /path/to/uniprot_afid_mapping.csv -o uniprot/outputs/parquet -r 2025_03 uniprot/data/uniprot_sprot.dat.gz uniprot/data/uniprot_trembl.dat.gz
   ```

   The command streams both flat files, writes the matching entries into
   `uniprot/outputs/parquet/entry.parquet`, shows a running count of matched
   accessions, and stops as soon as every accession in the mapping is found.


2. **Build the DuckDB cache**

   ```bash
   python3 uniprot/scripts/build_duckdb.py --parquet-dir uniprot/outputs/parquet --db uniprot/outputs/db/uniprot_2025_03.duckdb --force
   ```


3. **Export per-model metadata records (model-level Solr schema)**

   ```bash
   python3 uniprot/scripts/export_model_metadata.py \
     --model-entity-id AF-0000000000000004 \
     --db uniprot/outputs/db/uniprot_2025_03.duckdb \
     --config /path/to/dataset_config.json \
     --mapping /path/to/uniprot_afid_mapping.csv \
     --model-manifest /path/to/uniprot_model_metadata.csv \
     --out /path/to/per_accession/models/AF-0000000000000004.json
   ```

   - `--mapping` (chain manifest): `model_entity_id,entity_id,chain_id,uniprot_ac,sequence_start,sequence_end,is_fragment,is_isoform,entity_type,average_plddt,fraction_plddt_very_low,fraction_plddt_low,fraction_plddt_confident,fraction_plddt_very_high`
   - `--model-manifest` (model manifest, optional): `model_entity_id,ipTM,average_plddt,complexName,isAMdata`
   - Output: one JSON document per model.

4. **Export per-chain metadata records (chain-level Solr schema)**

   ```bash
   python3 uniprot/scripts/export_chain_metadata.py \
     --model-entity-id AF-0000000000000004 \
     --db uniprot/outputs/db/uniprot_2025_03.duckdb \
     --config /path/to/dataset_config.json \
     --mapping /path/to/uniprot_afid_mapping.csv \
     --model-manifest /path/to/uniprot_model_metadata.csv \
     --out /path/to/per_accession/chains/AF-0000000000000004.json
   ```

   Output: an array of JSON documents, one per chain for the requested model.

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
     --manifest examples/complexes/config/subset_uniprot_afid_mapping.csv \
     --db uniprot/outputs/db/uniprot_2025_03.duckdb \
     --template uniprot/templates/modelcif_metadata.json \
     --out examples/multimer_examples/homodimer_metadata_for_model_gen.json
   ```

   The manifest is a CSV with one row per chain:

   ```
   model_entity_id,entity_id,chain_id,uniprot_ac
   AF-0000000000000003,1,A,Q9TVL3   # monomer
   AF-0000000000000004,1,A,Q9TVL3   # homodimer
   AF-0000000000000004,1,B,Q9TVL3
   AF-0000000000000005,1,A,P12345   # heterotrimer example
   AF-0000000000000005,2,B,Q98765
   AF-0000000000000005,2,C,Q98765
   ```

   Entities that share the same UniProt accession reuse the same `entity_id`.
   The exporter deduplicates `_ma_target_ref_db_details` per entity while still
   emitting a `_ma_target_entity_instance` row for every chain.

   Copy `uniprot/templates/modelcif_metadata.json` to your workspace if you need
   to customise provider details, data-usage statements, or software fields
   before running large batches.


6. **Run the Nextflow workflows**

   Nextflow requires Java (e.g. `sudo apt-get install -y openjdk-17-jre`) and a
   local Nextflow installation (`curl -s https://get.nextflow.io | bash`). We
   provide two entry points:

   _AF metadata batches (per-accession + combine):_

   ```bash
  nextflow run workflow/af_metadata.nf \
    --db uniprot/outputs/db/uniprot_2025_03.duckdb \
    --config /path/to/dataset_config.json \
    --mapping /path/to/uniprot_afid_mapping.csv \
    --model_manifest /path/to/uniprot_model_metadata.csv \
    --per_outdir /path/to/per_accession/models \
    --batch_outdir /path/to/batches \
    -w "$PWD/work"
   ```

   _ModelCIF generator metadata (complex submissions):_

   ```bash
   nextflow run workflow/modelcif_metadata.nf \
     --db uniprot/outputs/db/uniprot_2025_03.duckdb \
     --manifest examples/complexes/config/subset_uniprot_afid_mapping.csv \
     --template uniprot/templates/modelcif_metadata.json \
     --output_dir examples/complexes/modelcif_metadata \
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
  package (e.g., one lives under `examples/complexes/` in this repo).  It holds
  dataset-level values such as `toolUsed`, `modelCreatedDate`, version lists, and
  placeholder pLDDT fractions (initially `null` so you can fill them in after
  batching).
- `uniprot_afid_mapping.csv` links model entity IDs to UniProt accessions.
  Extend or replace it with your own mapping before running the pipeline.
- Generated JSON files and DuckDB outputs are ignored by git, so you can rerun
  the pipeline without polluting the repository.

For further context on the AlphaFold integration toolkit, see the repository
root README.
