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
  scripts/
    extract_subset.py         # Stream UniProt releases into entry.parquet
    build_duckdb.py           # Materialise entry.parquet inside DuckDB
    export_metadata.py        # Emit a single metadata JSON entry
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


3. **Export per-accession metadata records**

   ```bash
   python3 uniprot/scripts/export_metadata.py --accession O00400 --db uniprot/outputs/db/uniprot_2025_03.duckdb --config /path/to/dataset_config.json --mapping /path/to/uniprot_afid_mapping.csv --out /path/to/per_accession/O00400.json
   ```


   This command joins UniProt data from DuckDB with the AF-ID mapping and dataset
   defaults to emit one JSON entry.  Run it once per accession (or in parallel
   via Nextflow) to populate `per_accession/` with the full set of metadata
   snippets.


4. **Combine per-accession JSON files**

   ```bash
   python3 uniprot/scripts/combine_metadata.py --input-dir /path/to/per_accession --output-dir /path/to/batches --output-prefix AF-metadata
   ```

   After all single-entry JSONs exist, combine them into batches ready for
   publication.  The script defaults to 10,000 records per file—adjust with
   `--chunk-size` if you need a different batch size.  Each output is named
   `<prefix>-N-of-M.json`.


5. **Run the Nextflow workflow (optional)**

   ```bash
   nextflow run workflow/af_metadata.nf --db uniprot/outputs/db/uniprot_2025_03.duckdb --config /path/to/dataset_config.json --mapping /path/to/uniprot_afid_mapping.csv --per_outdir /path/to/per_accession --batch_outdir /path/to/batches -w "$PWD/work"
   ```

   Nextflow requires Java (e.g. `sudo apt-get install -y openjdk-17-jre`) and a
   local Nextflow installation (`curl -s https://get.nextflow.io | bash`).

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
