# Post-processing Pipeline

Run on outputs through the full AFDB post-processing pipeline
(metadata, ModelCIF, DSSP, ipSAE, clash/interface analysis, BinaryCIF).

## Prerequisites

```bash
# Set up the environment (from the repo root)
conda env create -f environment.yml
conda activate afdb-toolkit
npm install -g molstar
```

The ipSAE C++ binary is compiled automatically on first run.
To build it manually: `cd ../afdb_integration_kit/ipsae && make`

## Two modes: homodimers inputs vs heterodimers files

Stages 1-16 are identical regardless of model type. The difference is
whether you need **Stage 0** to prepare your inputs first.

![Pipeline flow](../assets/pipeline_flow.png)

- **Prepared inputs** — you already have `{id}-model_v1.pdb` +
  `{id}-meta_v1.json` pairs and all config files (mapping, manifest,
  DuckDB, etc.). The pipeline runs Stages 1-16 directly.

- **Raw OF2 outputs** — you have `*_scores_*.json` +
  `*_unrelaxed_*.pdb` straight from OF2. Pass `--heterodimers`
  to enable **Stage 0**, which does two things:

  1. **File formatting** — scans for PDB+JSON pairs and symlinks them
     into the canonical AFDB layout (`{id}-model_v1.pdb` +
     `{id}-meta_v1.json`). This part is generic.

  2. **Chain-to-accession resolution** — parses model IDs to figure out
     the chain mapping. This IS model-type-aware: heterodimer IDs like
     `AF_XXXX_AF_YYYY` are split into two AF IDs (chain A + chain B
     get different accessions), while homodimer IDs like `AF_XXXX`
     replicate the same accession for both chains. It then resolves
     each AF ID to a real UniProt accession and builds all config files.

---

## Prepared inputs (no Stage 0)

This is the typical path for **homodimers** in the existing AFDB pipeline,
where inputs are already in the canonical layout and all config files
exist. No file renaming, no chain resolution — the pipeline starts
directly at Stage 1.

### What you need

| File | Description |
|------|-------------|
| **Prepared input dir** | Directory with `{id}-model_v1.pdb` + `{id}-meta_v1.json` pairs |
| **Mapping file** | TSV with one model ID per line |
| **Chain mapping CSV** | `model_entity_id,chain_id,uniprot_ac` mapping |
| **Dataset config JSON** | Dataset configuration (provider ID, tool, version) |
| **Provider JSON** | Provider metadata (name, URL, copyrights) |
| **UniProt DuckDB** | UniProt database for your accessions |

### Run

```bash
python production_pipeline.py \
    --output-dir /path/to/output \
    --input-dir /path/to/prepared_inputs \
    --mapping-file /path/to/mapping.tsv \
    --chain-mapping /path/to/manifest.csv \
    --dataset-config /path/to/dataset_config.json \
    --provider-json /path/to/provider.json \
    --uniprot-db /path/to/uniprot.duckdb \
    --cif-qa-metrics auto
```

---

## Raw Heterodimers outputs (Stage 0 prepares everything)

This is the typical path for **heterodimers** (or any new predictions
straight from ColabFold). The raw outputs need renaming and config
files need to be generated. Pass `--heterodimers` to enable Stage 0.

### What you need

| File | Description |
|------|-------------|
| **Input dir** | Directory with matched `*_scores_*.json` + `*_unrelaxed_*.pdb` pairs |
| **Chain mapping CSV** + **UniProt DuckDB** | Pre-built assets |

### Run

Provide a pre-built chain mapping and DuckDB. Config files (mapping TSV,
dataset config, provider JSON) are auto-generated if not provided. Model
IDs are derived from the chain mapping CSV.

```bash
python production_pipeline.py \
    --output-dir /path/to/output \
    --input-dir /path/to/raw_colabfold \
    --heterodimers \
    --chain-mapping /path/to/manifest.csv \
    --uniprot-db /path/to/uniprot.duckdb \
    --cif-qa-metrics auto
```

---

## Common options

```
--workers N                Parallel workers (default: all CPUs)
--resume                   Skip already-completed stages
--dry-run                  Show what would run without executing
--skip-stages S            Comma-separated stages to skip (e.g. stage_12,stage_13)
--dssp-algorithm ALG       psea | pydssp | tmalign (default: pydssp)
--pae-cutoff F             ipSAE PAE threshold (default: 10.0)
--dist-cutoff F            ipSAE distance threshold (default: 15.0)
--clash-cutoff F           VDW overlap threshold (default: 0.4 Å)
--interface-cutoff F       CA-CA interface distance (default: 8.0 Å)
--analysis-batch-size N    Batch size for clash/interface analysis (default: 4)
--cif-qa-metrics M         QA metrics to embed in mmCIF: 'auto' (default) or
                           comma-separated (e.g. 'ipsae_AB,iptm_af,N_clash_backbone')
--enrichment-metrics M...  Metric names for model/chain JSON enrichment (default: all)
--interface-clash-analysis Which analyses: interface backbone_clashes heavy_atom_clashes
--modelcif-template PATH   ModelCIF metadata template JSON
```

## Pipeline stages (execution order)

ipSAE and clash analysis run *before* metadata export so that quality
metrics are available for JSON enrichment and CIF embedding.

```
Stage  1  Prepare assets          (symlink PDB + meta JSON to staging)
Stage  2  Validate assets         (PDB/JSON consistency checks)
Stage  3  Convert ColabFold       (scores → AFDB confidence + PAE JSONs)
Stage  4  Merge manifests         (merge per-model chain/model manifests)
Stage 12  ipSAE                   (interface quality: ipSAE, pDockQ, LIS)
Stage 13  Clash & interface       (VDW clashes, interface residues)
Stage  5  Export model metadata   (per-model JSONs, enriched with iPSAE/clash)
Stage  6  Export chain metadata   (per-chain JSONs, enriched with iPSAE)
Stage  7  Combine model metadata  (batch into chunked JSONs)
Stage  8  Combine chain metadata  (batch into chunked JSONs)
Stage  9  Export ModelCIF input   (prepare ModelCIF metadata from template)
Stage 10  Generate ModelCIF       (PDB → mmCIF with metadata + QA metrics)
Stage 11  DSSP                    (3-state secondary structure annotation)
Stage 14  Enrich PDB              (add AFDB headers to PDB files)
Stage 15  CIF → BCIF              (BinaryCIF conversion)
Stage 16  Cleanup                 (currently no-op, keeps intermediates)
```

## Output directory layout

```
output/
├── staging/                    Stage 1 symlinks
├── scores/                     Stage 3 confidence + PAE JSONs
├── merged_manifests/           Stage 4 merged CSVs
├── ipsae/                      Stage 12 ipSAE summary CSV
├── clash_interface_analysis/   Stage 13 per-model clash/interface JSONs
├── model_jsons/                Stage 5 per-model metadata (enriched)
├── chain_jsons/                Stage 6 per-chain metadata (enriched)
├── model_batches/              Stage 7 batched model metadata
├── chain_batches/              Stage 8 batched chain metadata
├── modelcif_input/             Stage 9 ModelCIF metadata
├── modelcif/                   Stage 10 mmCIF files (with QA metrics)
├── dssp/                       Stage 11 final annotated mmCIF files
├── modelpdb/                   Stage 14 enriched PDB files
├── bcif/                       Stage 15 BinaryCIF files
├── config/                     Auto-generated configs (heterodimer mode)
├── logs/                       Pipeline + error logs
└── pipeline_results.json
```

## Resuming after failure

If a stage fails, fix the issue and re-run with `--resume`:

```bash
python production_pipeline.py \
    --output-dir /path/to/output \
    --input-dir /path/to/raw_colabfold \
    --heterodimers \
    --chain-mapping /path/to/manifest.csv \
    --uniprot-db /path/to/uniprot.duckdb \
    --resume
```

Completed stages are skipped automatically. Delete `.pipeline_cache.json`
in the output directory to force a full re-run.

## Preparing inputs standalone

`prepare_inputs.py` can be used independently to set up the input layout
without running the full pipeline:

```bash
python prepare_inputs.py \
    --input-dir /path/to/raw_colabfold \
    --output-dir /path/to/workdir \
    --chain-mapping /path/to/manifest.csv \
    --uniprot-db /path/to/uniprot.duckdb \
    --provider-id afcdb-heterodimers \
    --provider-name "AFCDB Heterodimers"
```

It prints a ready-to-paste command for the next step (running the pipeline
with the generated config files).
