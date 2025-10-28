# Examples Overview

This directory collects a small, end-to-end ColabFold run that we use when
manually exercising the metadata and conversion tooling. Three monomeric UniProt
entries (Q9TVL3, O00400, and O64637) were picked at random from UniProt release
2025_03 and predicted with the default ColabFold/AlphaFold2 monomer settings.
After the run finished we renamed the top-ranked model for each target to an
AlphaFold-style accession (`AF-000000000000000{1,2,3}`) so that the surrounding
files resemble the public AlphaFold DB layout.

## Directory layout

- **Root files** – treated as the canonical artifacts you might receive from a
  collaborator: per‑target structure files, quality metrics, alignment archives,
  and supporting provider metadata.
- **`tmp/`** – working files captured straight from ColabFold while we were
  iterating on the pipeline. They are useful for reproducing the workflow but
  can be discarded once the ModelCIF conversion succeeds.

## File guide

| Location / pattern | Description |
| --- | --- |
| `AF-000000000000000*-model_v1.cif` / `*.bcif` | Top-ranked models emitted by ColabFold, saved in text and binary mmCIF formats for the three UniProt accessions. |
| `AF-000000000000000*-model_v1.pdb` | PDB files generated from those mmCIF models via our post-processing. |
| `AF-000000000000000*-confidence_v1.json` | Per-residue pLDDT scores exported by ColabFold. |
| `AF-000000000000000*-predicted_aligned_error_v1.json` | Predicted aligned error (PAE) matrices from the same run. |
| `AF-000000000000000*-msa*.a3m` | The multiple-sequence alignments ColabFold constructed before prediction. |
| `AF-metadata-1-of-1.json` | Aggregate metadata record covering all three models (used by higher-level ingestion tests). |
| `DB1.json` | Example provider manifest referenced by the metadata JSON. |
| `sequences.fasta` | The FASTA batch submitted to ColabFold. |
| `tmp/AF-000000000000000*-metadata_for_model_gen.json` | Hand-written inputs that satisfy our ModelCIF schema and drive `run-modelcif-gen`. |
| `tmp/AF-000000000000000*.json` | Raw ColabFold JSON output (per-residue pLDDT array, PAE matrix, global `ptm`, and `max_pae`). |
| `tmp/AF-000000000000000*-model_v1.pdb` | The top-ranked PDB files straight from ColabFold, renamed to the AF accession. |

## Usage notes

1. To regenerate the dataset, run ColabFold on `sequences.fasta` with default
   monomer settings, copy the top-ranked outputs, and rename them to the desired
   AF accession numbers.
2. Update or regenerate the `*_metadata_for_model_gen.json` files if your schema
   changes. They were authored manually to exercise the validator.
3. Once you produce new ModelCIF files with `run-modelcif-gen`, the artifacts at
   the directory root should be all you need to share. The `tmp/` directory is
   purely auxiliary.
