# Pipeline Optimization Analysis

## Overview

This document analyzes the Python-level optimizations made to the AFDB Integration Kit and their actual impact on the Nextflow workflow (`workflow/end_to_end_with_validation_multibatch.nf`).

---

## Optimizations Implemented

### 1. JSON Parsing: `json` → `orjson`

**Location:** `afdb_integration_kit/validation/validators/pae.py`, `plddt.py`

**Before:**
```python
import json
data = json.loads(path.read_text(encoding="utf-8"))
```

**After:**
```python
import orjson
data = orjson.loads(path.read_bytes())
```

**Benefit:** `orjson` is 3-10x faster than stdlib `json` and uses `read_bytes()` to avoid encoding overhead.

---

### 2. Decimal Validation: `Decimal` → float math

**Location:** `afdb_integration_kit/validation/validators/pae.py`, `plddt.py`

**Before:**
```python
from decimal import Decimal, InvalidOperation

def _has_two_decimal_places(value: object) -> bool:
    try:
        dec = Decimal(str(value))
    except (InvalidOperation, ValueError):
        return False
    exponent = -dec.as_tuple().exponent if dec.as_tuple().exponent < 0 else 0
    return exponent <= 2
```

**After:**
```python
def _has_two_decimal_places(value: object) -> bool:
    try:
        f = float(value)
        return f == round(f, 2)
    except (TypeError, ValueError):
        return False
```

**Benefit:** Avoids expensive arbitrary-precision `Decimal` operations for simple 2-decimal-place checks.

---

### 3. PAE Matrix Validation: Nested Loops → NumPy Vectorized

**Location:** `afdb_integration_kit/validation/validators/pae.py`

**Before:**
```python
for row in matrix:
    if not isinstance(row, list) or len(row) != size:
        # error...
    for idx, value in enumerate(row):
        issue = _validate_value(value, path, f"predicted_aligned_error[{idx}]", enforce_decimal_places)
```

**After:**
```python
import numpy as np

arr = np.array(matrix, dtype=np.float64)  # One-shot conversion + type validation

if arr.ndim != 2 or arr.shape[0] != arr.shape[1]:
    # error...

if enforce_decimal_places:
    rounded = np.round(arr, 2)
    if not np.allclose(arr, rounded, rtol=0, atol=1e-9):
        # error...
```

**Benefit:** O(n²) Python loops replaced with vectorized NumPy operations. For a 1000×1000 PAE matrix, this is ~100x faster.

---

### 4. PDB Parsing: Line-by-Line → gemmi

**Location:** `afdb_integration_kit/colabfold/converter.py`

**Before:**
```python
def _chain_spans_from_pdb(pdb_path, ...):
    residues = OrderedDict()
    for chain_id, resseq, insertion_code in _iterate_pdb_residues(pdb_path):
        # Manual line-by-line parsing
```

**After:**
```python
import gemmi

def _chain_spans_from_pdb_gemmi(pdb_path, ...):
    structure = gemmi.read_structure(str(pdb_path))
    model = structure[0]
    for chain in model:
        for residue in chain:
            # Fast C++ parsing
```

**Benefit:** `gemmi` is a C++ library with Python bindings, providing 10-100x faster PDB parsing.

---

### 5. Parallel File Validation: Sequential → ProcessPoolExecutor

**Location:** `afdb_integration_kit/validation/validators/pae.py`, `plddt.py`, `_parallel.py`

**Before:**
```python
for path in sorted(candidates):
    # Sequential validation
    results.extend(validate_single_file(path))
```

**After:**
```python
from concurrent.futures import ProcessPoolExecutor

num_workers = min(len(candidates), os.cpu_count() or 4)

if num_workers > 1 and len(candidates) >= 10:
    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        batched = list(executor.map(_validate_single_pae, args))
    return [r for batch in batched for r in batch]
```

**Benefit:** Utilizes all CPU cores for I/O and CPU-bound validation when processing many files.

---

### 6. pLDDT Computation: Pure Python (Kept Simple)

**Location:** `afdb_integration_kit/validation/validators/plddt.py`

```python
mean_score = sum(scores) / len(scores)
```

**Rationale:** For small arrays (typical pLDDT lengths), pure Python `sum()/len()` is faster than NumPy due to array conversion overhead.

---

## Optimization Usage in Workflow

| Optimization | Applied To | Used in Workflow? | Reason |
|--------------|------------|-------------------|--------|
| **orjson** | `validators/pae.py`, `validators/plddt.py` | ❌ **NO** | Workflow uses inline Python with `import json` |
| **Decimal→float** | `validators/pae.py`, `validators/plddt.py` | ❌ **NO** | Validation module not called by workflow |
| **NumPy PAE** | `validators/pae.py` | ❌ **NO** | Validation module not called by workflow |
| **gemmi** | `colabfold/converter.py` | ✅ **YES** | `CONVERT_COLABFOLD` calls `afdb-colabfold-convert` |
| **ProcessPoolExecutor** | `validators/*.py` | ❌ **NO** | Nextflow handles parallelism; validation module not used |
| **Pure Python mean** | `validators/plddt.py` | ❌ **NO** | Validation module not called by workflow |

### Key Finding

**Only 1 of 6 optimizations (gemmi in converter.py) is actually utilized by the workflow.**

The workflow's `VALIDATE_ASSETS` process contains its own inline Python that doesn't use the optimized validation module.

---

## Current Bottlenecks

### Critical: Manifest Re-Parsing (O(n²) Behavior)

**Location:** `EXPORT_MODEL_METADATA`, `EXPORT_CHAIN_METADATA` processes

**Problem:**
```groovy
// For each model_id, both scripts:
// 1. Open a new Python process
// 2. Load the ENTIRE merged manifest CSV
// 3. Query DuckDB
// 4. Write one JSON file
```

**Impact:** For N models:
- N × CSV parses for model metadata
- N × CSV parses for chain metadata  
- = **2N full manifest loads**
- = **2N DuckDB connection opens**

For 10,000 models = 20,000 redundant manifest parses.

---

### Critical: Python Process Spawn Overhead

**Problem:** Every model spawns multiple Python processes:
- `EXPORT_MODEL_METADATA` (1 per model)
- `EXPORT_CHAIN_METADATA` (1 per model)
- `EXPORT_MODELCIF_INPUT` (1 per model)

**Impact:** ~200-500ms Python interpreter startup per process × 3 × N models.

---

### High: VALIDATE_ASSETS Uses Unoptimized Code

**Location:** Inline Python in `VALIDATE_ASSETS` process (lines 190-267)

**Problems:**
```python
import json  # Not orjson

def load_json(path: Path):
    with path.open() as fh:
        return json.load(fh)  # Slow stdlib json

# PDB parsing:
for line in fh:
    if line.startswith(("ATOM", "HETATM")):
        # Manual line-by-line parsing, not gemmi
```

---

### High: Scripts in `uniprot/scripts/` Use stdlib json

**Location:** `export_model_metadata.py`, `export_chain_metadata.py`, `combine_metadata.py`

```python
import json

# Output:
json.dump(record, handle, indent=2, ensure_ascii=False)  # ❌ Not orjson
```

---

### Low: Shell Loop Inefficiency in MERGE_MANIFESTS

**Location:** `MERGE_MANIFESTS` process

**Current:**
```bash
for f in ${chain_files.join(' ')}; do
    if [[ $first -eq 1 ]]; then
        head -n1 "$f"
        first=0
    fi
    tail -n +2 "$f"
done
```

**Better:**
```bash
awk 'FNR==1 && NR!=1 {next} {print}' ${chain_files.join(' ')} > output.csv
```

---

## Recommendations

### Priority 1: Batch Metadata Export

Create a single script that processes all models in one invocation:

```python
# batch_export_metadata.py
def export_all_models(model_ids: List[str], manifest_path: Path, db_path: Path):
    # Load manifest ONCE
    manifest = load_manifest(manifest_path)
    
    # Open DuckDB ONCE
    con = duckdb.connect(str(db_path), read_only=True)
    
    for model_id in model_ids:
        # Process each model without reloading
        ...
```

**Expected Impact:** 10,000x reduction in manifest parsing for large batches.

---

### Priority 2: Update VALIDATE_ASSETS Inline Python

Replace stdlib json with orjson and add gemmi for PDB parsing:

```python
import orjson
import gemmi

def load_json(path: Path):
    return orjson.loads(path.read_bytes())

def pdb_residue_count(path: Path):
    structure = gemmi.read_structure(str(path))
    # ...
```

---

### Priority 3: Update `uniprot/scripts/*.py` to Use orjson

```python
import orjson

# Reading:
data = orjson.loads(path.read_bytes())

# Writing:
with path.open("wb") as handle:
    handle.write(orjson.dumps(record, option=orjson.OPT_INDENT_2))
```

---

### Priority 4: Consider Using Validation Module

Instead of inline Python in `VALIDATE_ASSETS`, call the optimized validators:

```python
from afdb_integration_kit.validation.validators import pae, plddt
```

---

## Summary

| Issue | Impact | Fix Complexity | Priority |
|-------|--------|----------------|----------|
| Manifest re-parsing 2N times | O(n²) → O(n) | Medium | Critical |
| Python process spawn overhead | High latency | Medium | Critical |
| VALIDATE_ASSETS uses stdlib | Per-file slowdown | Easy | High |
| Scripts use stdlib json | I/O overhead | Easy | High |
| Shell loops in merge | Minor | Easy | Low |
