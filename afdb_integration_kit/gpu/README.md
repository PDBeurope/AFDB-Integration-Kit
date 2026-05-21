# Clashes - GPU-Accelerated Protein Structure Analysis

High-throughput analysis of protein structures for clashes and interfaces using GPU acceleration. Designed for processing millions of generated complexes.

## Features

- **Fast**: ~1500 proteins/s on a single GPU
- **Scalable**: Process 10M complexes in ~2 hours
- **Batched**: Efficient GPU memory utilization
- **Parallel parsing**: Multi-CPU PDB parsing with fastpdb
- **Complete analysis**: Clashes (backbone + heavy) and interface residues

## Installation

```bash
# Project production dependencies
uv pip install '.[production]'

# Check the installed PyTorch build
python -c "import torch; print(torch.__version__, torch.version.cuda)"

# torch_cluster must match the installed PyTorch version and CUDA runtime.
# Replace the URL suffix with the wheel index for your environment.
uv pip install torch_cluster -f https://data.pyg.org/whl/torch-2.8.0+cu128.html
```

Use `+cpu` for CPU-only PyTorch builds, or a CUDA suffix such as `+cu118`,
`+cu121`, `+cu124`, `+cu126`, or `+cu128` when it matches your installed
PyTorch wheel. Available wheels are listed at https://data.pyg.org/whl/.

## Quick Start

```python
from afdb_integration_kit.gpu import analyze_pdb_files_pipelined

# Analyze PDB files and write results to JSON
results = analyze_pdb_files_pipelined(
    ["complex1.pdb", "complex2.pdb", ...],
    output_dir="results",
    batch_size=512,
    device="cuda",
    n_workers=16,
)

# Each result contains:
for r in results:
    print(f"{r.path}:")
    print(f"  Backbone clashes: {r.n_backbone_clashes}")
    print(f"  Heavy atom clashes: {r.n_heavy_clashes}")
    print(f"  Clashing residues: {r.backbone_clashing_residues}")
    print(f"  Interface residues: {r.interface_residues}")
```

## Module Structure

```
gpu/
├── __init__.py      # Package exports
├── protein.py       # Protein dataclass
├── parse.py         # PDB parser (fastpdb + multiprocessing)
├── batch.py         # GPU batching utilities
├── clashes.py       # Clash detection
├── interface.py     # Interface residue detection
├── analyze.py       # Master analysis pipeline
└── README.md
```

## Output Format

```json
[
  {
    "path": "complex1.pdb",
    "n_residues": 512,
    "n_atoms": 4190,
    "n_backbone_clashes": 97,
    "n_heavy_clashes": 561,
    "backbone_clashing_residues": [
      {"res_id": 45, "chain_id": "A"},
      {"res_id": 52, "chain_id": "A"}
    ],
    "heavy_clashing_residues": [
      {"res_id": 45, "chain_id": "A"},
      {"res_id": 52, "chain_id": "A"},
      {"res_id": 78, "chain_id": "B"}
    ],
    "interface_residues": [
      {"res_id": 100, "chain_id": "A"},
      {"res_id": 200, "chain_id": "B"}
    ]
  }
]
```

## Definitions

### Clashes
Two atoms clash when their distance is less than the sum of their VDW radii minus a tolerance:
```
distance < (vdw_i + vdw_j) * (1 - clash_cutoff)
```

VDW radii (element-specific):
- C: 1.70 Å, N: 1.55 Å, O: 1.52 Å, S: 1.80 Å

Default parameters:
- `clash_cutoff=0.4` (40% overlap)
- `min_seq_sep=3` (ignore i, i+1, i+2 neighbors)

### Interface Residues
A residue is an interface residue if its CA atom is within `interface_cutoff` (default 8.0 Å) of a CA atom from a different chain.
