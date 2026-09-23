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
# PyTorch (adjust for your CUDA version)
pip install torch --index-url https://download.pytorch.org/whl/cu118

# torch_cluster
pip install torch_cluster -f https://data.pyg.org/whl/torch-2.0.0+cu118.html

# Other dependencies
pip install fastpdb biotite numpy tqdm
```

## Quick Start

```python
from gpu import analyze_pdb_files_pipelined

# Analyze PDB files and write results to JSON
results = analyze_pdb_files_pipelined(
    ["complex1.pdb", "complex2.pdb", ...],
    output_path="results",
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

