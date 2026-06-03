from __future__ import annotations

import importlib.util
import sys
from pathlib import Path


def _load_helper_module():
    repo_root = Path(__file__).resolve().parent.parent
    module_path = repo_root / "scripts" / "generate_colabfold_e2e_example.py"
    spec = importlib.util.spec_from_file_location("generate_colabfold_e2e_example", module_path)
    module = importlib.util.module_from_spec(spec)
    assert spec is not None
    assert spec.loader is not None
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_chain_manifest_rows_use_local_chain_ranges() -> None:
    module = _load_helper_module()

    models = [
        module.ExampleModel(
            category="homodimer",
            example_id="AF-0000000065760001",
            directory=Path("/tmp/homodimer"),
            chain_spans=[
                {
                    "chain_id": "A",
                    "sequence_start": 1,
                    "sequence_end": 256,
                    "residue_count": 256,
                    "uniprot_ac": "Q6GZX4",
                },
                {
                    "chain_id": "B",
                    "sequence_start": 257,
                    "sequence_end": 512,
                    "residue_count": 256,
                    "uniprot_ac": "Q6GZX4",
                },
            ],
        ),
        module.ExampleModel(
            category="heterodimer",
            example_id="AF-0000000300000101",
            directory=Path("/tmp/heterodimer"),
            chain_spans=[
                {
                    "chain_id": "A",
                    "sequence_start": 1,
                    "sequence_end": 30,
                    "residue_count": 30,
                    "uniprot_ac": "A0ABS2QMZ4",
                },
                {
                    "chain_id": "B",
                    "sequence_start": 31,
                    "sequence_end": 63,
                    "residue_count": 33,
                    "uniprot_ac": "A0ABS2QMF5",
                },
            ],
        ),
    ]

    assert module.chain_manifest_rows(models) == [
        {
            "model_entity_id": "AF-0000000065760001",
            "entity_id": "1",
            "chain_id": "A",
            "uniprot_ac": "Q6GZX4",
            "sequence_start": 1,
            "sequence_end": 256,
        },
        {
            "model_entity_id": "AF-0000000065760001",
            "entity_id": "1",
            "chain_id": "B",
            "uniprot_ac": "Q6GZX4",
            "sequence_start": 1,
            "sequence_end": 256,
        },
        {
            "model_entity_id": "AF-0000000300000101",
            "entity_id": "1",
            "chain_id": "A",
            "uniprot_ac": "A0ABS2QMZ4",
            "sequence_start": 1,
            "sequence_end": 30,
        },
        {
            "model_entity_id": "AF-0000000300000101",
            "entity_id": "2",
            "chain_id": "B",
            "uniprot_ac": "A0ABS2QMF5",
            "sequence_start": 1,
            "sequence_end": 33,
        },
    ]
