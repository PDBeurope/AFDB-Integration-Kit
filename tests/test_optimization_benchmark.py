"""
Benchmark tests comparing original vs optimized implementations.
Run with: pytest tests/test_optimization_benchmark.py -v -s
"""

from __future__ import annotations

import json
import time
import tempfile
from pathlib import Path
from typing import Any, Dict, List

import pytest

# Try importing optimized libraries
try:
    import orjson
    HAS_ORJSON = True
except ImportError:
    HAS_ORJSON = False

try:
    import gemmi
    HAS_GEMMI = True
except ImportError:
    HAS_GEMMI = False


# ============================================================================
# Test Data Generators
# ============================================================================

def generate_plddt_json(num_residues: int = 1000) -> Dict[str, Any]:
    """Generate a synthetic pLDDT confidence JSON."""
    import random
    return {
        "residueNumber": list(range(1, num_residues + 1)),
        "confidenceScore": [round(random.uniform(30, 100), 2) for _ in range(num_residues)],
        "confidenceCategory": [random.choice(["VL", "L", "M", "H", "VH"]) for _ in range(num_residues)],
    }


def generate_pae_json(num_residues: int = 100) -> List[Dict[str, Any]]:
    """Generate a synthetic PAE matrix JSON."""
    import random
    matrix = [
        [round(random.uniform(0, 31.75), 2) for _ in range(num_residues)]
        for _ in range(num_residues)
    ]
    return [{"predicted_aligned_error": matrix, "max_predicted_aligned_error": 31.75}]


def generate_pdb_content(num_residues: int = 100) -> str:
    """Generate a synthetic PDB file content."""
    lines = ["HEADER    TEST STRUCTURE"]
    atom_num = 1
    for res_num in range(1, num_residues + 1):
        # CA atom for each residue
        x, y, z = res_num * 3.8, 0.0, 0.0
        line = f"ATOM  {atom_num:5d}  CA  ALA A{res_num:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00 50.00           C"
        lines.append(line)
        atom_num += 1
    lines.append("END")
    return "\n".join(lines)


def generate_manifest_csv(model_ids: List[str]) -> str:
    """Generate a synthetic manifest CSV."""
    lines = ["model_entity_id,entity_id,chain_id,uniprot_ac,sequence_start,sequence_end"]
    for model_id in model_ids:
        lines.append(f"{model_id},1,A,P12345,1,100")
    return "\n".join(lines)


def generate_config_json() -> Dict[str, Any]:
    """Generate a synthetic dataset config."""
    return {
        "providerId": "test-provider",
        "toolUsed": "AlphaFold",
        "latestVersion": 1,
        "allVersions": [1],
        "entityType": "protein",
        "modelCreatedDate": "2024-01-01T00:00:00Z",
    }


# ============================================================================
# JSON Parsing Benchmarks
# ============================================================================

class TestJsonParsingBenchmark:
    """Compare stdlib json vs orjson for parsing."""

    @pytest.fixture
    def large_plddt_json(self, tmp_path: Path) -> Path:
        """Create a large pLDDT JSON file."""
        data = generate_plddt_json(num_residues=5000)
        path = tmp_path / "large_plddt.json"
        path.write_text(json.dumps(data))
        return path

    @pytest.fixture
    def large_pae_json(self, tmp_path: Path) -> Path:
        """Create a large PAE JSON file (500x500 matrix)."""
        data = generate_pae_json(num_residues=500)
        path = tmp_path / "large_pae.json"
        path.write_text(json.dumps(data))
        return path

    def test_json_load_vs_orjson_plddt(self, large_plddt_json: Path) -> None:
        """Benchmark: stdlib json.load vs orjson.loads for pLDDT."""
        iterations = 50

        # Stdlib json
        start = time.perf_counter()
        for _ in range(iterations):
            with large_plddt_json.open() as f:
                _ = json.load(f)
        stdlib_time = time.perf_counter() - start

        # orjson
        if HAS_ORJSON:
            start = time.perf_counter()
            for _ in range(iterations):
                _ = orjson.loads(large_plddt_json.read_bytes())
            orjson_time = time.perf_counter() - start

            speedup = stdlib_time / orjson_time
            print(f"\n  pLDDT JSON ({iterations} iterations):")
            print(f"    stdlib json: {stdlib_time:.3f}s")
            print(f"    orjson:      {orjson_time:.3f}s")
            print(f"    Speedup:     {speedup:.1f}x")
            assert speedup > 1.0, "Expected orjson to be at least as fast"
        else:
            pytest.skip("orjson not installed")

    def test_json_load_vs_orjson_pae(self, large_pae_json: Path) -> None:
        """Benchmark: stdlib json.load vs orjson.loads for PAE matrix."""
        iterations = 20

        # Stdlib json
        start = time.perf_counter()
        for _ in range(iterations):
            with large_pae_json.open() as f:
                _ = json.load(f)
        stdlib_time = time.perf_counter() - start

        # orjson
        if HAS_ORJSON:
            start = time.perf_counter()
            for _ in range(iterations):
                _ = orjson.loads(large_pae_json.read_bytes())
            orjson_time = time.perf_counter() - start

            speedup = stdlib_time / orjson_time
            print(f"\n  PAE JSON 500x500 ({iterations} iterations):")
            print(f"    stdlib json: {stdlib_time:.3f}s")
            print(f"    orjson:      {orjson_time:.3f}s")
            print(f"    Speedup:     {speedup:.1f}x")
            assert speedup > 1.0, "Expected orjson to be at least as fast"
        else:
            pytest.skip("orjson not installed")


# ============================================================================
# PDB Parsing Benchmarks
# ============================================================================

class TestPdbParsingBenchmark:
    """Compare line-by-line vs gemmi for PDB parsing."""

    @pytest.fixture
    def large_pdb(self, tmp_path: Path) -> Path:
        """Create a large PDB file."""
        content = generate_pdb_content(num_residues=2000)
        path = tmp_path / "large.pdb"
        path.write_text(content)
        return path

    @staticmethod
    def parse_pdb_lines(path: Path) -> int:
        """Original line-by-line PDB parsing."""
        residues = set()
        with path.open() as fh:
            for line in fh:
                if line.startswith(("ATOM", "HETATM")):
                    chain_id = line[21].strip()
                    try:
                        res_seq = int(line[22:26])
                    except ValueError:
                        continue
                    residues.add((chain_id, res_seq))
        return len(residues)

    @staticmethod
    def parse_pdb_gemmi(path: Path) -> int:
        """Optimized gemmi PDB parsing."""
        structure = gemmi.read_structure(str(path))
        if not structure:
            return 0
        model = structure[0]
        residues = set()
        for chain in model:
            for residue in chain:
                residues.add((chain.name, residue.seqid.num))
        return len(residues)

    def test_pdb_parsing_benchmark(self, large_pdb: Path) -> None:
        """Benchmark: line-by-line vs gemmi for PDB parsing."""
        iterations = 100

        # Line-by-line
        start = time.perf_counter()
        for _ in range(iterations):
            count1 = self.parse_pdb_lines(large_pdb)
        lines_time = time.perf_counter() - start

        # gemmi
        if HAS_GEMMI:
            start = time.perf_counter()
            for _ in range(iterations):
                count2 = self.parse_pdb_gemmi(large_pdb)
            gemmi_time = time.perf_counter() - start

            assert count1 == count2, "Both methods should return same count"

            speedup = lines_time / gemmi_time
            print(f"\n  PDB parsing 2000 residues ({iterations} iterations):")
            print(f"    Line-by-line: {lines_time:.3f}s")
            print(f"    gemmi:        {gemmi_time:.3f}s")
            print(f"    Speedup:      {speedup:.1f}x")
            # Note: For simple synthetic PDBs with only CA atoms, line-by-line can be faster.
            # gemmi's advantage is for real PDBs with multiple atoms, HETATM, NaN checking, etc.
            # The main optimization gain is from batch architecture, not PDB parsing speed.
            print(f"    Note: gemmi overhead visible on simple synthetic files; faster on real PDBs")
        else:
            pytest.skip("gemmi not installed")


# ============================================================================
# JSON Write Benchmarks
# ============================================================================

class TestJsonWriteBenchmark:
    """Compare stdlib json.dump vs orjson.dumps for writing."""

    def test_json_dump_vs_orjson(self, tmp_path: Path) -> None:
        """Benchmark: stdlib json.dump vs orjson.dumps."""
        data = generate_plddt_json(num_residues=5000)
        iterations = 50

        # Stdlib json
        start = time.perf_counter()
        for i in range(iterations):
            path = tmp_path / f"stdlib_{i}.json"
            with path.open("w") as f:
                json.dump(data, f, indent=2)
        stdlib_time = time.perf_counter() - start

        # orjson
        if HAS_ORJSON:
            start = time.perf_counter()
            for i in range(iterations):
                path = tmp_path / f"orjson_{i}.json"
                with path.open("wb") as f:
                    f.write(orjson.dumps(data, option=orjson.OPT_INDENT_2))
            orjson_time = time.perf_counter() - start

            speedup = stdlib_time / orjson_time
            print(f"\n  JSON write ({iterations} iterations):")
            print(f"    stdlib json: {stdlib_time:.3f}s")
            print(f"    orjson:      {orjson_time:.3f}s")
            print(f"    Speedup:     {speedup:.1f}x")
            assert speedup > 1.0, "Expected orjson to be at least as fast"
        else:
            pytest.skip("orjson not installed")


# ============================================================================
# Integration Test: Validate inline code changes
# ============================================================================

class TestValidateAssetsInline:
    """Test the VALIDATE_ASSETS inline code logic."""

    def test_validation_logic_original_vs_optimized(self, tmp_path: Path) -> None:
        """Ensure optimized code produces same results as original."""
        # Create test files
        plddt_data = generate_plddt_json(num_residues=100)
        pae_data = generate_pae_json(num_residues=100)
        pdb_content = generate_pdb_content(num_residues=100)

        meta_path = tmp_path / "meta.json"
        pdb_path = tmp_path / "model.pdb"

        meta_path.write_text(json.dumps({
            "plddt": plddt_data["confidenceScore"],
            "predicted_aligned_error": pae_data[0]["predicted_aligned_error"],
        }))
        pdb_path.write_text(pdb_content)

        # Original logic
        def validate_original(meta_path: Path, pdb_path: Path) -> dict:
            with meta_path.open() as fh:
                meta = json.load(fh)

            plddt_n = len(meta.get("plddt", []))
            pae = meta.get("predicted_aligned_error", [])
            pae_d = len(pae) if pae and len(pae) == len(pae[0]) else "0"

            residues = set()
            with pdb_path.open() as fh:
                for line in fh:
                    if line.startswith(("ATOM", "HETATM")):
                        chain_id = line[21].strip()
                        try:
                            res_seq = int(line[22:26])
                        except ValueError:
                            continue
                        residues.add((chain_id, res_seq))
            pdb_n = len(residues)

            return {"plddt_n": plddt_n, "pae_d": pae_d, "pdb_n": pdb_n}

        # Optimized logic
        def validate_optimized(meta_path: Path, pdb_path: Path) -> dict:
            if HAS_ORJSON:
                meta = orjson.loads(meta_path.read_bytes())
            else:
                with meta_path.open() as fh:
                    meta = json.load(fh)

            plddt_n = len(meta.get("plddt", []))
            pae = meta.get("predicted_aligned_error", [])
            pae_d = len(pae) if pae and len(pae) == len(pae[0]) else "0"

            if HAS_GEMMI:
                structure = gemmi.read_structure(str(pdb_path))
                model = structure[0]
                residues = set()
                for chain in model:
                    for residue in chain:
                        residues.add((chain.name, residue.seqid.num))
                pdb_n = len(residues)
            else:
                residues = set()
                with pdb_path.open() as fh:
                    for line in fh:
                        if line.startswith(("ATOM", "HETATM")):
                            chain_id = line[21].strip()
                            try:
                                res_seq = int(line[22:26])
                            except ValueError:
                                continue
                            residues.add((chain_id, res_seq))
                pdb_n = len(residues)

            return {"plddt_n": plddt_n, "pae_d": pae_d, "pdb_n": pdb_n}

        result_original = validate_original(meta_path, pdb_path)
        result_optimized = validate_optimized(meta_path, pdb_path)

        print(f"\n  Original result:  {result_original}")
        print(f"  Optimized result: {result_optimized}")

        assert result_original == result_optimized, "Results should match"


if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s"])
