import csv
import math
import shutil
import subprocess
from pathlib import Path

import pytest


IPSAE_DIR = Path(__file__).resolve().parents[1] / "afdb_integration_kit" / "ipsae"


def _find_eigen_dir() -> str | None:
    candidates = [
        IPSAE_DIR / "deps" / "eigen-3.4.0",
        Path("/usr/include/eigen3"),
        Path("/usr/local/include/eigen3"),
    ]
    for candidate in candidates:
        if (candidate / "Eigen").is_dir():
            return str(candidate)
    return None


def _write_tiny_batch_fixture(batch_dir: Path) -> None:
    pdb_path = batch_dir / "toy-model_v1.pdb"
    json_path = batch_dir / "toy-meta_v1.json"

    pdb_path.write_text(
        "\n".join(
            [
                "ATOM      1  CA  GLY A   1       0.000   0.000   0.000  1.00 90.00           C",
                "ATOM      2  CA  GLY B   1       4.000   0.000   0.000  1.00 80.00           C",
                "TER",
                "END",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    json_path.write_text(
        (
            '{"pae": [[0.0, 2.0], [3.0, 0.0]], '
            '"plddt": [90.0, 80.0], "iptm": 0.7}\n'
        ),
        encoding="utf-8",
    )


@pytest.mark.skipif(shutil.which("make") is None or shutil.which("g++") is None, reason="make/g++ not available")
def test_ipsae_cpp_makefile_batch_summary_smoke(tmp_path: Path) -> None:
    eigen_dir = _find_eigen_dir()
    if eigen_dir is None:
        pytest.skip("Eigen headers are not available locally")

    binary_path = tmp_path / "ipsae_cpp"
    build = subprocess.run(
        [
            "make",
            "dynamic",
            f"EIGEN_DIR={eigen_dir}",
            f"TARGET={binary_path}",
        ],
        cwd=IPSAE_DIR,
        capture_output=True,
        text=True,
    )
    assert build.returncode == 0, build.stderr or build.stdout
    assert binary_path.exists()

    batch_dir = tmp_path / "batch"
    batch_dir.mkdir()
    _write_tiny_batch_fixture(batch_dir)

    summary_path = tmp_path / "ipsae_summary.csv"
    run = subprocess.run(
        [
            str(binary_path),
            "--batch",
            str(batch_dir),
            "10",
            "5",
            "--summary",
            str(summary_path),
            "--workers",
            "1",
            "--quiet",
            "--no-individual",
        ],
        cwd=tmp_path,
        capture_output=True,
        text=True,
    )
    assert run.returncode == 0, run.stderr or run.stdout
    assert summary_path.exists()

    with summary_path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))

    assert len(rows) == 1
    row = rows[0]
    assert row["iptm_af"] == "0.7"
    assert row["n0chn"] == "2"
    assert row["d0chn"] == "1"
    assert math.isclose(float(row["ipsae_AB"]), 0.2, rel_tol=1e-6)
    assert math.isclose(float(row["ipsae_BA"]), 0.1, rel_tol=1e-6)
    assert math.isclose(float(row["LIS_AB"]), 10.0 / 12.0, rel_tol=1e-6)
    assert math.isclose(float(row["LIS_BA"]), 9.0 / 12.0, rel_tol=1e-6)
    assert float(row["pDockQ"]) > 0.0
