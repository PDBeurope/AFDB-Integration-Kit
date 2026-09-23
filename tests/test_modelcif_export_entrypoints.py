from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest


@pytest.mark.parametrize(
    "script_name",
    ["batch_export_modelcif_input.py", "export_modelcif_input.py"],
)
def test_modelcif_export_scripts_resolve_package_when_run_directly(
    script_name: str,
) -> None:
    repo_root = Path(__file__).resolve().parent.parent
    script = repo_root / "uniprot" / "scripts" / script_name

    result = subprocess.run(
        [sys.executable, "-I", str(script), "--help"],
        cwd=repo_root,
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 0, result.stderr
