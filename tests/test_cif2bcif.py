from __future__ import annotations

import gzip
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

from afdb_integration_kit.cif2bcif import convert


MINIMAL_CIF = """data_test
loop_
_test.id
_test.value
1 .
2 ?
3 7
"""


@pytest.fixture
def sample_cif(tmp_path: Path) -> Path:
    path = tmp_path / "input.cif"
    path.write_text(MINIMAL_CIF)
    return path


def test_run_cif2bcif_defaults_to_molstar_backend(
    monkeypatch: pytest.MonkeyPatch,
    sample_cif: Path,
    tmp_path: Path,
) -> None:
    output = tmp_path / "output.bcif"
    calls: list[str] = []

    def fake_molstar(input_file: Path, output_file: Path) -> bool:
        calls.append("molstar")
        output_file.write_bytes(b"molstar")
        return True

    def fake_biotite(*args, **kwargs) -> bool:
        calls.append("biotite")
        return False

    monkeypatch.setattr(convert, "_run_molstar_cif2bcif", fake_molstar)
    monkeypatch.setattr(convert, "_run_biotite_cif2bcif", fake_biotite)

    assert convert.run_cif2bcif(sample_cif, output) is True
    assert output.read_bytes() == b"molstar"
    assert calls == ["molstar"]


def test_run_molstar_cif2bcif_falls_back_to_npx(
    monkeypatch: pytest.MonkeyPatch,
    sample_cif: Path,
    tmp_path: Path,
) -> None:
    output = tmp_path / "output.bcif"
    commands: list[list[str]] = []

    monkeypatch.setattr(convert, "_molstar_cmd_path", lambda: ["cif2bcif"])
    monkeypatch.setattr(
        convert,
        "_molstar_cmd_npx",
        lambda: ["npx", "--yes", "-p", "molstar", "cif2bcif"],
    )

    def fake_run(cmd, capture_output, text, timeout):
        commands.append(cmd)
        if cmd[0] == "cif2bcif":
            return subprocess.CompletedProcess(cmd, 1, "", "path failure")
        output.write_bytes(b"npx")
        return subprocess.CompletedProcess(cmd, 0, "ok", "")

    monkeypatch.setattr(convert.subprocess, "run", fake_run)

    assert convert._run_molstar_cif2bcif(sample_cif, output) is True
    assert commands == [
        ["cif2bcif", str(sample_cif), str(output)],
        ["npx", "--yes", "-p", "molstar", "cif2bcif", str(sample_cif), str(output)],
    ]


def test_run_cif2bcif_auto_falls_back_to_biotite(
    monkeypatch: pytest.MonkeyPatch,
    sample_cif: Path,
    tmp_path: Path,
) -> None:
    output = tmp_path / "output.bcif"
    calls: list[str] = []

    def fake_molstar(input_file: Path, output_file: Path) -> bool:
        calls.append("molstar")
        return False

    def fake_biotite(input_file: Path, output_file: Path, tmpdir=None) -> bool:
        calls.append("biotite")
        output_file.write_bytes(b"biotite")
        return True

    monkeypatch.setattr(convert, "_run_molstar_cif2bcif", fake_molstar)
    monkeypatch.setattr(convert, "_run_biotite_cif2bcif", fake_biotite)
    monkeypatch.setattr(convert, "_biotite_version_ok", lambda: True)

    assert convert.run_cif2bcif(sample_cif, output, backend="auto") is True
    assert output.read_bytes() == b"biotite"
    assert calls == ["molstar", "biotite"]


def test_run_cif2bcif_rejects_invalid_backend(
    sample_cif: Path,
    tmp_path: Path,
) -> None:
    with pytest.raises(ValueError, match="Invalid backend"):
        convert.run_cif2bcif(sample_cif, tmp_path / "output.bcif", backend="bad")


def test_import_does_not_require_biotite() -> None:
    script = """
import sys
for name in list(sys.modules):
    if name == "biotite" or name.startswith("biotite."):
        sys.modules.pop(name)
import afdb_integration_kit.cif2bcif.convert
print(any(name == "biotite" or name.startswith("biotite.") for name in sys.modules))
"""
    result = subprocess.run(
        [sys.executable, "-c", script],
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == "False"


@pytest.mark.skipif(
    not convert._biotite_version_ok(),
    reason="Biotite is not installed in this test environment",
)
def test_biotite_backend_writes_bcif_and_preserves_missing_masks(
    sample_cif: Path,
    tmp_path: Path,
) -> None:
    from biotite.structure.io.pdbx import BinaryCIFFile

    output = tmp_path / "output.bcif"

    assert convert.run_cif2bcif(sample_cif, output, backend="biotite") is True
    assert output.exists()

    reread = BinaryCIFFile.read(str(output))
    column = reread["test"]["test"]["value"]
    assert column.as_array(str).tolist() == [".", "?", "7"]
    assert column.mask.array.tolist() == [1, 2, 0]


@pytest.mark.skipif(
    not convert._biotite_version_ok(),
    reason="Biotite is not installed in this test environment",
)
def test_biotite_backend_writes_bcif_gz(sample_cif: Path, tmp_path: Path) -> None:
    output = tmp_path / "output.bcif.gz"

    assert convert.run_cif2bcif(sample_cif, output, backend="biotite") is True
    assert output.exists()

    with gzip.open(output, "rb") as handle:
        assert handle.read()


def test_reserve_temp_output_path_is_unique(tmp_path: Path) -> None:
    first = convert._reserve_temp_output_path(str(tmp_path), Path("same.bcif"), ".tmp")
    second = convert._reserve_temp_output_path(str(tmp_path), Path("same.bcif"), ".tmp")

    try:
        assert first != second
        assert first.exists()
        assert second.exists()
    finally:
        first.unlink(missing_ok=True)
        second.unlink(missing_ok=True)


@pytest.mark.skipif(
    not convert._biotite_version_ok(),
    reason="Biotite is not installed in this test environment",
)
def test_biotite_backend_handles_cross_device_rename_fallback(
    monkeypatch: pytest.MonkeyPatch,
    sample_cif: Path,
    tmp_path: Path,
) -> None:
    output = tmp_path / "output.bcif"
    moved: list[tuple[str, str]] = []
    original_replace = Path.replace

    def fake_replace(self: Path, target: Path) -> Path:
        raise OSError("cross-device link")

    def fake_move(src: str, dst: str) -> str:
        moved.append((src, dst))
        shutil.copyfile(src, dst)
        Path(src).unlink()
        return dst

    monkeypatch.setattr(Path, "replace", fake_replace)
    monkeypatch.setattr(convert.shutil, "move", fake_move)

    try:
        assert convert._run_biotite_cif2bcif(sample_cif, output) is True
    finally:
        monkeypatch.setattr(Path, "replace", original_replace)

    assert output.exists()
    assert moved
