import importlib.util
import logging
import subprocess
import sys
from pathlib import Path
from unittest.mock import Mock

import pytest


SCRIPT_PATH = Path(__file__).parents[1] / "uniprot" / "scripts" / "batch_ipsae.py"
SPEC = importlib.util.spec_from_file_location("batch_ipsae", SCRIPT_PATH)
assert SPEC is not None and SPEC.loader is not None
batch_ipsae = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(batch_ipsae)


def test_run_ipsae_batch_invokes_cpp_batch_interface(monkeypatch, tmp_path):
    input_dir = tmp_path / "input"
    summary_csv = tmp_path / "ipsae_summary.csv"
    completed = subprocess.CompletedProcess([], 0, stdout="", stderr="")
    run = Mock(return_value=completed)

    monkeypatch.setattr(batch_ipsae, "IPSAE_BINARY", Path("/opt/ipsae_cpp"))
    monkeypatch.setattr(batch_ipsae, "ensure_ipsae_binary", lambda: True)
    monkeypatch.setattr(batch_ipsae.subprocess, "run", run)

    result = batch_ipsae.run_ipsae_batch(input_dir, summary_csv, 10.0, 8.0, 4)

    assert result is completed
    assert run.call_args.args[0] == [
        "/opt/ipsae_cpp",
        "--batch",
        str(input_dir),
        "10.0",
        "8.0",
        "--summary",
        str(summary_csv),
        "--workers",
        "4",
        "--quiet",
    ]


def test_main_stages_inputs_once_and_counts_summary_rows(
    monkeypatch, tmp_path, caplog
):
    pae_dir = tmp_path / "pae"
    pdb_dir = tmp_path / "pdb"
    output_dir = tmp_path / "output"
    pae_dir.mkdir()
    pdb_dir.mkdir()
    pae_file = pae_dir / "AF-TEST-meta_v1.json"
    pdb_file = pdb_dir / "AF-TEST-model_v1.pdb"
    pae_file.write_text("{}")
    pdb_file.write_text("END\n")
    calls = []

    def fake_run(input_dir, summary_csv, pae_cutoff, dist_cutoff, workers):
        calls.append((input_dir, summary_csv, pae_cutoff, dist_cutoff, workers))
        summary_csv.write_text("model,score\nAF-TEST,0.9\n")
        return subprocess.CompletedProcess([], 0, stdout="", stderr="")

    monkeypatch.setattr(batch_ipsae, "run_ipsae_batch", fake_run)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "batch_ipsae.py",
            "--pae-dir",
            str(pae_dir),
            "--pdb-dir",
            str(pdb_dir),
            "--output-dir",
            str(output_dir),
            "--pae-cutoff",
            "10",
            "--dist-cutoff",
            "8",
            "--workers",
            "4",
        ],
    )

    with caplog.at_level(logging.INFO):
        batch_ipsae.main()

    input_dir = output_dir / "input"
    assert (input_dir / pae_file.name).resolve() == pae_file.resolve()
    assert (input_dir / pdb_file.name).resolve() == pdb_file.resolve()
    assert calls == [
        (input_dir, output_dir / "ipsae_summary.csv", 10.0, 8.0, 4)
    ]
    assert "Summary rows: 1" in caplog.text


def test_main_propagates_ipsae_failure_code(monkeypatch, tmp_path):
    pae_dir = tmp_path / "pae"
    pdb_dir = tmp_path / "pdb"
    output_dir = tmp_path / "output"
    pae_dir.mkdir()
    pdb_dir.mkdir()
    (pae_dir / "AF-TEST-meta_v1.json").write_text("{}")
    (pdb_dir / "AF-TEST-model_v1.pdb").write_text("END\n")

    monkeypatch.setattr(
        batch_ipsae,
        "run_ipsae_batch",
        lambda *args: subprocess.CompletedProcess(
            [], 23, stdout="", stderr="batch failed"
        ),
    )
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "batch_ipsae.py",
            "--pae-dir",
            str(pae_dir),
            "--pdb-dir",
            str(pdb_dir),
            "--output-dir",
            str(output_dir),
        ],
    )

    with pytest.raises(SystemExit) as error:
        batch_ipsae.main()

    assert error.value.code == 23


def test_main_reuse_removes_only_stale_managed_staging_files(
    monkeypatch, tmp_path
):
    pae_dir = tmp_path / "pae"
    pdb_dir = tmp_path / "pdb"
    output_dir = tmp_path / "output"
    input_dir = output_dir / "input"
    pae_dir.mkdir()
    pdb_dir.mkdir()
    input_dir.mkdir(parents=True)
    (pae_dir / "AF-CURRENT-meta_v1.json").write_text("{}")
    (pdb_dir / "AF-CURRENT-model_v1.pdb").write_text("END\n")
    stale_target = tmp_path / "stale"
    stale_target.write_text("stale")
    (input_dir / "AF-STALE-meta_v1.json").symlink_to(stale_target)
    (input_dir / "AF-STALE-model_v1.pdb").write_text("stale")
    unrelated = input_dir / "notes.txt"
    unrelated.write_text("keep")
    managed_directory = input_dir / "AF-DIRECTORY-meta_v1.json"
    managed_directory.mkdir()

    seen_entries = []

    def fake_run(batch_dir, summary_csv, *args):
        seen_entries.extend(sorted(path.name for path in batch_dir.iterdir()))
        summary_csv.write_text("model,score\nAF-CURRENT,0.9\n")
        return subprocess.CompletedProcess([], 0, stdout="", stderr="")

    monkeypatch.setattr(batch_ipsae, "run_ipsae_batch", fake_run)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "batch_ipsae.py",
            "--pae-dir",
            str(pae_dir),
            "--pdb-dir",
            str(pdb_dir),
            "--output-dir",
            str(output_dir),
        ],
    )

    batch_ipsae.main()

    assert "AF-STALE-meta_v1.json" not in seen_entries
    assert "AF-STALE-model_v1.pdb" not in seen_entries
    assert unrelated.read_text() == "keep"
    assert managed_directory.is_dir()


@pytest.mark.parametrize("unsafe_kind", ["traversal", "absolute"])
def test_main_rejects_unsafe_model_ids_without_touching_outside_file(
    monkeypatch, tmp_path, unsafe_kind
):
    base_dir = tmp_path / "base"
    pae_dir = base_dir / "pae"
    pdb_dir = base_dir / "pdb"
    output_dir = base_dir / "output"
    pae_dir.mkdir(parents=True)
    pdb_dir.mkdir()
    (output_dir / "input").mkdir(parents=True)

    if unsafe_kind == "traversal":
        model_id = "../AF-EVIL"
        pae_file = base_dir / "AF-EVIL-meta_v1.json"
        pdb_file = base_dir / "AF-EVIL-model_v1.pdb"
        outside_guard = output_dir / "AF-EVIL-meta_v1.json"
    else:
        model_prefix = base_dir / "AF-EVIL"
        model_id = str(model_prefix)
        pae_file = Path(f"{model_prefix}-meta_v1.json")
        pdb_file = Path(f"{model_prefix}-model_v1.pdb")
        outside_guard = pae_file

    pae_file.write_text("{}")
    pdb_file.write_text("END\n")
    if outside_guard != pae_file:
        outside_guard.write_text("do not touch")
    original_guard = outside_guard.read_bytes()
    model_ids = tmp_path / "model_ids.txt"
    model_ids.write_text(f"{model_id}\n")
    run = Mock()
    monkeypatch.setattr(batch_ipsae, "run_ipsae_batch", run)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "batch_ipsae.py",
            "--pae-dir",
            str(pae_dir),
            "--pdb-dir",
            str(pdb_dir),
            "--output-dir",
            str(output_dir),
            "--model-ids",
            str(model_ids),
        ],
    )

    with pytest.raises(SystemExit) as error:
        batch_ipsae.main()

    assert error.value.code == 1
    assert outside_guard.read_bytes() == original_guard
    run.assert_not_called()


def test_main_rejects_duplicate_model_ids(monkeypatch, tmp_path):
    pae_dir = tmp_path / "pae"
    pdb_dir = tmp_path / "pdb"
    output_dir = tmp_path / "output"
    pae_dir.mkdir()
    pdb_dir.mkdir()
    (pae_dir / "AF-TEST-meta_v1.json").write_text("{}")
    (pdb_dir / "AF-TEST-model_v1.pdb").write_text("END\n")
    model_ids = tmp_path / "model_ids.txt"
    model_ids.write_text("AF-TEST\nAF-TEST\n")
    run = Mock()
    monkeypatch.setattr(batch_ipsae, "run_ipsae_batch", run)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "batch_ipsae.py",
            "--pae-dir",
            str(pae_dir),
            "--pdb-dir",
            str(pdb_dir),
            "--output-dir",
            str(output_dir),
            "--model-ids",
            str(model_ids),
        ],
    )

    with pytest.raises(SystemExit) as error:
        batch_ipsae.main()

    assert error.value.code == 1
    run.assert_not_called()


@pytest.mark.parametrize("summary_kind", ["missing", "stale", "short"])
def test_main_rejects_invalid_success_summary(
    monkeypatch, tmp_path, summary_kind
):
    pae_dir = tmp_path / "pae"
    pdb_dir = tmp_path / "pdb"
    output_dir = tmp_path / "output"
    pae_dir.mkdir()
    pdb_dir.mkdir()
    output_dir.mkdir()
    model_ids = ["AF-ONE", "AF-TWO"] if summary_kind == "short" else ["AF-ONE"]
    for model_id in model_ids:
        (pae_dir / f"{model_id}-meta_v1.json").write_text("{}")
        (pdb_dir / f"{model_id}-model_v1.pdb").write_text("END\n")

    summary_csv = output_dir / "ipsae_summary.csv"
    if summary_kind == "stale":
        summary_csv.write_text("model,score\nAF-ONE,old\n")

    def fake_run(input_dir, output_summary, *args):
        if summary_kind == "short":
            output_summary.write_text("model,score\nAF-ONE,0.9\n")
        return subprocess.CompletedProcess([], 0, stdout="", stderr="")

    monkeypatch.setattr(batch_ipsae, "run_ipsae_batch", fake_run)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "batch_ipsae.py",
            "--pae-dir",
            str(pae_dir),
            "--pdb-dir",
            str(pdb_dir),
            "--output-dir",
            str(output_dir),
        ],
    )

    with pytest.raises(SystemExit) as error:
        batch_ipsae.main()

    assert error.value.code == 1
    if summary_kind == "stale":
        assert not summary_csv.exists()
