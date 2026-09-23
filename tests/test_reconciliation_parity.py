from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import pytest


REPO_ROOT = Path(__file__).resolve().parent.parent
MODULE_PATH = REPO_ROOT / "scripts" / "run_homodimer_reconciliation.py"


def _load_runner():
    spec = importlib.util.spec_from_file_location("reconciliation_parity", MODULE_PATH)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


runner = _load_runner()


def _write(root: Path, relative_path: str, content: bytes) -> None:
    path = root / relative_path
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(content)


def test_inventory_output_hashes_all_retained_regular_files(tmp_path: Path) -> None:
    retained = [
        "modelcif/a.cif",
        "modelpdb/a.pdb",
        "bcif/a.bcif",
        "modelcif_input/a.json",
        "manifests/a.csv",
        "batches/a.txt",
        "jsons/a.json",
        "scores/a.json",
        "ipsae/a.csv",
        "clash/a.json",
        "DSSP/a.cif",
        "other/nested/a.dat",
    ]
    for relative_path in retained:
        _write(tmp_path, relative_path, relative_path.encode())
    for relative_path in [
        "input/raw.pdb",
        "logs/run.log",
        ".cache/state",
        "other/input/raw.pdb",
        "pipeline_results.json",
        "report.html",
        ".DS_Store",
        "other/Thumbs.db",
    ]:
        _write(tmp_path, relative_path, b"ignored")

    artifacts = runner.inventory_output(tmp_path)

    assert [item.relative_path for item in artifacts] == sorted(retained)
    assert all(item.size == len(item.relative_path.encode()) for item in artifacts)
    assert all(len(item.sha256) == 64 for item in artifacts)
    assert runner.counts_by_top_level(artifacts) == {
        "DSSP": 1,
        "batches": 1,
        "bcif": 1,
        "clash": 1,
        "ipsae": 1,
        "jsons": 1,
        "manifests": 1,
        "modelcif": 1,
        "modelcif_input": 1,
        "modelpdb": 1,
        "other": 1,
        "scores": 1,
    }


def test_compare_outputs_accepts_identical_inventories(tmp_path: Path) -> None:
    baseline = tmp_path / "baseline"
    candidate = tmp_path / "candidate"
    _write(baseline, "modelcif/model.cif", b"same")
    _write(candidate, "modelcif/model.cif", b"same")

    result = runner.compare_outputs(
        runner.inventory_output(baseline), runner.inventory_output(candidate)
    )

    assert result.identical
    assert result.added == ()
    assert result.removed == ()
    assert result.changed == ()


def test_compare_outputs_reports_added_path(tmp_path: Path) -> None:
    _write(tmp_path, "candidate/jsons/new.json", b"new")

    result = runner.compare_outputs(
        runner.inventory_output(tmp_path / "baseline"),
        runner.inventory_output(tmp_path / "candidate"),
    )

    assert not result.identical
    assert result.added == ("jsons/new.json",)


def test_compare_outputs_reports_removed_path(tmp_path: Path) -> None:
    _write(tmp_path, "baseline/scores/old.json", b"old")

    result = runner.compare_outputs(
        runner.inventory_output(tmp_path / "baseline"),
        runner.inventory_output(tmp_path / "candidate"),
    )

    assert not result.identical
    assert result.removed == ("scores/old.json",)


def test_compare_outputs_reports_rename_as_removed_and_added(tmp_path: Path) -> None:
    _write(tmp_path, "baseline/modelpdb/old.pdb", b"same bytes")
    _write(tmp_path, "candidate/modelpdb/new.pdb", b"same bytes")

    result = runner.compare_outputs(
        runner.inventory_output(tmp_path / "baseline"),
        runner.inventory_output(tmp_path / "candidate"),
    )

    assert not result.identical
    assert result.removed == ("modelpdb/old.pdb",)
    assert result.added == ("modelpdb/new.pdb",)


def test_compare_outputs_reports_byte_changed_path(tmp_path: Path) -> None:
    _write(tmp_path, "baseline/bcif/model.bcif", b"before")
    _write(tmp_path, "candidate/bcif/model.bcif", b"after")

    result = runner.compare_outputs(
        runner.inventory_output(tmp_path / "baseline"),
        runner.inventory_output(tmp_path / "candidate"),
    )

    assert not result.identical
    assert result.changed == ("bcif/model.bcif",)


def test_compare_outputs_accepts_only_approved_ipsae_operational_drift(
    tmp_path: Path,
) -> None:
    header = "pdb_path,scientific_score,processing_time_ms\n"
    _write(
        tmp_path,
        "baseline/ipsae/ipsae_summary.csv",
        (header + "/run/baseline/model.pdb,0.921113,33.2319\n").encode(),
    )
    _write(
        tmp_path,
        "candidate/ipsae/ipsae_summary.csv",
        (header + "/run/candidate/model.pdb,0.921113,31.3105\n").encode(),
    )

    result = runner.compare_outputs(
        runner.inventory_output(tmp_path / "baseline"),
        runner.inventory_output(tmp_path / "candidate"),
        tmp_path / "baseline",
        tmp_path / "candidate",
    )

    assert result.identical
    assert result.raw_changed == ("ipsae/ipsae_summary.csv",)
    assert result.changed == ()
    assert result.accepted_operational_drift == ("ipsae/ipsae_summary.csv",)
    assert result.policy_errors == ()


def test_compare_outputs_rejects_ipsae_scientific_column_change(
    tmp_path: Path,
) -> None:
    header = "pdb_path,scientific_score,processing_time_ms\n"
    _write(
        tmp_path,
        "baseline/ipsae/ipsae_summary.csv",
        (header + "/run/a/model.pdb,0.921113,33.2\n").encode(),
    )
    _write(
        tmp_path,
        "candidate/ipsae/ipsae_summary.csv",
        (header + "/run/b/model.pdb,0.921114,31.3\n").encode(),
    )

    result = runner.compare_outputs(
        runner.inventory_output(tmp_path / "baseline"),
        runner.inventory_output(tmp_path / "candidate"),
        tmp_path / "baseline",
        tmp_path / "candidate",
    )

    assert not result.identical
    assert result.raw_changed == ("ipsae/ipsae_summary.csv",)
    assert result.changed == ("ipsae/ipsae_summary.csv",)
    assert result.accepted_operational_drift == ()


@pytest.mark.parametrize(
    ("baseline_score", "candidate_score"),
    [
        ("1", '"1"'),
        ("1", " 1"),
        ("1", "1 "),
    ],
)
def test_compare_outputs_preserves_raw_nonapproved_cell_bytes(
    tmp_path: Path, baseline_score: str, candidate_score: str
) -> None:
    header = "pdb_path,scientific_score,processing_time_ms\n"
    _write(
        tmp_path,
        "baseline/ipsae/ipsae_summary.csv",
        (header + f"/run/a,{baseline_score},33.2\n").encode(),
    )
    _write(
        tmp_path,
        "candidate/ipsae/ipsae_summary.csv",
        (header + f"/run/b,{candidate_score},31.3\n").encode(),
    )

    result = runner.compare_outputs(
        runner.inventory_output(tmp_path / "baseline"),
        runner.inventory_output(tmp_path / "candidate"),
        tmp_path / "baseline",
        tmp_path / "candidate",
    )

    assert not result.identical
    assert result.changed == ("ipsae/ipsae_summary.csv",)
    assert result.accepted_operational_drift == ()


def test_compare_outputs_rejects_lf_vs_crlf(tmp_path: Path) -> None:
    _write(
        tmp_path,
        "baseline/ipsae/ipsae_summary.csv",
        b"pdb_path,score,processing_time_ms\n/run/a,1,2\n",
    )
    _write(
        tmp_path,
        "candidate/ipsae/ipsae_summary.csv",
        b"pdb_path,score,processing_time_ms\r\n/run/b,1,3\r\n",
    )

    result = runner.compare_outputs(
        runner.inventory_output(tmp_path / "baseline"),
        runner.inventory_output(tmp_path / "candidate"),
        tmp_path / "baseline",
        tmp_path / "candidate",
    )

    assert not result.identical
    assert result.changed == ("ipsae/ipsae_summary.csv",)
    assert any("header bytes differ" in error for error in result.policy_errors)


def test_compare_outputs_rejects_quoted_vs_unquoted_header(tmp_path: Path) -> None:
    _write(
        tmp_path,
        "baseline/ipsae/ipsae_summary.csv",
        b"pdb_path,score,processing_time_ms\n/run/a,1,2\n",
    )
    _write(
        tmp_path,
        "candidate/ipsae/ipsae_summary.csv",
        b'"pdb_path",score,processing_time_ms\n/run/b,1,3\n',
    )

    result = runner.compare_outputs(
        runner.inventory_output(tmp_path / "baseline"),
        runner.inventory_output(tmp_path / "candidate"),
        tmp_path / "baseline",
        tmp_path / "candidate",
    )

    assert not result.identical
    assert any("ambiguous header dialect" in error for error in result.policy_errors)


def test_compare_outputs_rejects_quoted_approved_cell(tmp_path: Path) -> None:
    header = b"pdb_path,score,processing_time_ms\n"
    _write(
        tmp_path,
        "baseline/ipsae/ipsae_summary.csv",
        header + b'"/run/a",1,2\n',
    )
    _write(
        tmp_path,
        "candidate/ipsae/ipsae_summary.csv",
        header + b'"/run/b",1,3\n',
    )

    result = runner.compare_outputs(
        runner.inventory_output(tmp_path / "baseline"),
        runner.inventory_output(tmp_path / "candidate"),
        tmp_path / "baseline",
        tmp_path / "candidate",
    )

    assert not result.identical
    assert any("ambiguous data-row dialect" in error for error in result.policy_errors)


@pytest.mark.parametrize(
    ("baseline_csv", "candidate_csv", "error_text"),
    [
        (
            "pdb_path,score,processing_time_ms\n/a,1,2\n",
            "pdb_path,renamed_score,processing_time_ms\n/b,1,3\n",
            "header bytes differ",
        ),
        (
            'pdb_path,score,processing_time_ms\n"/a,1,2\n',
            "pdb_path,score,processing_time_ms\n/b,1,3\n",
            "malformed CSV",
        ),
        (
            "pdb_path,score\n/a,1\n",
            "pdb_path,score\n/b,1\n",
            "missing canonicalized column",
        ),
        (
            "pdb_path,pdb_path,score,processing_time_ms\n/a,/a,1,2\n",
            "pdb_path,pdb_path,score,processing_time_ms\n/b,/b,1,3\n",
            "duplicate canonicalized column",
        ),
    ],
)
def test_compare_outputs_rejects_invalid_ipsae_policy_inputs(
    tmp_path: Path,
    baseline_csv: str,
    candidate_csv: str,
    error_text: str,
) -> None:
    _write(
        tmp_path,
        "baseline/ipsae/ipsae_summary.csv",
        baseline_csv.encode(),
    )
    _write(
        tmp_path,
        "candidate/ipsae/ipsae_summary.csv",
        candidate_csv.encode(),
    )

    result = runner.compare_outputs(
        runner.inventory_output(tmp_path / "baseline"),
        runner.inventory_output(tmp_path / "candidate"),
        tmp_path / "baseline",
        tmp_path / "candidate",
    )

    assert not result.identical
    assert result.changed == ("ipsae/ipsae_summary.csv",)
    assert any(error_text in error for error in result.policy_errors)


def test_compare_outputs_does_not_normalize_same_named_csv_elsewhere(
    tmp_path: Path,
) -> None:
    header = "pdb_path,score,processing_time_ms\n"
    _write(
        tmp_path,
        "baseline/other/ipsae_summary.csv",
        (header + "/run/a,1,2\n").encode(),
    )
    _write(
        tmp_path,
        "candidate/other/ipsae_summary.csv",
        (header + "/run/b,1,3\n").encode(),
    )

    result = runner.compare_outputs(
        runner.inventory_output(tmp_path / "baseline"),
        runner.inventory_output(tmp_path / "candidate"),
        tmp_path / "baseline",
        tmp_path / "candidate",
    )

    assert not result.identical
    assert result.changed == ("other/ipsae_summary.csv",)
    assert result.accepted_operational_drift == ()


def test_source_commit_reads_archive_provenance_marker(tmp_path: Path) -> None:
    (tmp_path / ".source_commit").write_text("0a89d48\n", encoding="utf-8")

    assert runner.source_commit(tmp_path) == "0a89d48"


def test_inventory_fixture_includes_input_files(tmp_path: Path) -> None:
    _write(tmp_path, "config/model_ids.txt", b"model\n")
    _write(tmp_path, "input/model.pdb", b"ATOM\n")

    records = runner.inventory_fixture(tmp_path)

    assert [item.relative_path for item in records] == [
        "config/model_ids.txt",
        "input/model.pdb",
    ]


def _main_args(tmp_path: Path, image: Path, digest: str) -> list[str]:
    for name in ("baseline", "candidate"):
        repo = tmp_path / name
        repo.mkdir()
        (repo / ".source_commit").write_text(f"{name}-commit\n", encoding="utf-8")
    fixture = tmp_path / "fixture"
    fixture.mkdir()
    return [
        "--baseline-repo",
        str(tmp_path / "baseline"),
        "--candidate-repo",
        str(tmp_path / "candidate"),
        "--fixture",
        str(fixture),
        "--work-dir",
        str(tmp_path / "work"),
        "--report-json",
        str(tmp_path / "report.json"),
        "--runtime-image-path",
        str(image),
        "--runtime-image-sha256",
        digest,
    ]


def test_main_verifies_matching_runtime_image_before_pipelines(
    tmp_path: Path, monkeypatch
) -> None:
    image = tmp_path / "runtime.sqsh"
    image.write_bytes(b"pinned image")
    expected = runner._sha256(image)
    calls = []
    monkeypatch.setattr(
        runner,
        "run_pipeline",
        lambda repo, fixture, output: calls.append((repo, fixture, output)),
    )

    exit_code = runner.main(_main_args(tmp_path, image, expected))

    report = json.loads((tmp_path / "report.json").read_text(encoding="utf-8"))
    assert exit_code == 0
    assert len(calls) == 2
    assert report["runtime_image"] == {
        "actual_sha256": expected,
        "expected_sha256": expected,
        "match": True,
        "path": str(image),
        "readable_regular_file": True,
        "size": len(b"pinned image"),
    }


def test_main_rejects_runtime_image_digest_mismatch_before_pipelines(
    tmp_path: Path, monkeypatch
) -> None:
    image = tmp_path / "runtime.sqsh"
    image.write_bytes(b"wrong bytes")
    calls = []
    monkeypatch.setattr(runner, "run_pipeline", lambda *args: calls.append(args))

    exit_code = runner.main(_main_args(tmp_path, image, "0" * 64))

    report = json.loads((tmp_path / "report.json").read_text(encoding="utf-8"))
    assert exit_code == 1
    assert calls == []
    assert report["runtime_image"]["match"] is False
    assert report["integrity_failure"] == "runtime image SHA-256 mismatch"


def test_main_rejects_missing_runtime_image_before_pipelines(
    tmp_path: Path, monkeypatch
) -> None:
    image = tmp_path / "missing.sqsh"
    calls = []
    monkeypatch.setattr(runner, "run_pipeline", lambda *args: calls.append(args))

    exit_code = runner.main(_main_args(tmp_path, image, "0" * 64))

    report = json.loads((tmp_path / "report.json").read_text(encoding="utf-8"))
    assert exit_code == 1
    assert calls == []
    assert report["runtime_image"]["readable_regular_file"] is False
    assert report["runtime_image"]["actual_sha256"] is None
    assert report["integrity_failure"] == "runtime image is not a readable regular file"


def test_main_rejects_malformed_runtime_image_digest_before_pipelines(
    tmp_path: Path, monkeypatch
) -> None:
    image = tmp_path / "runtime.sqsh"
    image.write_bytes(b"pinned image")
    calls = []
    monkeypatch.setattr(runner, "run_pipeline", lambda *args: calls.append(args))

    exit_code = runner.main(_main_args(tmp_path, image, "not-a-sha256"))

    report = json.loads((tmp_path / "report.json").read_text(encoding="utf-8"))
    assert exit_code == 1
    assert calls == []
    assert report["runtime_image"]["expected_sha256"] == "not-a-sha256"
    assert report["runtime_image"]["match"] is False
    assert report["integrity_failure"] == "invalid expected SHA-256 format"
