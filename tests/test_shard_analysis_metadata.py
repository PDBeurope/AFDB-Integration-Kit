import csv
import json

from afdb_integration_kit import analysis_metadata


def test_build_analysis_metadata_rows_from_ipsae_and_clash_outputs(tmp_path):
    work_dir = tmp_path / "work"
    ipsae_dir = work_dir / "ipsae"
    analysis_dir = work_dir / "clash_interface_analysis"
    ipsae_dir.mkdir(parents=True)
    analysis_dir.mkdir(parents=True)

    ipsae_csv = ipsae_dir / "ipsae_summary.csv"
    with ipsae_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "pdb_path",
                "pae_cutoff",
                "dist_cutoff",
                "iptm_af",
                "ipsae_AB",
                "ipsae_BA",
                "pDockQ2_AB",
                "pDockQ2_BA",
                "processing_time_ms",
            ],
        )
        writer.writeheader()
        writer.writerow({
            "pdb_path": "/tmp/AF-0001-model_v1.pdb",
            "pae_cutoff": "10",
            "dist_cutoff": "10",
            "iptm_af": "0.8",
            "ipsae_AB": "0.55",
            "ipsae_BA": "0.61",
            "pDockQ2_AB": "0.20",
            "pDockQ2_BA": "0.24",
            "processing_time_ms": "11",
        })

    (analysis_dir / "AF-0001-model_v1_clashes.json").write_text(
        json.dumps(
            {
                "sites": [
                    {
                        "label": "backbone_clashes",
                        "additional_site_annotations": {"n_clashes": 2},
                    },
                    {
                        "label": "heavy_atom_clashes",
                        "additional_site_annotations": {"n_clashes": 7},
                    },
                ]
            }
        )
    )
    (analysis_dir / "AF-0001-model_v1_interface.json").write_text(json.dumps({
        "sites": [
            {"additional_site_annotations": {"interactions": [{"a": 1}, {"a": 2}]}}
        ]
    }))

    rows = analysis_metadata._build_analysis_metadata_rows(
        work_dir=work_dir,
        run_name="test_run",
        original_ids=["AF-old"],
        model_ids=["AF-0001"],
        status_by_id={"AF-0001": "uploaded"},
        failure_reasons={},
        source_archive="heterodimers/chunk_0000.tar.lz4",
        shard_id=3,
        task_id=1,
        batch_id=2,
        batch_started_at="2026-04-29T00:00:00+00:00",
        batch_finished_at="2026-04-29T00:01:00+00:00",
        ipsae_threshold=0.6,
        pdockq2_threshold=0.23,
        include_scores=True,
    )

    assert len(rows) == 1
    row = rows[0]
    assert row["passes_quality_threshold"] == "true"
    assert row["ipsae_max"] == 0.61
    assert row["pdockq2_max"] == 0.24

    scores = json.loads(row["scores_json"])
    assert scores["iptm_af"] == 0.8
    assert scores["N_clash_backbone"] == 2
    assert scores["N_clash_heavyAtom"] == 7
    assert scores["N_interface_interactions"] == 2

    expected_files = json.loads(row["expected_output_files_json"])
    assert "AF-0001-model_v1.cif" in expected_files
    assert "AF-0001-model_v1_interface.json" in expected_files


def test_finalize_analysis_metadata_flattens_score_columns(tmp_path):
    tmp_path.mkdir(parents=True, exist_ok=True)
    csv_path = tmp_path / "analysis_metadata.csv"
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=analysis_metadata.ANALYSIS_METADATA_FIELDS
        )
        writer.writeheader()
        writer.writerow({
            "timestamp": "2026-04-29T00:00:00+00:00",
            "run_name": "test_run",
            "model_id": "AF-0001",
            "original_id": "AF-old",
            "upload_status": "uploaded",
            "passes_quality_threshold": "true",
            "ipsae_max": "0.61",
            "pdockq2_max": "0.24",
            "quality_ipsae_threshold": "0.6",
            "quality_pdockq2_threshold": "0.23",
            "source_archive": "batch-1.tar",
            "shard_id": "3",
            "task_id": "1",
            "batch_id": "2",
            "batch_started_at": "2026-04-29T00:00:00+00:00",
            "batch_finished_at": "2026-04-29T00:01:00+00:00",
            "failure_reason": "",
            "expected_output_files_json": "[]",
            "scores_json": json.dumps({
                "ipsae_AB": 0.61,
                "pDockQ2_BA": 0.24,
                "N_clash_backbone": 2,
            }),
        })

    score_keys = analysis_metadata._score_keys(csv_path)
    assert score_keys == ["N_clash_backbone", "ipsae_AB", "pDockQ2_BA"]

    with csv_path.open(newline="") as handle:
        row = next(csv.DictReader(handle))
    parsed = analysis_metadata._coerce_row(row, score_keys)

    assert parsed["score_ipsae_AB"] == 0.61
    assert parsed["score_pDockQ2_BA"] == 0.24
    assert parsed["score_N_clash_backbone"] == 2.0
