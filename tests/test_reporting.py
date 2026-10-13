from __future__ import annotations

import json

import pytest

from genome_amr_predictor import engine, reporting


@pytest.fixture
def completed(setup, tmp_path):
    fasta, config, cache, _ = setup
    engine.analyze_genome(fasta, tmp_path / "run", "sample", config, detection_cache=cache)
    return tmp_path / "run"


def test_load_run_verifies_artifacts(completed):
    run = reporting.load_run(completed)
    assert run["status"] == "completed"
    (completed / "predictions.csv").write_text("antibiotic,prediction\nampicillin,S\n")
    with pytest.raises(ValueError, match=r"predictions\.csv changed"):
        reporting.load_run(completed)


def test_load_run_rejects_unfinished_and_foreign(completed, tmp_path):
    manifest = json.loads((completed / "run_manifest.json").read_text())
    manifest["stage"] = "detection"
    (completed / "run_manifest.json").write_text(json.dumps(manifest))
    with pytest.raises(ValueError, match="not finished"):
        reporting.load_run(completed)
    with pytest.raises(ValueError, match="Not an integrated"):
        reporting.load_run(_foreign(tmp_path))


def _foreign(tmp_path):
    directory = tmp_path / "foreign"
    directory.mkdir()
    (directory / "run_manifest.json").write_text(json.dumps({"schema_version": "other"}))
    return directory


def test_find_runs_skips_detection_manifests(setup, tmp_path):
    fasta, config, cache, _ = setup
    engine.analyze_genome(fasta, tmp_path / "runs/a", "sample", config, detection_cache=cache)
    (tmp_path / "runs/a/nested").mkdir()
    (tmp_path / "runs/a/nested/run_manifest.json").write_text(
        (cache / "run_manifest.json").read_text()
    )
    assert reporting.find_runs(tmp_path / "runs", tmp_path / "missing") == [tmp_path / "runs/a"]


def test_detection_failure_is_unknown_not_absent(setup, tmp_path):
    fasta, config, cache, _ = setup
    (cache / "amrfinder.tsv").write_text("altered")
    engine.analyze_genome(fasta, tmp_path / "run", "sample", config, detection_cache=cache)
    run = reporting.load_run(tmp_path / "run")
    summary = reporting.detection_summary(run)
    assert "presence/absence is unknown" in summary
    assert "No known determinants" not in summary
    assert reporting.determinant_classes(run) == {}
    pdf = reporting.render_pdf(run, compress=False)
    assert pdf.startswith(b"%PDF")
    assert b"presence/absence is unknown" in pdf


def test_no_hits_does_not_claim_susceptibility():
    run = {"detection": {"status": "completed", "evidence": {"summary": {"total_hits": 0}}}}
    assert "does not establish susceptibility" in reporting.detection_summary(run)


@pytest.mark.parametrize(
    ("row", "note"),
    [
        ({"status": "unsupported", "resistance_score": None, "features_present": None}, "No valid"),
        ({"status": "withheld_qc", "resistance_score": None, "features_present": 5}, "Withheld"),
        ({"status": "predicted", "resistance_score": 0.2, "features_present": 0}, "baseline"),
        ({"status": "predicted", "resistance_score": 0.52, "features_present": 9}, "threshold"),
        ({"status": "predicted", "resistance_score": 0.9, "features_present": 9}, ""),
    ],
)
def test_prediction_notes(row, note):
    run = {"phenotype_prediction": {"result": {"predictions": [{"antibiotic": "x", **row}]}}}
    [result] = reporting.prediction_rows(run)
    assert note in result["note"] and (note or result["note"] == "")


def test_pdf_contains_both_outputs_and_limits(completed):
    run = reporting.load_run(completed)
    pdf = reporting.render_pdf(run, compress=False)
    assert pdf.startswith(b"%PDF")
    for text in (b"ampicillin", b"Known resistance determinants", b"laboratory confirmation"):
        assert text in pdf
    assert reporting.determinant_classes(run)


def test_pdf_for_rejected_input(setup, tmp_path):
    fasta, config, _, _ = setup
    fasta.write_text(">dup\nACGT\n>dup\nACGT\n")
    engine.analyze_genome(fasta, tmp_path / "run", "sample", config)
    run = reporting.load_run(tmp_path / "run")
    pdf = reporting.render_pdf(run, compress=False)
    assert b"Input was rejected" in pdf


def test_cli_export_pdf(completed, tmp_path):
    from genome_amr_predictor.cli import main

    output = tmp_path / "report.pdf"
    assert main(["export-pdf", "--run-dir", str(completed), "--output", str(output)]) == 0
    assert output.read_bytes().startswith(b"%PDF")
    (completed / "report.md").write_text("edited")
    with pytest.raises(SystemExit):
        main(["export-pdf", "--run-dir", str(completed), "--output", str(tmp_path / "x.pdf")])
    assert not (tmp_path / "x.pdf").exists()
