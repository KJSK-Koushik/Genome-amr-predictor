from __future__ import annotations

from pathlib import Path

import pytest

from genome_amr_predictor import dashboard, engine

pytest.importorskip("streamlit")
from streamlit.testing.v1 import AppTest

APP = str(Path(dashboard.__file__))


def _app(run_dir: Path | None = None) -> AppTest:
    app = AppTest.from_file(APP, default_timeout=60)
    if run_dir is not None:
        app.session_state["run_dir"] = str(run_dir)
    return app


@pytest.mark.parametrize(
    ("name", "expected"),
    [("genome.FNA", "input.fna"), ("x.fasta.gz", "input.fna.gz"), ("../../a.fa", "input.fna")],
)
def test_upload_path_never_uses_user_name(tmp_path, name, expected):
    path = dashboard.upload_path(name, tmp_path)
    assert path.name == expected
    assert path.parent.parent == tmp_path


def test_upload_rejects_non_fasta(tmp_path):
    with pytest.raises(ValueError, match="FASTA"):
        dashboard.upload_path("reads.fastq.gz", tmp_path)


def test_safe_name():
    assert dashboard.safe_name("SAMN 1/../x") == "SAMN_1_.._x"
    assert dashboard.safe_name("///") == "sample"


def test_landing_page(monkeypatch, tmp_path):
    monkeypatch.setenv(dashboard.RUN_ROOTS_ENV, str(tmp_path))
    app = _app().run()
    assert not app.exception
    assert any("not a clinical diagnosis" in item.value for item in app.info)


def test_completed_run_renders_with_downloads(setup, tmp_path, monkeypatch):
    fasta, config, cache, _ = setup
    engine.analyze_genome(fasta, tmp_path / "run", "S1", config, detection_cache=cache)
    monkeypatch.setenv(dashboard.RUN_ROOTS_ENV, str(tmp_path))
    app = _app(tmp_path / "run").run()
    assert not app.exception
    assert app.title[0].value == "Sample S1"
    assert any("completed" in item.value for item in app.success)
    assert [tab.label for tab in app.tabs][:2] == ["Phenotype predictions", "Known determinants"]
    assert len(app.get("download_button")) == 5


def test_failed_detection_shows_unknown(setup, tmp_path, monkeypatch):
    fasta, config, cache, _ = setup
    (cache / "detections.json").write_text("altered")
    engine.analyze_genome(fasta, tmp_path / "run", "S1", config, detection_cache=cache)
    monkeypatch.setenv(dashboard.RUN_ROOTS_ENV, str(tmp_path))
    app = _app(tmp_path / "run").run()
    assert not app.exception
    assert any("presence/absence is unknown" in item.value for item in app.error)
    assert any("partial" in item.value for item in app.warning)


def test_tampered_run_is_refused(setup, tmp_path, monkeypatch):
    fasta, config, cache, _ = setup
    engine.analyze_genome(fasta, tmp_path / "run", "S1", config, detection_cache=cache)
    (tmp_path / "run/detections.csv").write_text("edited")
    monkeypatch.setenv(dashboard.RUN_ROOTS_ENV, str(tmp_path))
    app = _app(tmp_path / "run").run()
    assert any("integrity check" in item.value for item in app.error)
    assert not app.title or app.title[0].value != "Sample S1"


def test_run_title(setup, tmp_path):
    fasta, config, cache, _ = setup
    engine.analyze_genome(fasta, tmp_path / "run", "S1", config, detection_cache=cache)
    assert dashboard.run_title(tmp_path / "run").startswith("S1 · completed · ")
    assert dashboard.run_title(tmp_path / "missing").endswith("missing")
