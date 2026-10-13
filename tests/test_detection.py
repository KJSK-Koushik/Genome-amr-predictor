from __future__ import annotations

import gzip
import json
import subprocess
from pathlib import Path
from types import SimpleNamespace

import pytest

from genome_amr_predictor.cli import main
from genome_amr_predictor.detection import (
    AMRFinderRuntime,
    detect_genome,
    parse_amrfinder,
    write_evidence,
)

HEADER = [
    "Protein id",
    "Contig id",
    "Start",
    "Stop",
    "Strand",
    "Element symbol",
    "Element name",
    "Scope",
    "Type",
    "Subtype",
    "Class",
    "Subclass",
    "Method",
    "Target length",
    "Reference sequence length",
    "% Coverage of reference",
    "% Identity to reference",
    "Alignment length",
    "Closest reference accession",
    "Closest reference name",
    "HMM accession",
    "HMM description",
    "Hierarchy node",
]
ROW = [
    "NA",
    "contig1",
    "1",
    "12",
    "-",
    "blaTEM-1",
    "TEM-1 beta-lactamase",
    "core",
    "AMR",
    "AMR",
    "BETA-LACTAM",
    "BETA-LACTAM",
    "ALLELEX",
    "4",
    "4",
    "100.00",
    "100.00",
    "4",
    "WP_fixture",
    "fixture reference",
    "NA",
    "NA",
    "blaTEM",
]


def table(path: Path, rows: list[list[str]], header: list[str] | None = None) -> Path:
    path.write_text("\n".join("\t".join(row) for row in [header or HEADER, *rows]) + "\n")
    return path


def test_allele_and_reverse_coordinates(tmp_path):
    hit = parse_amrfinder(table(tmp_path / "hits.tsv", [ROW]), "isolate1")[0]
    assert hit["allele_symbol"] == "blaTEM-1"
    assert (hit["start"], hit["stop"], hit["strand"]) == (1, 12, "-")
    assert hit["drug_class"] == "BETA-LACTAM"
    assert hit["hierarchy_nodes"] == ["blaTEM"]
    assert hit["mechanism"] is None
    assert hit["raw"]["Closest reference accession"] == "WP_fixture"


@pytest.mark.parametrize(
    "method,warning",
    [
        ("PARTIALX", "partial_hit"),
        ("PARTIAL_CONTIG_ENDX", "partial_hit"),
        ("INTERNAL_STOP", "internal_stop"),
        ("NEW_METHOD", "unrecognized_method"),
    ],
)
def test_uncertain_hits_are_retained_with_warning(tmp_path, method, warning):
    row = ROW.copy()
    row[12] = method
    hit = parse_amrfinder(table(tmp_path / "hits.tsv", [row]), "s")[0]
    assert warning in hit["warnings"]
    assert hit["allele_symbol"] is None


def test_point_and_non_amr_elements_are_distinct(tmp_path):
    point = ROW.copy()
    point[5], point[9], point[12] = "gyrA_S83L", "POINT", "POINTX"
    other = ROW.copy()
    other[7], other[8], other[9] = "plus", "STRESS", "METAL"
    hits = parse_amrfinder(table(tmp_path / "hits.tsv", [point, other]), "s")
    assert [hit["evidence_kind"] for hit in hits] == ["mutation", "other_element"]
    assert "plus_or_unrecognized_scope" in hits[1]["warnings"]


def test_legacy_headers_hmm_na_and_header_only(tmp_path):
    header = HEADER.copy()
    header[5], header[6] = "Gene symbol", "Sequence name"
    header[18], header[19] = "Accession of closest sequence", "Name of closest sequence"
    row = ROW.copy()
    row[12] = "HMM"
    row[14:20] = ["NA"] * 6
    hit = parse_amrfinder(table(tmp_path / "hits.tsv", [row], header), "s")[0]
    assert hit["identity_percent"] is None
    assert "Gene symbol" in hit["raw"]
    assert parse_amrfinder(table(tmp_path / "none.tsv", []), "s") == []


@pytest.mark.parametrize(
    "problem",
    [
        "empty",
        "truncated",
        "duplicate",
        "nan",
        "percent",
        "zero_start",
        "reverse_start",
        "missing_stop",
    ],
)
def test_invalid_table_is_not_a_no_hit_result(tmp_path, problem):
    header, row = HEADER.copy(), ROW.copy()
    if problem == "empty":
        path = tmp_path / "bad.tsv"
        path.write_text("")
    else:
        if problem == "truncated":
            row.pop()
        elif problem == "duplicate":
            header[-1] = header[0]
        elif problem == "nan":
            row[16] = "nan"
        elif problem == "percent":
            row[16] = "101"
        elif problem == "zero_start":
            row[2] = "0"
        elif problem == "reverse_start":
            row[2] = "13"
        elif problem == "missing_stop":
            row[3] = "NA"
        path = table(tmp_path / "bad.tsv", [row], header)
    with pytest.raises(ValueError):
        parse_amrfinder(path, "s")


class FakeRuntime(AMRFinderRuntime):
    """A deterministic subprocess boundary for software tests, never a biology validation."""

    def __init__(self, root: Path, rows=None, fail=False):
        executable = root / "amrfinder"
        executable.write_text("test executable")
        database = root / "db"
        database.mkdir()
        (database / "version.txt").write_text("fixture-db")
        super().__init__(executable, database)
        self.rows = [ROW] if rows is None else rows
        self.fail = fail

    def invoke(self, arguments, timeout):
        if arguments == ["--version"]:
            return subprocess.CompletedProcess(arguments, 0, "fixture-version\n", "")
        if "-V" in arguments:
            return subprocess.CompletedProcess(arguments, 0, "fixture-db\n", "")
        if self.fail:
            raise RuntimeError("fixture execution failure")
        assert arguments[arguments.index("-O") + 1] == "Escherichia"
        assert "--print_node" in arguments
        assert "--plus" not in arguments
        input_path = Path(arguments[arguments.index("-n") + 1])
        assert input_path.read_text().startswith(">contig1")
        table(Path(arguments[arguments.index("-o") + 1]), self.rows)
        return subprocess.CompletedProcess(arguments, 0, "", "fixture run log")


def test_compressed_fasta_run_and_provenance(tmp_path):
    fasta = tmp_path / "genome.fna.gz"
    with gzip.open(fasta, "wt") as stream:
        stream.write(">contig1\nACGTACGTACGT\n")
    output = tmp_path / "run"
    manifest = detect_genome(fasta, output, "s", FakeRuntime(tmp_path))
    assert manifest["status"] == "completed"
    assert manifest["software_version"] == "fixture-version"
    assert "version.txt" in manifest["database_files_sha256"]
    evidence = json.loads((output / "detections.json").read_text())
    assert evidence["summary"]["amr_hits"] == 1
    assert "phenotypic" in evidence["warning"]
    assert not any(key.startswith("phenotype") for key in evidence)


def test_failed_run_has_failed_manifest_and_no_success_evidence(tmp_path):
    fasta = tmp_path / "genome.fa"
    fasta.write_text(">contig1\nACGTACGTACGT\n")
    output = tmp_path / "run"
    with pytest.raises(RuntimeError, match="execution failure"):
        detect_genome(fasta, output, "s", FakeRuntime(tmp_path, fail=True))
    assert json.loads((output / "run_manifest.json").read_text())["status"] == "failed"
    assert not (output / "detections.json").exists()


def test_rejects_duplicate_contigs_and_stale_output(tmp_path):
    fasta = tmp_path / "genome.fa"
    fasta.write_text(">contig1\nACGT\n>contig1\nACGT\n")
    runtime = FakeRuntime(tmp_path)
    with pytest.raises(ValueError, match="unique"):
        detect_genome(fasta, tmp_path / "run", "s", runtime)
    fasta.write_text(">contig1\nACGTACGTACGT\n")
    output = tmp_path / "run"
    output.mkdir()
    (output / "stale.tsv").touch()
    with pytest.raises(ValueError, match="empty"):
        detect_genome(fasta, output, "s", runtime)


def test_no_hit_evidence_and_offline_cli(tmp_path):
    source = table(tmp_path / "empty.tsv", [])
    output = tmp_path / "evidence.json"
    assert (
        main(
            [
                "parse-detection",
                "--input",
                str(source),
                "--sample-id",
                "s",
                "--software-version",
                "fixture",
                "--database-version",
                "fixture",
                "--output",
                str(output),
            ]
        )
        == 0
    )
    document = json.loads(output.read_text())
    assert document["status"] == "completed"
    assert document["summary"]["total_hits"] == 0
    assert "susceptibility is unknown" in document["no_hit_interpretation"]
    assert document["provenance"]["version_source"] == "user_supplied_not_runtime_verified"
    with pytest.raises(SystemExit):
        main(
            [
                "parse-detection",
                "--input",
                str(source),
                "--sample-id",
                "s",
                "--software-version",
                "fixture",
                "--database-version",
                "fixture",
                "--output",
                str(output),
            ]
        )


def test_native_timeout_is_an_execution_failure(tmp_path, monkeypatch):
    runtime = FakeRuntime(tmp_path)

    def timeout(*args, **kwargs):
        raise subprocess.TimeoutExpired("amrfinder", 1)

    monkeypatch.setattr(subprocess, "run", timeout)
    with pytest.raises(RuntimeError, match="could not execute"):
        AMRFinderRuntime.invoke(runtime, ["--version"], 1)


def test_write_evidence_keeps_no_hit_warning(tmp_path):
    output = tmp_path / "evidence.json"
    write_evidence(output, [], {"software_version": "fixture", "database_version": "fixture"})
    assert json.loads(output.read_text())["hits"] == []


def test_official_release_expected_output():
    fixture = Path(__file__).parent / "fixtures" / "amrfinder_v4.2.7_dna.expected.tsv"
    hits = parse_amrfinder(fixture, "ncbi-release-test")
    assert len(hits) == 24
    methods = {hit["method"] for hit in hits}
    assert {"POINTN", "POINTX", "ALLELEX", "PARTIAL_CONTIG_ENDX", "INTERNAL_STOP"} <= methods
    assert any(hit["element_symbol"] == "pmrB_C84R" for hit in hits)


def test_wsl_filesystem_path_conversion_and_distribution_guard(tmp_path):
    runtime = FakeRuntime(tmp_path)
    runtime.wsl_distribution = "Ubuntu"
    path = SimpleNamespace(resolve=lambda: "\\\\wsl.localhost\\Ubuntu\\opt\\AMR tools\\amrfinder")
    assert runtime.runtime_path(path) == "/opt/AMR tools/amrfinder"
    runtime.wsl_distribution = "Debian"
    with pytest.raises(ValueError, match="different distribution"):
        runtime.runtime_path(path)
