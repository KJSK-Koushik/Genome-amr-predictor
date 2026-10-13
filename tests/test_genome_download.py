import gzip
import io
import json
import zipfile

import pytest

from genome_amr_predictor.data.genome_download import (
    _extract_safely,
    _fetch_ftp_fasta,
    _ftp_url,
    _inspect_downloaded_accession,
    _prefetch_ftp,
    download_genome_shards,
)


def test_completed_shard_is_reused_without_subprocess(tmp_path):
    executable = tmp_path / "datasets"
    executable.write_text("placeholder", encoding="utf-8")
    shard = tmp_path / "downloads" / "shard-0001"
    shard.mkdir(parents=True)
    row = {
        "assembly_accession": "GCF_1",
        "download_status": "verified",
        "sequence_sha256": "abc",
    }
    (shard / "genome_metrics.jsonl").write_text(json.dumps(row) + "\n", encoding="utf-8")
    (shard / "complete.json").write_text("{}\n", encoding="utf-8")

    rows = download_genome_shards(["GCF_1"], tmp_path / "downloads", executable, shard_size=1)

    assert rows == [row]
    assert (tmp_path / "downloads" / "genome_manifest.csv").is_file()


def test_cached_absolute_fasta_path_is_made_portable(tmp_path):
    executable = tmp_path / "datasets"
    executable.write_text("placeholder", encoding="utf-8")
    downloads = tmp_path / "downloads"
    shard = downloads / "shard-0001"
    fasta = shard / "package" / "ncbi_dataset" / "data" / "GCF_1" / "x.fna.gz"
    fasta.parent.mkdir(parents=True)
    row = {
        "assembly_accession": "GCF_1",
        "download_status": "verified",
        "sequence_sha256": "abc",
        "path": str(fasta.resolve()),
    }
    (shard / "genome_metrics.jsonl").write_text(json.dumps(row) + "\n", encoding="utf-8")
    (shard / "complete.json").write_text("{}\n", encoding="utf-8")

    rows = download_genome_shards(["GCF_1"], downloads, executable, shard_size=1)

    assert rows[0]["path"] == "shard-0001/package/ncbi_dataset/data/GCF_1/x.fna.gz"


def test_unsafe_zip_member_is_rejected(tmp_path):
    archive_path = tmp_path / "unsafe.zip"
    with zipfile.ZipFile(archive_path, "w") as archive:
        archive.writestr("../escape.txt", "unsafe")

    with pytest.raises(ValueError, match="Unsafe path"):
        _extract_safely(archive_path, tmp_path / "output")


def test_invalid_worker_count_is_rejected(tmp_path):
    executable = tmp_path / "datasets"
    executable.write_text("placeholder", encoding="utf-8")

    with pytest.raises(ValueError, match="between 1 and 30"):
        download_genome_shards(["GCF_1"], tmp_path / "output", executable, max_workers=0)


def test_truncated_gzip_is_reported_instead_of_crashing(tmp_path):
    import gzip

    assembly = tmp_path / "ncbi_dataset" / "data" / "GCF_1"
    assembly.mkdir(parents=True)
    fasta = assembly / "sample_genomic.fna.gz"
    fasta.write_bytes(gzip.compress(b">sample\nACGT\n")[:-5])
    row = _inspect_downloaded_accession(("GCF_1", tmp_path, tmp_path))
    assert row["download_status"] == "invalid_fasta"


def test_changed_accessions_cannot_reuse_old_shard(tmp_path):
    executable = tmp_path / "datasets"
    executable.write_text("placeholder", encoding="utf-8")
    shard = tmp_path / "downloads" / "shard-0001"
    shard.mkdir(parents=True)
    (shard / "complete.json").write_text("{}", encoding="utf-8")
    (shard / "genome_metrics.jsonl").write_text(
        json.dumps({"assembly_accession": "GCF_1"}) + "\n", encoding="utf-8"
    )
    with pytest.raises(ValueError, match="Cached accession list differs"):
        download_genome_shards(["GCF_2"], tmp_path / "downloads", executable, shard_size=1)


def test_ftp_download_is_verified_and_installed_atomically(tmp_path, monkeypatch):
    payload = gzip.compress(b">software-test\nACGT\n")
    monkeypatch.setattr(
        "genome_amr_predictor.data.genome_download.urllib.request.urlopen",
        lambda *args, **kwargs: io.BytesIO(payload),
    )
    target = tmp_path / "sample_genomic.fna.gz"
    assert _fetch_ftp_fasta(("https://ftp.ncbi.nlm.nih.gov/example", target))
    assert gzip.decompress(target.read_bytes()) == b">software-test\nACGT\n"
    assert not target.with_name(target.name + ".download.gz").exists()


def test_ftp_failure_does_not_install_invalid_fasta(tmp_path, monkeypatch):
    monkeypatch.setattr(
        "genome_amr_predictor.data.genome_download.urllib.request.urlopen",
        lambda *args, **kwargs: io.BytesIO(b"invalid gzip"),
    )
    monkeypatch.setattr("genome_amr_predictor.data.genome_download.time.sleep", lambda _: None)
    target = tmp_path / "sample_genomic.fna.gz"
    assert not _fetch_ftp_fasta(("https://ftp.ncbi.nlm.nih.gov/example", target))
    assert not target.exists()


def test_ftp_url_uses_ncbi_accession_directory_layout():
    assert _ftp_url("GCA_012603675.1", "GCA_012603675.1_PDT000644902.1_genomic.fna") == (
        "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/012/603/675/"
        "GCA_012603675.1_PDT000644902.1/"
        "GCA_012603675.1_PDT000644902.1_genomic.fna.gz"
    )


def test_ftp_prefetch_marks_complete_only_for_all_requested_accessions(tmp_path, monkeypatch):
    dataset = tmp_path / "ncbi_dataset"
    dataset.mkdir()
    (dataset / "fetch.txt").write_text(
        "https://example.org\t0\tdata/GCA_012603675.1/GCA_012603675.1_PDT000644902.1_genomic.fna\n",
        encoding="utf-8",
    )
    monkeypatch.setattr(
        "genome_amr_predictor.data.genome_download._fetch_ftp_fasta", lambda _: True
    )
    assert _prefetch_ftp(tmp_path, ["GCA_012603675.1"], 1)
    assert not _prefetch_ftp(tmp_path, ["GCA_012603675.1", "GCA_012603676.1"], 1)
