import gzip

import pytest

from genome_amr_predictor.data.fasta import inspect_fasta


def test_inspect_gzipped_fasta(tmp_path):
    path = tmp_path / "genome.fna.gz"
    with gzip.open(path, "wt", encoding="ascii") as stream:
        stream.write(">contig1\nACGTNN\n>contig2\nGGCC\n")

    metrics = inspect_fasta(path)

    assert metrics.total_length == 10
    assert metrics.contig_count == 2
    assert metrics.contig_n50 == 6
    assert metrics.gc_percent == 60.0
    assert metrics.ambiguous_bases == 2


def test_invalid_fasta_symbol_is_rejected(tmp_path):
    path = tmp_path / "bad.fna"
    path.write_text(">contig\nACGTZ\n", encoding="ascii")

    with pytest.raises(ValueError, match="Invalid FASTA symbols"):
        inspect_fasta(path)
