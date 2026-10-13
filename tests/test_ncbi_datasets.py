import json

from genome_amr_predictor.data.ncbi_datasets import (
    AssemblyMetadata,
    NCBIDatasetsClient,
    resolve_assemblies,
    select_preferred_assembly,
)


def _report(accession: str, source: str = "SOURCE_DATABASE_REFSEQ") -> dict:
    return {
        "accession": accession,
        "paired_accession": "GCA_000001.1",
        "source_database": source,
        "organism": {"tax_id": 562, "organism_name": "Escherichia coli"},
        "assembly_info": {
            "assembly_level": "Scaffold",
            "assembly_name": "test",
            "assembly_status": "current",
            "bioproject_accession": "PRJNA1",
            "release_date": "2026-01-01",
            "biosample": {"accession": "SAMN1"},
        },
        "assembly_stats": {
            "total_sequence_length": "5000000",
            "number_of_contigs": 20,
            "contig_n50": 100000,
            "gc_percent": 50.5,
        },
        "checkm_info": {"completeness": 99.0, "contamination": 0.2},
        "average_nucleotide_identity": {"match_status": "species_match"},
    }


def test_assembly_metadata_and_refseq_selection():
    refseq = AssemblyMetadata.from_report(_report("GCF_000001.1"))
    genbank = AssemblyMetadata.from_report(
        _report("GCA_000001.1", source="SOURCE_DATABASE_GENBANK")
    )

    assert refseq.total_sequence_length == 5_000_000
    assert select_preferred_assembly([genbank, refseq]) == refseq


def test_resolution_cache_is_resumable(tmp_path, monkeypatch):
    cache_path = tmp_path / "assemblies.jsonl"
    cached = AssemblyMetadata.from_report(_report("GCF_000001.1"))
    cache_path.write_text(
        json.dumps({"query_biosample": "SAMN1", "assembly": cached.as_dict()}) + "\n",
        encoding="utf-8",
    )

    def fail_if_called(*args, **kwargs):
        raise AssertionError("network should not be called for a complete cache")

    monkeypatch.setattr(NCBIDatasetsClient, "assembly_reports", fail_if_called)
    resolved = resolve_assemblies(["SAMN1"], NCBIDatasetsClient(), cache_path, batch_size=50)

    assert resolved["SAMN1"] == cached
