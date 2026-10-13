from genome_amr_predictor.cohort_config import AssemblyQCConfig, Phase2Config, SplitConfig
from genome_amr_predictor.data.cohort import build_cohort, normalize_training_labels
from genome_amr_predictor.data.ncbi_datasets import AssemblyMetadata
from genome_amr_predictor.data.schema import ASTRecord


def _record(biosample: str, antibiotic: str, label: int) -> ASTRecord:
    return ASTRecord(
        biosample_accession=biosample,
        organism="Escherichia coli",
        taxonomy_id="562",
        bioproject_accession="",
        sra_accession="",
        isolate=biosample,
        collection_date="",
        geo_loc_name="",
        isolation_source="",
        host="",
        antibiotic=antibiotic,
        phenotype_raw="resistant" if label else "susceptible",
        phenotype_normalized="resistant" if label else "susceptible",
        binary_label=label,
        measurement_sign="",
        measurement="",
        measurement_units="",
        laboratory_method="",
        laboratory_platform="",
        vendor="",
        method_version_or_reagent="",
        testing_standard="CLSI",
    )


def _assembly(biosample: str, accession: str, project: str = "PRJNA1") -> AssemblyMetadata:
    return AssemblyMetadata(
        biosample_accession=biosample,
        assembly_accession=accession,
        paired_accession="",
        source_database="SOURCE_DATABASE_REFSEQ",
        organism="Escherichia coli",
        tax_id=562,
        assembly_level="Scaffold",
        assembly_name="test",
        assembly_status="current",
        bioproject_accession=project,
        release_date="2026-01-01",
        total_sequence_length=5_000_000,
        total_ungapped_length=5_000_000,
        number_of_contigs=20,
        contig_n50=100_000,
        gc_percent=50.5,
        sequencing_technology="Illumina",
        assembly_method="test",
        checkm_completeness=99.0,
        checkm_contamination=0.2,
        ani_taxonomy_check_status="OK",
        ani_match_status="species_match",
        annotation_name="annotation",
    )


def _config() -> Phase2Config:
    return Phase2Config(
        species="Escherichia coli",
        antibiotics=("ampicillin",),
        assembly_request_batch_size=50,
        assembly_qc=AssemblyQCConfig(
            minimum_length=4_000_000,
            maximum_length=6_500_000,
            maximum_contigs=500,
            minimum_contig_n50=20_000,
            minimum_checkm_completeness=95.0,
            maximum_checkm_contamination=5.0,
            require_checkm=False,
            reject_ani_mismatch=True,
        ),
        splits=SplitConfig(
            seed="test",
            train_fraction=0.7,
            validation_fraction=0.15,
            test_fraction=0.15,
        ),
    )


def test_conflicts_and_duplicate_assemblies_are_removed():
    normalized = normalize_training_labels(
        [
            _record("SAMN1", "ampicillin", 1),
            _record("SAMN1", "ampicillin", 1),
            _record("SAMN2", "ampicillin", 0),
            _record("SAMN3", "ampicillin", 0),
            _record("SAMN3", "ampicillin", 1),
        ],
        ("ampicillin",),
    )
    bundle = build_cohort(
        normalized,
        {
            "SAMN1": _assembly("SAMN1", "GCF_1"),
            "SAMN2": _assembly("SAMN2", "GCF_1"),
        },
        _config(),
    )

    assert normalized.duplicate_rows_removed == 2
    assert normalized.conflicting_pairs_removed == 1
    assert len(bundle.cohort_rows) == 1
    assert bundle.summary["exclusion_reasons"]["duplicate_assembly_accession"] == 1


def test_all_samples_in_a_project_receive_the_same_frozen_split():
    records = [_record(f"SAMN{index}", "ampicillin", index % 2) for index in range(1, 10)]
    normalized = normalize_training_labels(records, ("ampicillin",))
    assemblies = {
        f"SAMN{index}": _assembly(
            f"SAMN{index}",
            f"GCF_{index}",
            project=f"PRJNA{(index - 1) // 3}",
        )
        for index in range(1, 10)
    }

    first = build_cohort(normalized, assemblies, _config())
    second = build_cohort(normalized, assemblies, _config())

    project_splits = {}
    for row in first.assembly_rows:
        project_splits.setdefault(row["bioproject_accession"], set()).add(row["split"])
    assert all(len(splits) == 1 for splits in project_splits.values())
    assert [row["split"] for row in first.assembly_rows] == [
        row["split"] for row in second.assembly_rows
    ]
