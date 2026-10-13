import pytest

from genome_amr_predictor.data.labels import normalize_antibiotic, normalize_phenotype


@pytest.mark.parametrize(
    ("raw", "expected"),
    [
        ("R", ("resistant", 1)),
        (" resistant ", ("resistant", 1)),
        ("Sensitive", ("susceptible", 0)),
        ("S", ("susceptible", 0)),
        ("I", ("intermediate", None)),
        ("nonsusceptible", ("excluded", None)),
        (None, ("excluded", None)),
    ],
)
def test_normalize_phenotype(raw, expected):
    assert normalize_phenotype(raw) == expected


def test_normalize_antibiotic_collapses_whitespace_and_dash():
    assert normalize_antibiotic("  Trimethoprim\N{EN DASH}Sulfamethoxazole ") == (
        "trimethoprim-sulfamethoxazole"
    )
