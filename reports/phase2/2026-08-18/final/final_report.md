# Phase 2 Sequence-Verified Cohort Report

Generated: 2026-10-12T18:24:28.340827+00:00

## Outcome

- Metadata-QC assemblies: 7,911
- Verified downloaded genomes: 7,911
- Final unique genomes: 7,911
- Final isolate-antibiotic rows: 47,485
- Exact duplicate sequence groups: 0
- Duplicate-sequence phenotype conflicts flagged: 0
- Downloaded sequence files reread: True
- Sequence quality thresholds applied: True

## Genome exclusions

- None

## Frozen splits after exact-sequence deduplication

- train: 5,592 BioSamples
- validation: 749 BioSamples
- test: 1,570 BioSamples

## AST metadata coverage

| Field | Value | Isolate-antibiotic rows |
|---|---|---:|
| laboratory_method | not_reported | 47,485 |
| laboratory_platform | Etest | 3 |
| laboratory_platform | KB panel | 6 |
| laboratory_platform | Microscan | 5 |
| laboratory_platform | Phoenix | 14,478 |
| laboratory_platform | Phoenix NMIC-203 card | 7 |
| laboratory_platform | Scan 500 | 13 |
| laboratory_platform | Sensititer | 1,071 |
| laboratory_platform | Sensititre | 24,975 |
| laboratory_platform | Sensititre microbroth dilution using EUVSEC and EUVSEC2 plates | 5 |
| laboratory_platform | Vitek | 4,818 |
| laboratory_platform | not_reported | 2,104 |
| testing_standard | CLSI | 44,804 |
| testing_standard | EUCAST | 2,561 |
| testing_standard | SFM | 96 |
| testing_standard | not_reported | 24 |

Missing method or testing-standard details remain explicit. These counts do not establish that
laboratory breakpoints and methods are comparable across all submitters; model development must
review these strata and the available measurement metadata.

The final cohort uses verified genome manifest entries. Exact duplicate sequences retain
one deterministic BioSample representative; conflicting antibiotic labels among duplicate sequences
are excluded and reported separately for scientific review. Sequence quality thresholds are applied
when configured. Closely related strains may remain across splits; genomic-distance validation is
still required before reporting model performance.
