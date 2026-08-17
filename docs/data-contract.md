# Phase 1 Data Contract

## Unit of observation

One normalized row represents one BioSample isolate tested against one antibiotic. A single
isolate can therefore contribute several rows.

## Normalized fields

| Field | Required | Description |
|---|---:|---|
| `biosample_accession` | yes | Stable BioSample accession, such as `SAMN...` |
| `organism` | yes | Organism name supplied by BioSample |
| `taxonomy_id` | no | NCBI taxonomy identifier |
| `bioproject_accession` | no | Linked BioProject accession when present |
| `sra_accession` | no | Linked SRA accession when present |
| `isolate` | no | Submitter-provided isolate identifier |
| `collection_date` | no | Submitter-provided collection date |
| `geo_loc_name` | no | Submitter-provided geographic location |
| `isolation_source` | no | Submitter-provided isolation source |
| `host` | no | Submitter-provided host |
| `antibiotic` | yes | Lowercase normalized antibiotic name |
| `phenotype_raw` | yes | Original submitted phenotype string |
| `phenotype_normalized` | yes | `resistant`, `susceptible`, `intermediate`, or `excluded` |
| `binary_label` | no | `1` for resistant, `0` for susceptible, empty otherwise |
| `measurement_sign` | no | Comparison operator associated with the measurement |
| `measurement` | no | MIC or disk-diffusion value as submitted |
| `measurement_units` | no | Measurement units as submitted |
| `laboratory_method` | no | Laboratory typing method |
| `laboratory_platform` | no | Platform used for testing |
| `vendor` | no | Test vendor |
| `method_version_or_reagent` | no | Method version or reagent |
| `testing_standard` | no | Testing standard, such as CLSI |

## Version 1 label policy

- `resistant` and `R` map to binary label `1`.
- `susceptible`, `sensitive`, and `S` map to binary label `0`.
- `intermediate`, `I`, nonsusceptible, susceptible-dose dependent, undefined, and missing
  outcomes do not receive binary labels.
- Conflicting duplicate labels are flagged for later removal or review.
- Synthetic phenotype labels are prohibited.

## Provenance

Each audit writes a manifest containing the query, retrieval timestamp, reported NCBI result
count, processed record count, source URLs, configuration, and SHA-256 hash of that configuration.
