# Phase 3 Completion Note

Completed on 2026-10-13 (Asia/Calcutta).

## Outcome

The package now runs AMRFinderPlus on plain or compressed assembled *E. coli* FASTA and imports
existing detector TSVs. It writes structured gene, allele, and mutation evidence with native
coordinates, match methods, identity/coverage metrics, curated class/subclass labels, hierarchy
nodes, warnings, and the original rows. Every live run records executable/database versions and
hashes, input/sequence hashes, adapter code hash, command, and output fingerprints.

Runs use explicit database release 2026-08-07.1 with AMRFinderPlus 4.2.7 through Ubuntu WSL.
BLAST+ 2.12.0 and HMMER 3.4 are installed. The default analysis uses nucleotide mode,
`--organism Escherichia`, and native curated thresholds. No database update occurs during detection.

`detect` provides the executable adapter; `parse-detection` provides offline normalization.
Header-only no-hit output remains distinct from an empty or malformed output or a failed command.
Existing output directories containing files are rejected to prevent stale results being reused.

## Validation

All **66 package tests passed** and **Ruff passed**. Tests include current and legacy TSV headers,
allele and mutation calls, HMM missing values, partial/internal-stop warnings, reverse-strand
coordinates, truncated rows, invalid metrics, duplicate FASTA identifiers, gzip input, command
failure/timeout, no-hit interpretation, and NCBI's unchanged release expected-output fixture.

Three real genomes from the frozen Phase 2 cohort completed live detection:

| BioSample | Assembly | Frozen split | Total hits | AMR hits | Mutation hits |
|---|---|---|---:|---:|---:|
| SAMN03075588 | GCF_000770275.1 | train | 18 | 17 | 7 |
| SAMN03177615 | GCF_000819375.1 | validation | 6 | 6 | 0 |
| SAMN11638310 | GCA_014625455.1 | test | 37 | 35 | 7 |

An independent verifier passed raw-to-structured equality, summary counts, input/output hashes,
executable/database fingerprints, and version provenance. These are integration checks on three
genomes, not an assessment of population detection accuracy or phenotype prediction.

All **four NCBI installation reference tests passed** in the native Linux directory with one
thread: protein-only, nucleotide-only, combined protein/nucleotide, and disruption calls. The
upstream test script was changed only from six threads to one; expected outputs were unchanged.
The earlier six-thread attempt encountered an out-of-memory kill of BLASTX in WSL. A one-thread retry
exposed an upstream StxTyper path-quoting issue in the original installation directory. The active
runtime was moved to `/opt/genome-amr/amrfinderplus/4.2.7` to avoid spaces and shell metacharacters.
The earlier failures are retained in local logs.

## Evidence and reproduction

- [Detection contract and CLI example](phase-3-detection-contract.md)
- [Live smoke summary](../reports/phase3/2026-10-13/smoke-linux-runtime/smoke_summary.json)
- [Independent verification](../reports/phase3/2026-10-13/smoke-linux-runtime/verification.json)
- [Train genome evidence](../reports/phase3/2026-10-13/smoke-linux-runtime/SAMN03075588/detections.json)
- [Validation genome evidence](../reports/phase3/2026-10-13/smoke-linux-runtime/SAMN03177615/detections.json)
- [Test genome evidence](../reports/phase3/2026-10-13/smoke-linux-runtime/SAMN11638310/detections.json)
- [Official fixture provenance](../tests/fixtures/amrfinder-provenance.md)
- [NCBI reference-suite log](../reports/phase3/2026-10-13/ncbi-installation-test.log)
- [Reference-suite verification](../reports/phase3/2026-10-13/ncbi-installation-verification.json)

```powershell
.venv/Scripts/python.exe scripts/run_phase3_smoke.py --output-dir reports/phase3/new-smoke
.venv/Scripts/python.exe scripts/verify_phase3.py --report-dir reports/phase3/new-smoke
```

The Linux release archive SHA-256 is
`68045a8bccdbe3c5dcdf941bebe2352ed419758a9914c41f48f0bbbd6fbade56`.
Genomes, executable/database downloads, and temporary decompressed sequence files remain outside
version control. Saved report files contain detection evidence and compact provenance.

## Scientific boundary and next phase

Detection evidence remains separate from phenotype prediction. Neither a gene hit nor a no-hit
result establishes a laboratory R/S phenotype. Partial and internal-stop calls may not represent
functional genes; method and scope remain visible for downstream review.

Ontology normalization preserves the curated NCBI classes and hierarchy nodes. A separate
mechanism annotation is not supplied by the native TSV, so the mechanism field remains explicitly
unreported. Mechanism enrichment requires an additional sourced, versioned mapping. No broad
class is converted into a drug-specific phenotype assertion. Species is assumed from the
supported cohort/input scope, not inferred by this adapter.

Phase 4 covers cohort-wide detector features, model training, calibration, evaluation, and model
cards. Before final performance claims, genomic-similarity sensitivity analysis and laboratory
metadata limitations documented in Phase 2 must be addressed.
