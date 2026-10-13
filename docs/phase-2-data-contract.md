# Phase 2 data contract

Phase 2 turns the NCBI AST export into a reproducible, sequence-linked *Escherichia coli*
cohort. The outputs are research data, not clinical advice.

## Inputs

- NCBI Pathogen Detection AST Browser CSV, filtered to *E. coli*.
- NCBI Datasets assembly reports resolved from each BioSample accession.
- NCBI genome FASTA files downloaded by assembly accession.
- Versioned settings in `configs/phase2_cohort.json`.

Raw AST and FASTA files remain outside Git. Reports record input fingerprints, configuration,
software version, retrieval time, and source endpoints so the cohort can be rebuilt.

## Label rules

- Only the seven configured antibiotics are retained.
- Resistant observations map to binary label `1`; susceptible maps to `0`. Intermediate and
  ambiguous outcomes remain available for audit but are excluded from this binary endpoint.
- Repeated identical BioSample-antibiotic labels are collapsed.
- A BioSample-antibiotic pair containing both binary outcomes is excluded as a conflict.
- Synthetic phenotype labels are prohibited. Synthetic data may later test software, but may not
  be mixed into training or evaluation.

The 2026-08-18 source export has laboratory platform and testing-standard columns, but no
separate laboratory-method or standard-year column. The final report audits missing metadata
and the observed standard/platform strata. Reported R/S labels remain the ground truth; no
breakpoint harmonization or reclassification from MIC is claimed. The frozen table's measurement
units are unpopulated for this export; the original MIC (mg/L) and disk-diffusion (mm) columns must
be consulted before any measurement-based analysis.

## Assembly selection and quality control

One assembly is selected deterministically per BioSample, preferring current RefSeq and higher
assembly quality. The current boundaries are 4.0--6.5 Mb total length, at most 500 contigs, and
contig N50 of at least 20 kb. Reported ANI species mismatches are rejected. CheckM completeness
and contamination are applied when NCBI reports them; their absence is retained as a warning
rather than treated as failure because coverage is incomplete.

Every retained FASTA must decompress, contain at least one record, use the accepted IUPAC DNA
alphabet, and produce consistent length, contig, N50, GC, ambiguity, file SHA-256, and normalized
sequence SHA-256 measurements.

## Split and leakage controls

The split is frozen with the seed `genome-amr-phase2-v1`. All BioSamples from the same BioProject
stay in one split. Assignment balances the overall number of genomes and the available
antibiotic/class strata while preserving whole projects. Large BioProjects make exact 70/15/15
fractions impossible, so the achieved counts are reported rather than hidden.

After download, exact normalized genome-sequence duplicates are collapsed to one deterministic
representative. Any phenotype disagreement inside an exact-duplicate group is written to a
separate conflict report and excluded even on the selected representative. This controls exact
duplication and project leakage; it does not claim
that all closely related strains have been detected. Stronger genomic-distance clustering is a
documented refinement for model evaluation.

## Principal outputs

- `cohort.csv`: one metadata-QC-passing BioSample-antibiotic label per row.
- `assemblies.csv`: one selected, metadata-QC-passing assembly per BioSample.
- `split_summary.csv`: antibiotic and class counts for each frozen split.
- `data_manifest.json`: input/configuration provenance and SHA-256 fingerprints.
- `final/final_cohort.csv`: sequence-verified, exact-deduplicated training rows.
- `final/final_assemblies.csv`: sequence and FASTA metrics for retained genomes.
- `final/genome_exclusions.csv`: explicit post-download exclusion reasons.
- `final/duplicate_sequence_label_conflicts.csv`: disagreements requiring scientific review.

The portable `fasta_path` is relative to the local genome download directory. A consumer must
join on `assembly_accession`; row order is never an identifier.

## Resuming the frozen cohort

Use the project virtual environment. To resume downloads without rebuilding the frozen cohort:

```powershell
./scripts/run_phase2_cohort.ps1 -ResumeExistingCohort -OutputDate 2026-08-18
```

Finalization requires `--genome-dir` and `--config`. It rereads each downloaded FASTA, checks
both stored hashes, applies the configured length/contig/N50 thresholds to measured sequence
metrics, and writes final per-antibiotic split counts. Damaged gzip files are preserved with an
`.invalid-*` suffix and downloaded again once. Native command failures stop the phase wrapper.

Datasets rehydration is the default sequence transport. The optional `download-genomes
--prefer-ftp` flag first tries compressed FASTA downloads over HTTPS from the official NCBI
Genomes FTP archive, using the accession and filename from the Datasets fetch manifest. Files
are validated in a temporary location before installation; Datasets supplies any remaining files.
The download manifest records the requested transport. Resumed caches may use either transport;
the sequence hashes use the same normalization for both.
