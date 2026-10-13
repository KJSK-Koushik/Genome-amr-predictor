# Phase 2 Assembly-Linked Cohort Report

Generated: 2026-08-18T06:21:44.570654+00:00

## Outcome

- Species: *Escherichia coli*
- Normalized binary AST labels: 49,702
- BioSamples with usable labels: 8,273
- BioSamples with a selected assembly: 7,931
- BioSamples passing metadata QC: 7,911
- Frozen isolate-antibiotic rows: 47,485

## Label cleaning

- Duplicate label rows removed: 29
- Conflicting isolate-antibiotic pairs removed: 1
- Binary rows without BioSample accessions removed: 0

## Metadata-QC exclusions

- contig_n50_too_low: 2
- no_assembly: 342
- too_many_contigs: 20

## Frozen group-aware split

- train: 5,592 BioSamples (70.7%; target 70.0%)
- validation: 749 BioSamples (9.5%; target 15.0%)
- test: 1,570 BioSamples (19.8%; target 15.0%)

The cohort contains 46 indivisible BioProject/BioSample groups. The five
largest are:

- PRJNA292663: 2,638 BioSamples
- PRJNA809394: 2,075 BioSamples
- PRJNA514245: 1,545 BioSamples
- PRJNA278886: 508 BioSamples
- PRJNA885502: 335 BioSamples

All rows from the same BioProject are assigned to one split. If BioProject is unavailable, the
BioSample accession is used as the group. A deterministic allocator uses the versioned seed while
balancing sample counts and antibiotic/class strata as closely as whole-project grouping permits.

## Important boundary

This report freezes the assembly-linked, metadata-QC cohort. The genome download stage verifies
FASTA integrity and exact sequence hashes. No resistance detector or phenotype model is trained
in Phase 2.

## Reproducibility

See `data_manifest.json` for source hashes, NCBI endpoint, retrieval time, configuration hash,
software version, and cohort counts. The raw AST and genome files remain outside version control.
