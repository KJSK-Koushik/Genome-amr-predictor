# Phase 1 NCBI AST Data Audit

Generated: 2026-08-18T04:42:03.399090+00:00

## Scope

- Species: *Escherichia coli*
- Query: `file import`
- Unique BioSamples after the configured species filter: 10,797
- BioSamples processed: 10,797
- AST rows parsed: 169,274
- Unique antibiotics observed: 93

This is an availability and label-quality audit. It does not train a model and does not establish
clinical performance.

## Binary-label usability

- Resistant or susceptible rows: 130,101
- Rows excluded from the binary endpoint: 39,173
- Duplicate isolate-antibiotic pairs: 52
- Conflicting binary-label pairs: 6

## Candidate antibiotic decision

| Antibiotic | R | S | Binary total | Minority fraction | Decision |
| --- | --- | --- | --- | --- | --- |
| ampicillin | 3838 | 4204 | 8042 | 0.477 | retain |
| ceftriaxone | 1580 | 5636 | 7216 | 0.219 | retain |
| ciprofloxacin | 1516 | 6541 | 8057 | 0.188 | retain |
| gentamicin | 1305 | 6428 | 7733 | 0.169 | retain |
| tetracycline | 3221 | 4154 | 7375 | 0.437 | retain |
| trimethoprim-sulfamethoxazole | 1484 | 6018 | 7502 | 0.198 | retain |
| cefepime | 746 | 3061 | 3807 | 0.196 | retain |
| meropenem | 217 | 6348 | 6565 | 0.033 | defer |

**Retain for Phase 2 cohort construction:** ampicillin, ceftriaxone, ciprofloxacin, gentamicin, tetracycline, trimethoprim-sulfamethoxazole, cefepime.

**Defer pending more or different data:** meropenem.

These decisions use only the configured count and class-balance thresholds. Phase 2 must still
verify assembly availability, duplicate genomes, lineage leakage, testing-standard consistency,
and genome quality before a final training cohort is frozen.

## Highest-coverage antibiotics

| Antibiotic | Unique samples | R | S | I | Excluded |
| --- | --- | --- | --- | --- | --- |
| ciprofloxacin | 8975 | 1516 | 6541 | 159 | 765 |
| ampicillin | 9635 | 3838 | 4204 | 16 | 1581 |
| gentamicin | 9902 | 1305 | 6428 | 103 | 2072 |
| trimethoprim-sulfamethoxazole | 9055 | 1484 | 6018 | 1 | 1557 |
| tetracycline | 8862 | 3221 | 4154 | 47 | 1443 |
| ceftriaxone | 7244 | 1580 | 5636 | 11 | 21 |
| meropenem | 7290 | 217 | 6348 | 14 | 712 |
| amoxicillin-clavulanic acid | 9039 | 961 | 5578 | 491 | 2012 |
| cefoxitin | 5593 | 1107 | 4242 | 130 | 114 |
| nalidixic acid | 4777 | 384 | 4389 | 4 | 0 |
| chloramphenicol | 6908 | 485 | 4138 | 70 | 2215 |
| azithromycin | 4626 | 55 | 4557 | 1 | 13 |
| streptomycin | 4150 | 1487 | 2653 | 0 | 10 |
| sulfisoxazole | 4129 | 1466 | 2615 | 4 | 44 |
| cefepime | 4665 | 746 | 3061 | 90 | 770 |

## Important limitations

- BioSample antibiograms are submitted by data providers and are not fully harmonized by NCBI.
- Recent or highly represented projects can dominate the records.
- This audit does not yet join BioSamples to GenBank assemblies.
- Counts do not guarantee that genomes will pass Phase 2 sequence-quality checks.
- `intermediate`, nonsusceptible, undefined, and other non-binary outcomes are excluded from the
  initial R/S endpoint.
- Laboratory method, platform, testing standard, source, time, and geography require cohort-level
  review before modeling.

## Reproducibility

The manifest records the source filename and SHA-256 fingerprint, configuration hash, retrieval time, and processed counts. The raw dataset remains outside version control.
