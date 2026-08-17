# Phase 1 NCBI AST Data Audit

Generated: 2026-08-17T17:09:43.203592+00:00

## Scope

- Species: *Escherichia coli*
- Query: `Escherichia coli[Organism] AND antibiogram[filter]`
- NCBI matching BioSamples at retrieval: 11,075
- BioSamples processed: 11,075
- AST rows parsed: 171,471
- Unique antibiotics observed: 96

This is an availability and label-quality audit. It does not train a model and does not establish
clinical performance.

## Binary-label usability

- Resistant or susceptible rows: 131,831
- Rows excluded from the binary endpoint: 39,640
- Duplicate isolate-antibiotic pairs: 52
- Conflicting binary-label pairs: 6

## Candidate antibiotic decision

| Antibiotic | R | S | Binary total | Minority fraction | Decision |
| --- | --- | --- | --- | --- | --- |
| ampicillin | 3917 | 4267 | 8184 | 0.479 | retain |
| ceftriaxone | 1646 | 5665 | 7311 | 0.225 | retain |
| ciprofloxacin | 1546 | 6587 | 8133 | 0.190 | retain |
| gentamicin | 1321 | 6467 | 7788 | 0.170 | retain |
| tetracycline | 3250 | 4169 | 7419 | 0.438 | retain |
| trimethoprim-sulfamethoxazole | 1500 | 6053 | 7553 | 0.199 | retain |
| cefepime | 768 | 3084 | 3852 | 0.199 | retain |
| meropenem | 225 | 6400 | 6625 | 0.034 | defer |

**Retain for Phase 2 cohort construction:** ampicillin, ceftriaxone, ciprofloxacin, gentamicin, tetracycline, trimethoprim-sulfamethoxazole, cefepime.

**Defer pending more or different data:** meropenem.

These decisions use only the configured count and class-balance thresholds. Phase 2 must still
verify assembly availability, duplicate genomes, lineage leakage, testing-standard consistency,
and genome quality before a final training cohort is frozen.

## Highest-coverage antibiotics

| Antibiotic | Unique samples | R | S | I | Excluded |
| --- | --- | --- | --- | --- | --- |
| ampicillin | 9785 | 3917 | 4267 | 16 | 1589 |
| ciprofloxacin | 9083 | 1546 | 6587 | 161 | 795 |
| gentamicin | 9991 | 1321 | 6467 | 108 | 2101 |
| trimethoprim-sulfamethoxazole | 9113 | 1500 | 6053 | 1 | 1564 |
| tetracycline | 8910 | 3250 | 4169 | 47 | 1447 |
| ceftriaxone | 7344 | 1646 | 5665 | 12 | 25 |
| amoxicillin-clavulanic acid | 9187 | 976 | 5685 | 499 | 2030 |
| meropenem | 7385 | 225 | 6400 | 18 | 743 |
| cefoxitin | 5649 | 1118 | 4283 | 133 | 115 |
| nalidixic acid | 4811 | 390 | 4416 | 5 | 0 |
| chloramphenicol | 6964 | 488 | 4162 | 71 | 2243 |
| azithromycin | 4648 | 55 | 4579 | 1 | 13 |
| streptomycin | 4171 | 1499 | 2660 | 0 | 12 |
| sulfisoxazole | 4147 | 1471 | 2628 | 4 | 44 |
| cefepime | 4740 | 768 | 3084 | 91 | 799 |

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

See `data_manifest.json` for source URLs, query translation, configuration hash, retrieval time,
and processed counts. Raw XML was parsed in batches and was not retained.
