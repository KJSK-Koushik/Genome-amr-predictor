# ADR 0003: Project-grouped frozen data splits

- Status: Accepted for Phase 2
- Date: 2026-08-18

## Decision

Assign each BioProject, and therefore every BioSample belonging to it, to exactly one of train,
validation, or test. Use a fixed seed and a deterministic balancing procedure across total genomes
and antibiotic/class strata. Freeze the resulting assignment before feature engineering.

## Reason

Random row splitting can place measurements from one isolate in several sets. Random isolate
splitting can still put samples collected and sequenced in the same project into several sets,
making evaluation overly optimistic. Keeping projects intact provides a practical leakage barrier
using metadata that is available for nearly every selected assembly.

## Consequences

- All antibiotic rows for a BioSample share one split.
- All BioSamples in a BioProject share one split.
- Very large projects prevent exact target fractions; achieved per-class counts are reported.
- Exact sequence duplicates are collapsed after FASTA verification.
- Before final model claims, a genomic-similarity sensitivity analysis should measure how results
  change under a stricter distance-based grouping rule.
