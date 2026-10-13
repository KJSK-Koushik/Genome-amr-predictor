# ADR 0002: Separate ARG detection from phenotype prediction

- Status: Accepted for Phase 1
- Date: 2026-08-17

## Decision

The product will expose two separate result types:

1. known ARG and supported mutation evidence from a curated detector; and
2. calibrated phenotype probabilities from one model per antibiotic.

## Reason

A gene hit and a laboratory phenotype are related but not equivalent. Treating them as one output
would overstate certainty and make model evaluation difficult to interpret.

## Consequences

- AMRFinderPlus integration is scheduled before phenotype-model development.
- Model features may use detector results, but reports retain the distinction.
- The interface includes research-use warnings and may abstain on unsupported inputs.
