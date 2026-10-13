# ADR 0001: Species-specific first release

- Status: Accepted for Phase 1
- Date: 2026-08-17

## Decision

The first phenotype-prediction release will target *Escherichia coli*. ARG detection can later
support additional organisms, but phenotype models must remain species-specific unless evidence
supports a broader design.

## Reason

Resistance relationships, intrinsic mechanisms, mutation catalogs, genome structure, and
available AST labels vary by organism. A universal first model would create avoidable biological
confounding and make evaluation harder to trust.

## Consequences

- Phase 1 audits *E. coli* labels and candidate antibiotics.
- Later organisms require their own data audit and model approval.
- Unsupported organisms must not silently receive *E. coli* phenotype predictions.
