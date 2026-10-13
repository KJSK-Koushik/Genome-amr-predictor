# Genome AMR analysis: mismatch-check

Run status: **partial**

Detection status: failed
Prediction status: completed

Empty CSV files from failed stages do not indicate absent genes or susceptibility.

## Limitations

- Research and educational use only; laboratory confirmation is required.
- Known determinant detection and predicted phenotypes are separate outputs.
- No detected determinant does not establish susceptibility.
- E. coli species is user-supplied; automatic identification is unavailable.
- External models; not validated on our cohort.
- Training overlap with our NCBI/PATRIC cohort is unknown.
- Scores are not locally calibrated probabilities.
- Cefepime and trimethoprim/sulfamethoxazole unsupported.
- Detection failed; gene presence/absence is unknown.

## Phenotype predictions

| Antibiotic | Prediction | Resistance score | Status |
| --- | --- | --- | --- |
| ampicillin | R | 0.9852 | predicted |
| ceftriaxone | S | 0.4908 | predicted |
| ciprofloxacin | S | 0.0418 | predicted |
| gentamicin | S | 0.2543 | predicted |
| tetracycline | S | 0.1103 | predicted |
| cefepime | — | — | unsupported |
| trimethoprim/sulfamethoxazole | — | — | unsupported |

Scores are external model outputs, not locally calibrated confidence.

## Errors

- detection: Detection result belongs to a different input genome
