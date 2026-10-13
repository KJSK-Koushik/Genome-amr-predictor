# Genome AMR analysis: SAMN03075588

Run status: **completed**

Detection status: completed
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

## Phenotype predictions

| Antibiotic | Prediction | Resistance score | Status |
| --- | --- | --- | --- |
| ampicillin | R | 0.9877 | predicted |
| ceftriaxone | R | 0.9476 | predicted |
| ciprofloxacin | R | 0.9977 | predicted |
| gentamicin | S | 0.2139 | predicted |
| tetracycline | R | 0.9493 | predicted |
| cefepime | — | — | unsupported |
| trimethoprim/sulfamethoxazole | — | — | unsupported |

Scores are external model outputs, not locally calibrated confidence.

## Known determinant evidence

Total detected elements: 18; AMR elements: 17.

| Element | Type | Evidence | Class | Method |
| --- | --- | --- | --- | --- |
| gyrA_D87N | AMR | mutation | QUINOLONE | POINTX |
| gyrA_S83L | AMR | mutation | QUINOLONE | POINTX |
| ptsI_V25I | AMR | mutation | FOSFOMYCIN | POINTX |
| uhpT_E350Q | AMR | mutation | FOSFOMYCIN | POINTX |
| parC_E84V | AMR | mutation | QUINOLONE | POINTX |
| parC_S80I | AMR | mutation | QUINOLONE | POINTX |
| parE_I529L | AMR | mutation | QUINOLONE | POINTX |
| blaCTX-M-27 | AMR | gene | BETA-LACTAM | ALLELEX |
| mph(A) | AMR | gene | MACROLIDE | EXACTX |
| mrx(A) | AMR | gene | MACROLIDE | EXACTX |
| sul1 | AMR | gene | SULFONAMIDE | EXACTX |
| qacEdelta1 | STRESS | other_element | QUATERNARY AMMONIUM | ALLELEX |
| aadA5 | AMR | gene | AMINOGLYCOSIDE | EXACTX |
| dfrA17 | AMR | gene | TRIMETHOPRIM | EXACTX |
| tet(A) | AMR | gene | TETRACYCLINE | EXACTX |
| aph(6)-Id | AMR | gene | AMINOGLYCOSIDE | EXACTX |
| aph(3'')-Ib | AMR | gene | AMINOGLYCOSIDE | EXACTX |
| sul2 | AMR | gene | SULFONAMIDE | EXACTX |

## Errors

