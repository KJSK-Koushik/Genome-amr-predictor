# Genome AMR analysis: SAMN11638310

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
| ampicillin | R | 0.9956 | predicted |
| ceftriaxone | R | 0.7719 | predicted |
| ciprofloxacin | R | 0.9975 | predicted |
| gentamicin | R | 0.9616 | predicted |
| tetracycline | R | 0.8666 | predicted |
| cefepime | — | — | unsupported |
| trimethoprim/sulfamethoxazole | — | — | unsupported |

Scores are external model outputs, not locally calibrated confidence.

## Known determinant evidence

Total detected elements: 37; AMR elements: 35.

| Element | Type | Evidence | Class | Method |
| --- | --- | --- | --- | --- |
| gyrA_D87N | AMR | mutation | QUINOLONE | POINTX |
| gyrA_S83L | AMR | mutation | QUINOLONE | POINTX |
| ftsI_I336IKYRI | AMR | mutation | BETA-LACTAM | POINTX |
| parC_S80I | AMR | mutation | QUINOLONE | POINTX |
| parE_S458A | AMR | mutation | QUINOLONE | POINTX |
| ompF_Y253EfsTer3 | AMR | mutation | BETA-LACTAM | POINTX |
| cyaA_S352T | AMR | mutation | FOSMIDOMYCIN | POINTX |
| blaTEM-1 | AMR | gene | BETA-LACTAM | ALLELEX |
| blaNDM-5 | AMR | gene | BETA-LACTAM | ALLELEX |
| ble | AMR | gene | BLEOMYCIN | EXACTX |
| blaCTX-M-15 | AMR | gene | BETA-LACTAM | ALLELEX |
| tet(B) | AMR | gene | TETRACYCLINE | EXACTX |
| aadA2 | AMR | gene | AMINOGLYCOSIDE | EXACTX |
| dfrA12 | AMR | gene | TRIMETHOPRIM | EXACTX |
| catB3 | AMR | gene | PHENICOL | PARTIAL_CONTIG_ENDX |
| blaOXA-1 | AMR | gene | BETA-LACTAM | ALLELEX |
| aac(6')-Ib-cr5 | AMR | gene | AMINOGLYCOSIDE/QUINOLONE | ALLELEX |
| aac(3)-IId | AMR | gene | AMINOGLYCOSIDE | EXACTX |
| aadA5 | AMR | gene | AMINOGLYCOSIDE | EXACTX |
| dfrA17 | AMR | gene | TRIMETHOPRIM | EXACTX |
| sul1 | AMR | gene | SULFONAMIDE | EXACTX |
| qacEdelta1 | STRESS | other_element | QUATERNARY AMMONIUM | ALLELEX |
| blaNDM-5 | AMR | gene | BETA-LACTAM | ALLELEX |
| ble | AMR | gene | BLEOMYCIN | EXACTX |
| blaTEM-1 | AMR | gene | BETA-LACTAM | ALLELEX |
| blaCTX-M-15 | AMR | gene | BETA-LACTAM | ALLELEX |
| tet(B) | AMR | gene | TETRACYCLINE | EXACTX |
| aac(6')-Ib-cr5 | AMR | gene | AMINOGLYCOSIDE/QUINOLONE | ALLELEX |
| blaOXA-1 | AMR | gene | BETA-LACTAM | ALLELEX |
| catB3 | AMR | gene | PHENICOL | PARTIALX |
| dfrA12 | AMR | gene | TRIMETHOPRIM | EXACTX |
| aadA2 | AMR | gene | AMINOGLYCOSIDE | EXACTX |
| dfrA17 | AMR | gene | TRIMETHOPRIM | EXACTX |
| aadA5 | AMR | gene | AMINOGLYCOSIDE | EXACTX |
| aac(3)-IId | AMR | gene | AMINOGLYCOSIDE | EXACTX |
| qacEdelta1 | STRESS | other_element | QUATERNARY AMMONIUM | ALLELEX |
| sul1 | AMR | gene | SULFONAMIDE | EXACTX |

## Errors

