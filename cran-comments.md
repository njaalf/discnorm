## Test environments
* local R installation, R 3.6.1
* x86_64-apple-darwin15.6.0 (64-bit)
* macOS High Sierra 10.13.6

## R CMD check results

0 errors | 0 warnings | 0 notes

## 0.1.0 Submission

Received feedback from CRAN submission officer Jelena Saf on 12th may
to fix cat(" B = ") to only show with verbose = TRUE. This has been done. 

## 0.1.1 Submission

In computeT, change option use_pbv from 0 to 2 in sirt::polychoric2
In bootTest, change the polychoric estimator from sirt::polychoric2 to lavaan::lavCor