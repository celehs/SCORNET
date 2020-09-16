## Overview

The Censor-Time Current Status Kernel Survival Estimator (CCSK) is a consistent, semi-supervised, non-parametric survival curve estimator optimized for efficient use of Electronic Health Record (EHR) data with a limited number of current status labels. Derived from van der Laan and Robins' Inverse Probability of Censoring Weighted (IPCW) estimator, it achieves locally efficient survival curve estimation using current status labels - binary indicators of phenotype status at censoring time - rather than more expensive event time labels. CCSK boosts efficiency over IPCW in the typical EHR setting by (1) utilizing unlabeled patients in a semi-supervised fashion, and (2) leveraging information-dense engineered EHR features to maximize imputation precision in the unlabeled set.

See Ahuja et al. (2020) for details.

## References

Ahuja Y, Liang L, Huang S, Cai T (2020). Censor-Time Current Status Kernel Estimator (CCSK): A Locally Efficient Survival Curve Estimator Using Electronic Health Record Data. BioArxiv.

van der Laan MJ and Robins JM (1998). Locally efficient estimation with current status data and time-dependent covariates. Journal of the American Statistical Association, 93, 693-701.