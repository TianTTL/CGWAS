# C-GWAS
A whole solution for illustrating multi-trait effect of all SNPs with a set of GWASs summary data.

## Overview
C-GWAS begins with GWASs summary as inputs and outputs a single vector of combined p-values testing if the null is deviated. For each SNP, the null is the absence of any effect on all traits, and the alternative is that its
effect deviates from 0 for at least one trait.

C-GWAS integrates two different statistical methods to ensure the optimal power under various and complex scenarios while keeping a stable study-wide type-I error rate. The first method uses an iterative effect based inversed
covariance weighting (i-EbICoW), which appears the most powerful when the assumption (all SNPs share the same variance–covariance matrix of effect sizes across traits) tends to be satisfied or moderately violated. The second
method is a truncated Wald test (TWT) which is more powerful than i-EbICoW when the assumption tends to be severely violated. For each SNP, TWT proposes the best subset of phenotypes by applying the Wald test to all subsets
determined under a series of preset thresholds.

C-GWAS controls the study-wide type-I error rate in an empirical manner via simulations and adjust the resultant p-values in such a way that they are directly comparable with those from traditional GWAS of a single trait.

## Download & Install
```
install.packages('devtools')
library(devtools)
install_github('https://github.com/TianTTL/CGWAS')
```

