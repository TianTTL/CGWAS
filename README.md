# C-GWAS

Combining GWAS summary statistics of multiple potentially related traits.

## Overview

C-GWAS is a powerful method for combining GWAS summary statistics of multiple potentially related traits and detect SNPs with multi-trait effects. 

C-GWAS begins with GWASs summary as inputs and outputs a single vector of combined p-values testing if the null is deviated. For each SNP, the null is the absence of any effect on all traits, and the alternative is that its effect deviates from 0 for at least one trait. C-GWAS integrates two different statistical methods with complementary statistical features to ensure the optimal power under various and complex scenarios while keeping a stable study-wide type-I error rate. The first method is called iterative effect based inverse covariance weighting (`i-EbICoW`) and the second method is called truncated Wald test (`TWT`). 

C-GWAS controls the study-wide type-I error rate in an empirical manner via simulations and adjust the resultant p-values in such a way that they are directly comparable with those from traditional GWAS of a single trait.

## System Requirements

**Depends**

R (>= 3.1.0)

**Imports**

data.table (>= 1.13.0)

foreach (>= 1.5.0)

MASS (>= 7.3-51)

## Download & Installation

First, we need to install R package `devtools`:

```R
install.packages('devtools')
library(devtools)
```

Then we just call

```R
install_github('https://github.com/FanLiuLab/CGWAS')
library(CGWAS)
```

Typically, this process takes between 20 seconds and 1 minute, depending on the network conditions.

## Example

example of the input GWAS file

```R
f1 <- read.table(file.path(system.file("extdata", package = "CGWAS"), 'Y1.assoc'), header=TRUE)
head(f1)
```

example of the input mean reference allele frequency file

```R
f2 <- read.table(file.path(system.file("extdata", package = "CGWAS"), 'MRAF'), header=TRUE)
head(f2)
```

example of the input SNP information file

```R
f3 <- read.table(file.path(system.file("extdata", package = "CGWAS"), 'SnpInfo'), header=TRUE)
head(f3)
```

demo of whole C-GWAS procedure implementing

```R
outputPath <- getwd() # the output files are in the current working directory
ExDataDir <- system.file("extdata", package = "CGWAS")
gwasFileName <- c("Y1.assoc", "Y2.assoc", "Y3.assoc",
                  "Y4.assoc", "Y5.assoc", "Y6.assoc",
                  "Y7.assoc", "Y8.assoc", "Y9.assoc",
                  "Y10.assoc", "Y11.assoc", "Y12.assoc")
gwasFilePath <- file.path(ExDataDir, gwasFileName)
snpFilePath <- file.path(ExDataDir, 'SnpInfo')
traitName <- c("Y1", "Y2", "Y3",
               "Y4", "Y5", "Y6",
               "Y7", "Y8", "Y9",
               "Y10", "Y11L", "Y12")
mrafFilePath <- file.path(ExDataDir, 'MRAF')
indSNPN = 1e5

cgwas(gwasFilePath, snpFilePath, outputPath,
      traitName = traitName, mrafFilePath = mrafFilePath, indSNPN = indSNPN)
```

This example takes 202.6 seconds on a laptop with 12 CPU cores (i7-9750H CPU @ 2.60 GHz). Since CGWAS uses half of the available cores for parallel computing by default, the number of parallel threads is 6.

## Instructions for Use

CGWAS implements whole procedure into a single function `cgwas`. 

The input files of `cgwas` contains GWAS summary statistics files,  SNP information file and mean reference allele frequency (MRAF) file. 

**GWAS summary statistics files**	These files should be in the space or tab delimited format. Each file contains two and only two columns. The first column is the regression betas of all SNPs and the second column is the P-values of all SNPs. These files require a header line with two items: BETA and P. Note that for each SNP, all betas must be based on the same reference allele. A reference to the GWAS summary statistics file is provided in the example.

**SNP information file**	A space delimited file consisting of three columns representing chromosome, base pair and SNP identifier. The file requires a header line with three items CHR, BP and SNP in the specified order. A reference to the SNP information file is provided in the example.

**MRAF file**	This file requires a header MRAF, and contains one column of mean frequency of the reference allele of each SNP. Note that the reference allele must correspond to the beta in the input GWAS. Reference allele frequencies are used to estimate the weights of input GWASs and all intermediate GWASs. This MRAF file is recommended but not obligatory with default `NULL`. If `NULL`, weights will still be estimated by approximations. A reference to the MRAF file is provided in the example.

The final and intermediate results of `cgwas` are saved in the output directory specified by user. Two new folders will be created in this directory: `Results/` and `Details/`.

The `Results/` folder contains final results

1. Manhattan plots `CGWAS-GWAS.jpg`, 
2. Q-Q plots `CGWASminpQQ.jpg` of C-GWAS and MinGWAS, 
3. study-wide suggestively significant SNPs table `SummarySugSigSNP.csv` 
4. p-values of all SNPs `C-GWAS.p`. 

The `Details/` folder contains all intermediate results including

1. A table of results of the `getI` function `SummaryGetI.txt`; 

2. A table of results of `getPsi` function `SummaryGetPsi.txt`, 
3. A table of summary results of `i-EbICoW` `EbICoW.txt`, 
4. A table of results of `getPsi` function applied to all EbICoW GWAS pairs `SummaryEbICoWGetPsi.txt`, 
5. A table of records of all `i-EbICoW` iterations `Summaryi-EbICoW.txt`, 
6. A table of C-GWAS LOESS model samples simulated in `getCoef` function `NullCGWAScorrection.txt`, 
7. A figure for the performance of `getCoef` based correction in C-GWAS simulation `NullCGWASdistribution.jpg` 
8. A table of MinGWAS LOESS model samples simulated in `getCoef` function `NullMinpcorrection.txt`, 
9. A figure for the performance of `getCoef` based correction in MinGWAS simulation `NullMinpdistribution.jpg`. 
10. If user choose to keep `i-EbICoW` output (`keepIEb = TRUE`), all summary statistics (effect sizes, *.beta and test statistics, *.stat) of EbICoW GWASs are saved in `Details/i-EbICoW`.

The detailed description of arguments and the legends of output tables and figures can be found in the manual of `cgwas` by calling:

```R
?CGWAS::cgwas
```

***Note:*** For normal user, setting the internal arguments to default values will run C-GWAS analysis correctly in most cases.

## Future Development Plan

Standing along functions facilitating the C-GWAS analysis will be implemented in a later version, e.g., `getI`, `getPsi`, `getPi`, `getCoef`, `i-EbICoW`, `TWT`.

