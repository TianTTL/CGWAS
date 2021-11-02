#' combine GWAS (C-GWAS)
#'
#' A whole solution for illustrating multi-trait effect of all SNPs with a set
#' of GWASs summary data.
#'
#' C-GWAS begins with GWASs summary as inputs and outputs a single vector of
#' combined p-values testing if the null is deviated. For each SNP, the null is
#' the absence of any effect on all traits, and the alternative is that its
#' effect deviates from 0 for at least one trait.
#'
#' To estimate the effect of each SNP more accurately, several data organization
#' standard should be noticed. (1) The effect size in each GWAS should be in the
#' same scale. (2) In one GWAS, the minimum samples size of each SNP should
#' be no less than 60% of the median.
#'
#' C-GWAS integrates two different statistical methods to ensure the optimal
#' power under various and complex scenarios while keeping a stable study-wide
#' type-I error rate. The first method uses an iterative effect based inversed
#' covariance weighting (i-EbICoW), which appears the most powerful when the
#' assumption (all SNPs share the same variance–covariance matrix of effect
#' sizes across traits) tends to be satisfied or moderately violated. The second
#' method is a truncated Wald test (TWT) which is more powerful than i-EbICoW
#' when the assumption tends to be severely violated. For each SNP, TWT proposes
#' the best subset of phenotypes by applying the Wald test to all subsets
#' determined under a series of preset thresholds.
#'
#' C-GWAS controls the study-wide type-I error rate in an empirical manner via
#' simulations and adjust the resultant p-values in such a way that they are
#' directly comparable with those from traditional GWAS of a single trait.
#'
#' @param gwasFilePath A string list of path to GWASs summary data files.
#' @param assocColInx A number list indicates column number of SNP index (CHR &
#' BP & SNP) and association result (BETA & P).
#' @param mafFilePath Path to allele frequency file.
#' @param threadN The number of threads for parallel computing. The default value
#' is half of the maximum number of cores on the hardware.
#' @param outputPath Result output directory.
#' @param traitName A list of trait names.
#' @param indSNPN Independent SNPs number.
#' @param simulTime Number of simulations for multiple testing correction. This
#' parameter greatly affects the operating efficiency.
#' @param ebicoPwrInc Minimum EbICo power increase ratio.
#' @param sampleNInc Equivalent sample size in EbICo increase ratio.
#'
#' @return Several .csv files and figures will be generated in \code{outputPath}.
#'
#' \code{MutipleTestingCorrection.txt}.
#'
#' \code{GWAShits.csv}.
#'
#' \code{EbICohits.csv}.
#'
#' \code{SWaThits.csv}.
#'
#' \code{Summaryhits.csv}.
#'
#' \code{Summary.csv}.
#'
#' @examples
#'
#' gwasFileDir <- system.file("extdata", package = "CGWAS")
#' gwasFileName <- c("phe1.assoc", "phe2.assoc", "phe3.assoc",
#'                   "phe4.assoc", "phe5.assoc")
#' gwasFilePath <- file.path(gwasFileDir, gwasFileName)
#' assocColInx <- c(1, 2, 3, 4, 5)
#' outputPath <- getwd()
#' traitName <- c("phe1", "phe2", "phe3", "phe4", "phe5")
#' simulTime <- 200
#'
#' cgwas(gwasFilePath, assocColInx, outputPath,
#'       traitName = traitName, simulTime = simulTime)
#'
#' @export
cgwas <- function(gwasFilePath,
                  assocColInx,
                  outputPath,
                  mafFilePath = NULL,
                  threadN = ceiling(parallel::detectCores() / 2),
                  traitName = NULL,
                  indSNPN = 1e6,
                  simulTime = NULL,
                  ebicoPwrInc = 1.1,
                  sampleNInc = 0.7) {

  print("---- C-GWAS version 1.8.2 start ----")

  cgwasenv <- list()

  # A string list of path to GWASs summary data files.
  # export
  cgwasenv$.GWAS_FILE_PATH <- gwasFilePath

  # Column number of SNP index (CHR & BP & SNP) and association data (BETA & P)
  # export
  cgwasenv$.ASSOC_COLUMN_INDEX <- assocColInx

  # Path to allele frequency file.
  # export
  cgwasenv$.MAF_FILE_PATH <- mafFilePath

  # whether MAF file exist.
  # hide
  cgwasenv$.MAF_FILE_EXIST <- !is.null(mafFilePath)

  # The number of threads for parallel computing.
  # export
  chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")

  if (nzchar(chk) && chk == "TRUE") {
    # use 2 cores in CRAN/Travis/AppVeyor
    cgwasenv$.PARAL_NUM <- min(parallel::detectCores() - 1, 2)
  } else {
    cgwasenv$.PARAL_NUM <- min(parallel::detectCores() - 1, threadN)
  }


  # Result output directory.
  # export
  cgwasenv$.CGWAS_DIR <- outputPath

  # A list of trait names.
  # export
  if (is.null(traitName)) {
    cgwasenv$.TRAIT_NAME <- basename(gwasFilePath)
  } else {
    cgwasenv$.TRAIT_NAME <- traitName
  }

  # Trait numbers.
  # hide
  cgwasenv$.TRAIT_NUM <- length(cgwasenv$.TRAIT_NAME)

  # SNP numbers.
  # hide
  cgwasenv$.SNP_N <- length(readLines(cgwasenv$.GWAS_FILE_PATH[1])) - 1

  # Independent SNP number.
  # export
  cgwasenv$.IND_SNP_N <- indSNPN

  # FDR test set.
  # hide
  cgwasenv$.FDR_SET <- c(seq(0.0001, 0.0009, 0.0001), seq(0.001, 0.05, 0.001))

  # main effect threshold
  # hide
  cgwasenv$.P_THRD <- 1 / cgwasenv$.IND_SNP_N

  # minimum EbICo power increase ratio.
  # export
  cgwasenv$.MIN_EbiCo_POWER_INC <- ebicoPwrInc

  # minimum correlation difference.
  # hide
  cgwasenv$.MIN_CORR_DIFF <- 0.05

  # high correlation restriction (squared).
  # hide
  cgwasenv$.HIGH_CORR2_RES <- 0.7

  # high correlation restriction.
  # hide
  cgwasenv$.HIGH_CORR_RES <- sqrt(cgwasenv$.HIGH_CORR2_RES)

  # Equivalent sample size increase ratio.
  # export
  cgwasenv$.SAMPLE_SIZE_INC <- sampleNInc

  # SWaT stratification cutoff.
  # hide
  cgwasenv$.SWaT_STRAT_CUT <- 10^(seq(0,-5,-1/3))[-1]

  # Simulated SNP number.
  # hide
  cgwasenv$.SIMUL_SNP_N <- 5e4

  # Number of simulations for multiple testing correction.
  # export
  if (is.null(simulTime)) {
    cgwasenv$.SIMUL_N <- 100 * cgwasenv$.IND_SNP_N / cgwasenv$.SIMUL_SNP_N
  } else {
    cgwasenv$.SIMUL_N <- simulTime
  }

  # LOESS span vector.
  # hide
  cgwasenv$.LOESS_SPAN_V <- c(0.03, 0.1, 0.75)

  # LOESS interval p.
  # hide
  cgwasenv$.LOESS_INTER_P <- c(0.05, 0.001)

  # Loci interval.
  # hide
  cgwasenv$.LOCI_INTER <- 2.5e5

  # create directory of intermediate result
  cgwasenv$.CGWAS_TEMPDATA_PATH <- file.path(cgwasenv$.CGWAS_DIR, "Tempdata")
  cgwasenv$.CGWAS_COLDATA_PATH <- file.path(cgwasenv$.CGWAS_DIR, "Tempdata", "Coldata")
  cgwasenv$.CGWAS_RESULT_PATH <- file.path(cgwasenv$.CGWAS_DIR, "Result")
  dir.create(cgwasenv$.CGWAS_TEMPDATA_PATH, recursive = T, showWarnings = F)
  dir.create(cgwasenv$.CGWAS_COLDATA_PATH, recursive = T, showWarnings = F)
  dir.create(cgwasenv$.CGWAS_RESULT_PATH, recursive = T, showWarnings = F)

  print(cgwasenv)

  # Step1 Statistic Summary Extraction
  ts <- Sys.time()
  step1(cgwasenv)
  print("Step1 finished")
  print(difftime(Sys.time(), ts))

  # Step2 Inflation Correlation Estimation
  ts <- Sys.time()
  step2(cgwasenv)
  print("Step2 finished")
  print(difftime(Sys.time(), ts))

  # Step3 EbICo
  ts <- Sys.time()
  step3(cgwasenv)
  print("Step3 finished")
  print(difftime(Sys.time(), ts))

  # Step4 SWaT
  ts <- Sys.time()
  step4(cgwasenv)
  print("Step4 finished")
  print(difftime(Sys.time(), ts))

  # Step5 Summary
  ts <- Sys.time()
  step5(cgwasenv)
  print("Step5 finished")
  print(difftime(Sys.time(), ts))
}
