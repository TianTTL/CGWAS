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
#' @expose
cgwas <- function(gwasFilePath,
                  assocColInx,
                  traitName = basename(gwasFilePath[1]),
                  exNa = TRUE,
                  mafFilePath = NULL,
                  outputPath,
                  keepIEb = FALSE,
                  threadN = ceiling(parallel::detectCores()/2),
                  indSNPN = 1e6,
                  simulDep = 100,
                  pStudy = 0.05/indSNPN,
                  pSugst = 1/indSNPN,
                  pMain = 3/indSNPN,
                  ebicoPwrInc = 1,
                  minCorDiff = 0.05,
                  hiCorRestc = 0.5,
                  sampleNInc = 0.5,
                  twtStCo = 10^(seq(0,1-ceiling(log10(indSNPN)),-1/3))[-1],
                  LoInrP = c(0.05,0.001),
                  LoSpanV = c(0.03,0.1,0.75),
                  lociInr = 2.5e5) {

  cat("
   ______       ______ _       __ ___    _____
  / ____/      / ____/| |     / //   |  / ___/
 / /   ______ / / __  | | /| / // /| |  \\__ \\
/ /___/_____// /_/ /  | |/ |/ // ___ | ___/ /
\\____/       \\____/   |__/|__//_/  |_|/____/      version 1.9.1")

  cgwasenv <- list()

  # A string list of path to GWASs summary data files.
  # expose
  cgwasenv$.GWAS_FILE_PATH <- gwasFilePath

  # Column number of SNP index (CHR & BP & SNP) and association data (BETA & P)
  # expose
  cgwasenv$.ASSOC_COLUMN_INDEX <- assocColInx

  # A list of trait names.
  # expose
  cgwasenv$.TRAIT_NAME <- traitName

  # Trait numbers.
  # hide
  cgwasenv$.TRAIT_NUM <- length(cgwasenv$.TRAIT_NAME)

  # SNP numbers.
  # hide
  cgwasenv$.SNP_N <- length(readLines(cgwasenv$.GWAS_FILE_PATH[1])) - 1

  # Exclude NA (exclude SNP with NA in at least one GWAS, otherwise replace NA with BETA=0 and P=1)
  # expose
  cgwasenv$.EXCLUDE_NA = exNa

  # Path to allele frequency file.
  # expose
  cgwasenv$.MAF_FILE_PATH <- mafFilePath

  # whether MAF file exist.
  # hide
  cgwasenv$.MAF_FILE_EXIST <- !is.null(mafFilePath)

  # Result output directory.
  # expose
  cgwasenv$.CGWAS_DIR <- outputPath

  # Keep i-EbICoW output (keep BETA and STAT of all i-EbICoW combinations or drop them)
  # expose
  cgwasenv$.KEEP_EbICoW = keepIEb

  # The number of threads for parallel computing.
  # expose
  chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
  if (nzchar(chk) && chk == "TRUE") {
    # use 2 cores in CRAN/Travis/AppVeyor
    cgwasenv$.PARAL_NUM <- min(parallel::detectCores() - 1, 2)
  } else {
    cgwasenv$.PARAL_NUM <- min(parallel::detectCores() - 1, threadN)
  }

  # Independent SNP number.
  # integer, >=1e4, multiple of 500, default 1e6
  # expose
  if (indSNPN < 1e4 | indSNPN %% 500 != 0) {
    return(message("ERROR: parament indSNP invalid!"))
  } else {
    cgwasenv$.IND_SNP_N <- indSNPN
  }

  # Quantile simulation depth
  # integer, >=100, multiple of 100, default 100
  # expose
  if (simulDep < 100 | simulDep %% 100 != 0) {
    return(message("ERROR: parament simulDep invalid!"))
  } else {
    cgwasenv$.SIMUL_DEP <- simulDep
  }

  # Simulated SNP number.
  # hide
  cgwasenv$.SIMUL_SNP_N <- 5e4

  # Number of simulations for multiple testing correction.
  # hide
  cgwasenv$.SIMUL_N <- cgwasenv$.SIMUL_DEP * cgwasenv$.IND_SNP_N / cgwasenv$.SIMUL_SNP_N

  # Study-wide significant threshold
  # >0 & <=5e-6, default 0.05/indSNPN
  # expose
  if (pStudy <= 0 | pStudy > 5e-6) {
    return(message("ERROR: parament pStudy invalid!"))
  } else {
    cgwasenv$.P_THRD_STUDY <- pStudy
  }

  # Suggestive significant threshold
  # >0 & <=1e-4, default 1/indSNPN
  # expose
  if (pSugst <= 0 | pSugst > 1e-4) {
    return(message("ERROR: parament pSugst invalid!"))
  } else {
    cgwasenv$.P_THRD_SUGST <- pSugst
  }

  # main effect threshold
  # >=3/1e4, default 3/indSNPN
  # expose
  if (pMain < 3/1e4) {
    return(message("ERROR: parament pMain invalid!"))
  } else {
    cgwasenv$.P_MAIN_EFFECT <- pMain
  }

  # min EbICoW power increase ratio
  # >0, default 1
  # expose
  if (ebicoPwrInc <= 0) {
    return(message("ERROR: parament ebicoPwrInc invalid!"))
  } else {
    cgwasenv$.MIN_EbICo_POWER_INC <- ebicoPwrInc
  }

  # min correlation difference
  # >0, default 0.05
  # expose
  if (minCorDiff <= 0) {
    return(message("ERROR: parament ebicoPwrInc invalid!"))
  } else {
    cgwasenv$.MIN_CORR_DIFF <- minCorDiff
  }

  # high correlation restriction (squared)
  # >0 & <1, default 0.5
  # expose
  if (hiCorRestc <= 0 | hiCorRestc >= 1) {
    return(message("ERROR: parament hiCorRestc invalid!"))
  } else {
    cgwasenv$.HIGH_CORR2_RES <- hiCorRestc
  }

  # high correlation restriction
  # hide
  cgwasenv$.HIGH_CORR_RES <- sqrt(cgwasenv$.HIGH_CORR2_RES)

  # Equivalent sample size increase ratio
  # >0 & <1, default 0.5
  # expose
  if (sampleNInc <= 0 | sampleNInc >= 1) {
    return(message("ERROR: parament sampleNInc invalid!"))
  } else {
    cgwasenv$.SAMPLE_SIZE_INC <- sampleNInc
  }

  # TWT stratification cutoff
  # >0 & <1, default 10^(seq(0,ceiling(log10(indSNPN))-1,-1/3))[-1]
  # expose
  if (min(twtStCo) <= 0 | max(twtStCo) >= 1) {
    return(message("ERROR: parament twtStCo invalid!"))
  } else {
    cgwasenv$.TWT_STRAT_CUT <- twtStCo
  }

  # LOESS interval p
  # decreasing order, >=10/indsnpn & <=0.1, default c(0.05,0.001)
  # expose
  if (min(LoInrP) < 10/cgwasenv$.IND_SNP_N | max(LoInrP) > 0.1) {
    return(message("ERROR: parament LoInrP invalid!"))
  } else {
    cgwasenv$.LOESS_INTER_P <- sort(LoInrP, decreasing = TRUE)
  }

  # LOESS span vector
  # length=length(interp)+1, >0 & <1, default c(0.03,0.1,0.75)
  # expose
  if (length(LoSpanV) != length(LoInrP)+1 | min(LoSpanV) <= 0 | max(LoSpanV) >= 1) {
    return(message("ERROR: parament LoSpanV invalid!"))
  } else {
    cgwasenv$.LOESS_SPAN_V <- LoSpanV
  }

  # Loci interval
  # integer, >=1e5 & <=1e6, default 2.5e5
  # expose
  if (lociInr < 1e5 | lociInr > 1e6) {
    return(message("ERROR: parament lociInr invalid!"))
  } else {
    cgwasenv$.LOCI_INTER <- lociInr
  }

  # create directory of intermediate result
  cgwasenv$.CGWAS_DETAIL_PATH <- file.path(cgwasenv$.CGWAS_DIR, "Details")
  cgwasenv$.CGWAS_iEbICoW_PATH <- file.path(cgwasenv$.CGWAS_DIR, "Details", "i-EbICoW")
  cgwasenv$.CGWAS_RESULT_PATH <- file.path(cgwasenv$.CGWAS_DIR, "Result")
  dir.create(cgwasenv$.CGWAS_TEMPDATA_PATH, recursive = T, showWarnings = F)
  dir.create(cgwasenv$.CGWAS_COLDATA_PATH, recursive = T, showWarnings = F)
  dir.create(cgwasenv$.CGWAS_RESULT_PATH, recursive = T, showWarnings = F)

  # Send R Output to a File
  sink(file.path(cgwasenv$.CGWAS_RESULT_PATH, 'cgwas.log'),
       append=FALSE, split=TRUE)
  # output basic paramenters
  paramOutput(cgwasenv)

  # Step1 Data formating
  ts <- Sys.time()
  step1(cgwasenv)
  print(paste0('Step 1 takes ',
               as.numeric(difftime(Sys.time(), ts, units = 'secs')),
               'Secs.'))

  # Step2 GetI
  ts <- Sys.time()
  step2(cgwasenv)
  print(paste0('Step 2 takes ',
               as.numeric(difftime(Sys.time(), ts, units = 'secs')),
               'Secs.'))

  # Step3 GetPsi
  ts <- Sys.time()
  step3(cgwasenv)
  print(paste0('Step 3 takes ',
               as.numeric(difftime(Sys.time(), ts, units = 'secs')),
               'Secs.'))

  # Step4 i-EbICoW
  ts <- Sys.time()
  step4(cgwasenv)
  print(paste0('Step 4 takes ',
               as.numeric(difftime(Sys.time(), ts, units = 'secs')),
               'Secs.'))

  # Step5 Truncated Wald Test
  ts <- Sys.time()
  step5(cgwasenv)
  print(paste0('Step 5 takes ',
               as.numeric(difftime(Sys.time(), ts, units = 'secs')),
               'Secs.'))

  # Step6 Summary of C-GWAS
  ts <- Sys.time()
  step6(cgwasenv)
  print(paste0('Step 6 takes ',
               as.numeric(difftime(Sys.time(), ts, units = 'secs')),
               'Secs.'))

  # R Output restoring default values
  sink()
}
