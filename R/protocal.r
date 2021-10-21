#' packagename: combine GWAS
#'
#' This package provides a toolkit for combine GWAS
#'
#' Integrating multiple single trait GWAS results
#'
#' \code{CGWAS} combine multiple single trait GWAS results which indicated by its
#' arguments to calculate integrate association result.
#'
#' \code{CGWAS} combine multiple single trait GWAS results to calculate integrate
#' association result from a set of trait combinations for each SNP specifically.
#' The method increase power detecting SNP having effect (may be weak in GWAS
#' for one trait) on multiple traits and type I error was controlled through
#' a large scale simulation analyses. Significant CGWAS result contains both
#' single trait GWAS hits and multiple traits combined GWAS hits after multiple test.
#'
#' \code{snpFilename} shoule be a white-space (space or tab) delimited file
#' with header line. This file desctribes the SNP information and the first three
#' columun must be Chromesome, Base-pair position, SNP id.
#'
#' Each data indicated by \code{gwasFilename} should be a white-space (space or tab)
#' delimited file with header line. The first two columns are BETA, P-vale and
#' the order of each line should align with the order of the \code{snpFilename}
#' file SNPs.
#'
#' If \code{traitName} is NA, it will be filled by \code{gwasFilename}.
#'
#' If \code{Gt} is NA, genome-wide threshold will be set as max(5e-8, 0.05/Rsnpnum).
#'
#' If \code{St} is NA, suggestive threshold will be set as max(1.6e-6, 1/Rsnpnum).
#'
#' \code{Rsnpnum} is set as the total number of SNPs.
#'
#' \code{paral} set the function running in parallel mode and limit the maximum
#' number of unsable CPUs if it is set as an integer greater than 1. Otherwise,
#' the function will run in serial mode.
#'
#' @param gwasFileDir blablabla
#' @param gwasFileName blablabla
#' @param snpColInx blablabla
#' @param assocColInx blablabla
#' @param mafFilePath blablabla
#' @param threadN blablabla
#' @param outputPath blablabla
#' @param traitName blablabla
#' @param pt blablabla
#' @param simulN blablabla
#'
#' @examples
#'
#' gwasFileDir <- system.file("extdata", package = "CGWAS")
#' gwasFileName <- c("phe1.assoc", "phe2.assoc", "phe3.assoc", "phe4.assoc", "phe5.assoc")
#' snpColInx <- c(1, 2, 3)
#' assocColInx <- c(4, 5)
#' mafFilePath <- NULL
#' threadN <- 4
#' outputPath <- getwd()
#' traitName <- c("phe1", "phe2", "phe3", "phe4", "phe5")
#' pt <- NULL
#' simulN <- 200
#'
#' CGWAS(gwasFileDir, gwasFileName, snpColInx, assocColInx,
#'       mafFilePath, threadN, outputPath, traitName, pt, simulN)
#'
#' @export
CGWAS <- function(gwasFileDir,
                  gwasFileName,
                  snpColInx,
                  assocColInx,
                  mafFilePath = NULL,
                  threadN = ceiling(doparallel::detectCores() / 2),
                  outputPath = getwd(),
                  traitName = NULL,
                  pt = NULL,
                  simulN = 2e3) {

  print("---- C-GWAS version 1.8 start ----")

  cgwasenv <- list()

  # --------------------- Basic parameters --------------------- #

  # Dirctory containing single trait association file
  cgwasenv$.GWAS_FILE_DIR <- gwasFileDir

  # Single trait association file name
  cgwasenv$.GWAS_FILE_NAME <- gwasFileName

  # Column number of SNP index (CHR & BP & SNP)
  cgwasenv$.SNP_COLUMN_INDEX <- snpColInx

  # Column number of association data (BETA & P)
  cgwasenv$.ASSOC_COLUMN_INDEX <- assocColInx

  # MAF file path
  cgwasenv$.MAF_FILE_PATH <- mafFilePath

  # MAF file exist?
  cgwasenv$.MAF_FILE_EXIST <- !is.null(mafFilePath)

  # Parallel task number
  cgwasenv$.PARAL_NUM <- min(parallel::detectCores() - 1, threadN)

  # Work dirctory of C-GWAS
  cgwasenv$.CGWAS_DIR <- outputPath

  # Trait name
  if (is.null(traitName)) {
    cgwasenv$.TRAIT_NAME <- gwasFileName
  } else {
    cgwasenv$.TRAIT_NAME <- traitName
  }

  # Traits number
  cgwasenv$.TRAIT_NUM <- length(cgwasenv$.TRAIT_NAME)

  # main effect threshold
  # 建议 pt>1 / 独立SNP数量，否则出错概率增大
  # 建议 pt<10 / 独立SNP数量，否则结果可能变差
  # independent SNP size = 1e6 * (snp.N / 1e7)
  if (is.null(pt)) {
    snpn.gwas.tmp <- length(count.fields(file.path(cgwasenv$.GWAS_FILE_DIR, cgwasenv$.GWAS_FILE_NAME[1]),
                                         '\n')) - 1
    cgwasenv$.P_THRD <- 1 / (1e6 * snpn.gwas.tmp / 1e7)
  } else {
    cgwasenv$.P_THRD <- pt
  }

  # Simulation times
  cgwasenv$.SIMUL_N <- simulN

  # --------------------- Advanced parameters --------------------- #

  # min SNP number threshold
  cgwasenv$.SNP_N_THD <- 5e5

  # FDR test set
  cgwasenv$.FDR_SET <- c(seq(0.0001, 0.0009, 0.0001), seq(0.001, 0.05, 0.001))

  # minimum EbICo power increase ratio
  cgwasenv$.MIN_EbiCo_POWER <- 1.1

  # minimum correlation difference
  cgwasenv$.MIN_CORR_DIFF <- 0.05

  # high correlation restriction (squared)
  cgwasenv$.HIGH_CORR2_RES <- 0.7
  cgwasenv$.HIGH_CORR_RES <- sqrt(cgwasenv$.HIGH_CORR2_RES)

  # Equivalent sample size increase ratio
  cgwasenv$.SAMPLE_SIZE_INC <- 0.7

  # SWaT stratification cutoff
  cgwasenv$.SWaT_STRAT_CUT <- 10^(seq(0,-5,-1/3))[-1]

  # Simulated SNP number
  cgwasenv$.SIMUL_SNP_N <- 5e4

  # Independent SNP number
  cgwasenv$.IND_SNP_N <- 1e6

  # LOESS span vector
  cgwasenv$.LOESS_SPAN_V <- c(0.03, 0.1, 0.75)

  # LOESS interval p
  cgwasenv$.LOESS_INTER_P <- c(0.05, 0.001)

  # Loci interval
  cgwasenv$.LOCI_INTER <- 2.5e5

  cgwasenv$.CGWAS_TEMPDATA_PATH <- file.path(cgwasenv$.CGWAS_DIR, "Tempdata")
  cgwasenv$.CGWAS_COLDATA_PATH <- file.path(cgwasenv$.CGWAS_DIR, "Tempdata", "Coldata")
  cgwasenv$.CGWAS_RESULT_PATH <- file.path(cgwasenv$.CGWAS_DIR, "Result")

  # create directory of intermediate result
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
