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
#' @param cgwas_dir blablabla
#' @param cgwas_data_dir blablabla
#' @param trait_data_name blablabla
#' @param trait_name blablabla
#' @param trait_num blablabla
#' @param snp_column_index blablabla
#' @param assoc_column_index blablabla
#' @param paral_num blablabla
#' @param heteroprop blablabla
#' @param Ccp blablabla
#' @param Manco blablabla
#' @param Locisep blablabla
#'
#' @examples
#'
#' cgwas_dir <- getwd()
#' cgwas_data_dir <- system.file("extdata", package = "CGWAS")
#' trait_data_name <- c("phe1.assoc", "phe2.assoc", "phe3.assoc", "phe4.assoc", "phe5.assoc")
#' trait_name <- c("phe1", "phe2", "phe3", "phe4", "phe5")
#' trait_num <-c(10015, 10015, 10015, 10015, 10015)
#' snp_column_index <- c(1,2,3)
#' assoc_column_index <- c(4,5)
#' paral_num <- 1
#' cgwas(cgwas_dir, cgwas_data_dir, trait_data_name, trait_name, trait_num,
#'       snp_column_index, assoc_column_index, paral_num)
#'
#' @export
cgwas <- function(cgwas_dir, cgwas_data_dir,
                  trait_data_name, trait_name, trait_num, snp_column_index,
                  assoc_column_index, paral_num,
                  heteroprop=1, Ccp=0.1, Manco=0.01, Locisep=2.5e5) {
  print("CGWAS version 2.3.0")

  cgwasenv <- list()

  cgwasenv$.CGWAS_DIR <- cgwas_dir
  cgwasenv$.CGWAS_DATA_DIR <- cgwas_data_dir
  cgwasenv$.TRAIT_DATA_NAME <- trait_data_name
  cgwasenv$.TRAIT_NAME <- trait_name
  cgwasenv$.TRAIT_EFFECT_SIZE <- trait_num
  cgwasenv$.SNP_COLUMN_INDEX <- snp_column_index
  cgwasenv$.ASSOC_COLUMN_INDEX <- assoc_column_index
  cgwasenv$.PARAL_NUM <- min(parallel::detectCores() - 1, paral_num)
  cgwasenv$.TRAIT_NUM <- length(cgwasenv$.TRAIT_NAME)
  cgwasenv$.HETEROPROP <- heteroprop
  cgwasenv$.CCP <- Ccp
  cgwasenv$.MANCO <- Manco
  cgwasenv$.LOCISEP <- Locisep

  cgwasenv$.CGWAS_TEMPDATA_PATH <- file.path(cgwasenv$.CGWAS_DIR, "Tempdata")
  cgwasenv$.CGWAS_COLDATA_PATH <- file.path(cgwasenv$.CGWAS_DIR, "Tempdata", "Coldata")
  cgwasenv$.CGWAS_INFCOR_PATH <- file.path(cgwasenv$.CGWAS_DIR, "Tempdata", "Infcor")
  cgwasenv$.CGWAS_GMA_PATH <- file.path(cgwasenv$.CGWAS_DIR, "Tempdata", "GMA")
  cgwasenv$.CGWAS_RESULT_PATH <- file.path(cgwasenv$.CGWAS_DIR, "Result")
  cgwasenv$.CGWAS_BARPLOTS_PATH <- file.path(cgwasenv$.CGWAS_DIR, "Result", "Barplots")

  cgwasenv

  dir.create(cgwasenv$.CGWAS_TEMPDATA_PATH, recursive = T, showWarnings = F)
  dir.create(cgwasenv$.CGWAS_COLDATA_PATH, recursive = T, showWarnings = F)
  dir.create(cgwasenv$.CGWAS_INFCOR_PATH, recursive = T, showWarnings = F)
  dir.create(cgwasenv$.CGWAS_GMA_PATH, recursive = T, showWarnings = F)
  dir.create(cgwasenv$.CGWAS_RESULT_PATH, recursive = T, showWarnings = F)
  dir.create(cgwasenv$.CGWAS_BARPLOTS_PATH, recursive = T, showWarnings = F)

  # Step1 StatExtraction
  ts <- Sys.time()
  step1(cgwasenv)
  print("Step1 finished")
  print(difftime(Sys.time(),ts))

  # Step2 InflationCorr
  ts <- Sys.time()
  step2(cgwasenv)
  print("Step2 finished")
  print(difftime(Sys.time(),ts))

  # Step3 CorEstimation
  ts <- Sys.time()
  step3(cgwasenv)
  print("Step3 finished")
  print(difftime(Sys.time(),ts))

  # Step4 preGMA
  ts <- Sys.time()
  step4(cgwasenv)
  print("Step4 finished")
  print(difftime(Sys.time(),ts))

  # Step5 GMA
  ts <- Sys.time()
  step5(cgwasenv)
  print("Step5 finished")
  print(difftime(Sys.time(),ts))

  # Step6 GMA+QFSW-2
  ts <- Sys.time()
  step6(cgwasenv)
  print("Step6 finished")
  print(difftime(Sys.time(),ts))

  # Step7 QFSW
  ts <- Sys.time()
  step7(cgwasenv)
  print("Step7 finished")
  print(difftime(Sys.time(),ts))

  # Step8 Simulation
  ts <- Sys.time()
  step8(cgwasenv)
  print("Step8 finished")
  print(difftime(Sys.time(),ts))

  # Step9 Visualization
  ts <- Sys.time()
  step9(cgwasenv)
  print("Step9 finished")
  print(difftime(Sys.time(),ts))

  # Step10 Summary
  ts <- Sys.time()
  step10(cgwasenv)
  print("Step10 finished")
  print(difftime(Sys.time(),ts))

  # Step11 Sumbarplots
  ts <- Sys.time()
  step11(cgwasenv)
  print("Step11 finished")
  print(difftime(Sys.time(),ts))
}
