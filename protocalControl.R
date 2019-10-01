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
#' Each data indicated by \code{ST_Filename} should be a white-space (space or tab)
#' delimited file with header line. The first five columns are mandatory:
#' Chromesome, Base-pair position, SNP id, BETA, P-vale.
#'
#' If \code{traitName} is NA, it will be filled by \code{ST_Filename}.
#' If \code{Gt} is NA, genome-wide threshold will be set as max(5e-8, 0.05/Rsnpnum).
#' If \code{St} is NA, suggestive threshold will be set as max(1.6e-6, 1/Rsnpnum).
#' \code{Rsnpnum} is set as the total number of SNPs.
#'
#' \code{paral} set the function running in parallel mode and limit the maximum
#' number of unsable CPUs if it is set as an integer greater than 1. Otherwise,
#' the function will run in serial mode.
#'
#' @param ST_Datapath The path of single trait GWAS results.
#' @param ST_Filename A character vector of single trait GWAS result file name.
#' @param CGPath The work path of CGWAS.
#' @param traitName A character vector of trait names.
#' @param simulNum Simulation times. Not less than 1000.
#' @param Gt Genome-wide threshold.
#' @param St Suggestive threshold.
#' @param Ht Display threshold.
#' @param Clv Inflation Factors.
#' @param Bdd BETA decimal digits.
#' @param SimuReg Simulation estimation proportion.
#' @param paral Number of parallel.
#'
#' @examples
#'
#' ST_Datapath <- dirname(system.file("extdata", "result_P_1_2.assoc.linear.arranged.part", package = "CGWAS"))
#' ST_Filename <- c("result_P_1_2.assoc.linear.arranged.part", "result_P_1_3.assoc.linear.arranged.part", "result_P_1_4.assoc.linear.arranged.part", "result_P_1_5.assoc.linear.arranged.part", "result_P_1_6.assoc.linear.arranged.part")
#' CGPath <- getwd()
#' CGWAS(ST_Datapath, ST_Filename, CGPath, paral = 1)
#'

#' @export
CGWAS <- function(ST_Datapath, ST_Filename, CGPath, traitName = NA, simulNum = 1000, Gt = NA, St = NA, Ht = 1e-2, Clv = NA, Bdd = 4, SimuReg = 0.8, paral = 1){
  #set global variants
  cgwasenv <<- new.env()

  cgwasenv$.cgwas_dir <- CGPath
  cgwasenv$.traitName <- traitName
  cgwasenv$.simulNum <- simulNum
  cgwasenv$.Gt <- Gt
  cgwasenv$.St <- St
  cgwasenv$.paral <- paral
  cgwasenv$.cgwas_dir <- paste0(cgwasenv$.cgwas_dir, '/')
  ST_Datapath <- paste0(ST_Datapath, '/')
  if (sum(is.na(cgwasenv$.traitName)))
    cgwasenv$.traitName <- ST_Filename
  cgwasenv$.pheNum <- length(cgwasenv$.traitName)
  if (!is.na(Bdd)) Fdd <- T

  cgwasenv$.cgwas_tempdata = file.path(cgwasenv$.cgwas_dir, "Tempdata")
  cgwasenv$.cgwas_rowdata = file.path(cgwasenv$.cgwas_dir, "Tempdata", "Rowdata")
  cgwasenv$.cgwas_simuldata = file.path(cgwasenv$.cgwas_dir, "Tempdata", "Simudata")
  cgwasenv$.cgwas_allp = file.path(cgwasenv$.cgwas_dir, "Tempdata", "Simudata", "Allp")
  cgwasenv$.cgwas_minp = file.path(cgwasenv$.cgwas_dir, "Tempdata", "Simudata", "Minp")
  cgwasenv$.cgwas_result = file.path(cgwasenv$.cgwas_dir, "Result")
  cgwasenv$.cgwas_correction = file.path(cgwasenv$.cgwas_dir, "Result", "Correction")

  dir.create(cgwasenv$.cgwas_dir, recursive = T, showWarnings = F)
  dir.create(cgwasenv$.cgwas_tempdata, recursive = T, showWarnings = F)
  dir.create(cgwasenv$.cgwas_rowdata, recursive = T, showWarnings = F)
  dir.create(cgwasenv$.cgwas_simuldata, recursive = T, showWarnings = F)
  dir.create(cgwasenv$.cgwas_allp, recursive = T, showWarnings = F)
  dir.create(cgwasenv$.cgwas_minp, recursive = T, showWarnings = F)
  dir.create(cgwasenv$.cgwas_result, recursive = T, showWarnings = F)
  dir.create(cgwasenv$.cgwas_correction, recursive = T, showWarnings = F)
  print('subfolder build')

  timeProtocalStart <- Sys.time()

  # Data preparation
  step1(ST_Datapath, ST_Filename, Clv, Bdd, Fdd)

  # simulation
  step2_parallel()

  # Simulation FDR (note: TaskAmount>TraitAmount)
  step3_parallel(SimuReg)

  # Actual data CGWAS
  step4_parallel()

  # Optimization inflation correction & Raw hits
  step5(Ht)

  # Merged Manhattan & QQ plots
  step6(m_wid=3500,m_hei=3500,sing_col=c('chartreuse2','chartreuse4'),comb_col=c('darkgoldenrod2','darkgoldenrod4'),speeduppoints=T,q_wid=3500,q_hei=2400)

  # Suggestive & Genomewide & Multiple test correction & TopSNPs
  step7()
  timeProtocalEnd <- Sys.time()
  print(paste0("whole protocal use time:"))
  print(timeProtocalEnd - timeProtocalStart)
}
