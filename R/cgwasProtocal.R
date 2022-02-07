 #' A Reference Class to save global variables
 #' @noRd
CGWAS_ENV <- setRefClass("CGWAS_ENV",
                         fields = c(".GWAS_FILE_PATH",".SNP_FILE_PATH",
                                    ".TRAIT_NAME",".TRAIT_NUM",".SNP_N",
                                    ".EXCLUDE_NA",".MAF_FILE_PATH",
                                    ".MAF_FILE_EXIST",".CGWAS_DIR",
                                    ".KEEP_EbICoW",".PARAL_NUM",".IND_SNP_N",
                                    ".SIMUL_DEP",".SIMUL_SNP_N",".SIMUL_N",
                                    ".P_THRD_STUDY",".P_THRD_SUGST",
                                    ".P_MAIN_EFFECT",".MIN_EbICo_POWER_INC",
                                    ".MIN_CORR_DIFF",".HIGH_CORR2_RES",
                                    ".HIGH_CORR_RES",".SAMPLE_SIZE_INC",
                                    ".TWT_STRAT_CUT",".LOESS_INTER_P",
                                    ".LOESS_SPAN_V",".LOCI_INTER",
                                    ".CGWAS_DETAIL_PATH",".CGWAS_iEbICoW_PATH",
                                    ".CGWAS_RESULT_PATH"),
                         inheritPackage = TRUE)

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
#' The summary and immediate results of C-GWAS analysis are saved in the output
#' directory \code{outputPath} defined by user. There are two directorys in this
#' directory: \code{Result/} and \code{Details/}. The former contains the
#' summary results, including 1) Manhattan plots \code{CGWAS-GWAS.jpg} and 2) Q-Q
#' plots \code{CGWASminpQQ.jpg} of C-GWAS and MinGWAS, 3) suggestively significant
#' SNPs table \code{SummarySugSigSNP.csv}, and 4) p-values of all SNPs \code{C-GWAS.p}.
#' The latter contains the immediate result, including 1) result of getI
#' function applied to all GWASs \code{SummaryGetI.txt}; 2) result of getPsi function
#' applied to all GWAS pairs \code{SummaryGetPsi.txt}, 3) summary of result of
#' i-EbICoW \code{EbICoW.txt}, 4) result of getPsi function applied to all EbICoW
#' GWAS pairs \code{SummaryEbICoWGetPsi.txt}, 5) records of all i-EbICoW iterations
#' \code{Summaryi-EbICoW.txt}, 6) C-GWAS LOESS model samples simulated in getCoef
#' function \code{NullCGWAScorrection.txt}, 7) performance of getCoef based
#' correction in C-GWAS simulation \code{NullCGWASdistribution.jpg}, 8) MinGWAS LOESS
#' model samples simulated in getCoef function \code{NullMinpcorrection.txt}, 9)
#' performance of getCoef based correction in MinGWAS simulation
#' \code{NullMinpdistribution.jpg}. If user choose to keep i-EbICoW output (Kio), all
#' summary statistics (effect sizes, *.beta and test statistics, *.stat) of
#' EbICoW GWASs would be kept in directory \code{Details/i-EbICoW}.
#'
#' @param gwasFilePath a string list containg the paths to GWASs summary files.
#' @param snpFilePath a string indicating path to SNP information file. The SNP
#' information file has three columns, which represent CHR, BP and SNP
#' respectively.
#' @param outputPath a string indicating path to result output directory.
#' @param traitName a string list of trait names.
#' @param exNa logical. If \code{TRUE}, SNPs with NA in at least one GWAS will
#' be removed; if If \code{FALSE}, the NA is replaced with BETA=0 and P=1.
#' @param mafFilePath a string indicating path to MAF file.
#' @param keepIEb logical. If \code{TRUE},, BETA and STAT of all i-EbICoW
#' combinations will be saved in \code{Details/i-EbICoW}.
#' @param threadN number of threads to be used for parallel computing. The
#' default value is half of the maximum number of cores
#' on the hardware.
#' @param indSNPN independent SNPs number, which should not be smaller than
#' 10,000 and a multiple of 500. The default value is 1e6.
#' @param simulDep an integer specifying quantile simulation depth, which should
#' not be smaller than 100 and a multiple of 100. The default value is 100.
#' @param pStudy study-wide significant threshold, which should
#' take values in \code{(0,5e-6]}. The default value is 0.05/indSNPN
#' @param pSugst suggestive significant threshold, which should
#' take values in \code{(0,1e-4]}. The default value is 1/indSNPN
#' @param pMain main effect threshold, which should take values in \code{(0,3/1e4]}.
#' The default value is 3/indSNPN.
#' @param ebicoPwrInc minimum EbICoW power increase ratio, which should greater
#' than 0. The default value is 1.
#' @param minCorDiff minimum correlation difference, which should greater than
#' 0. The default value is 0.05.
#' @param hiCorRestc high correlation restriction, which should take values in
#' \code{(0,1)}. The default value is 0.05.
#' @param sampleNInc equivalent sample size in EbICo increase ratio, which
#' should take values in \code{(0,1)}. The default value is 0.05.
#' @param twtStCo a vector of TWT stratification cutoffs, which should take
#' values in \code{(0,1)}. The default value is
#' \code{10^(seq(0,1-ceiling(log10(indSNPN)),-1/3))[-1]}.
#' @param LoInrP a vector of LOESS interval p, which should be in decreased
#' order and take values in \code{[10/indSNPN,0.1]}. The default value is
#' \code{c(0.05,0.001)}.
#' @param LoSpanV a vector of LOESS span, which should take values in \code{(0,1)}.
#' The default value is \code{c(0.03,0.1,0.75)}. The length of \code{LoSpanV}
#' should equal to \code{length(LoInrP)+1}.
#' @param lociInr an integer indicating loci interval, which should take values
#' in \code{[1e5,1e6]}. The default value is 2.5e5.
#'
#' @return Several files and figures will be generated in \code{outputPath}.
#'
#' \code{CGWAS-GWAS.jpg}: Combined Manhattan plots displays adjusted C-GWAS
#' p-value (upper parts) and adjusted MinGWAS p-value (lower parts) of all SNPs
#' in –log10 scale. The study-wide significance is indicated using dashed lines
#' and the suggestive significance is indicated using solid lines. Loci with the
#' lead SNP passing the study-wide significant line on both sides are
#' highlighted using different colors and shapes. One side’s unique loci are
#' indicated in orange diamond. Two sides’ overlapped loci are indicated in
#' green cross.
#'
#' \code{CGWASminpQQ.jpg}: Adjusted p-values in –log10 scale (y-axis) of MinGWAS
#' (blue) and C-GWAS (orange) for all SNPs are plotted against their expected
#' uniformly distributed quantile in –log10 scale (x-axis) under the null.
#'
#' \code{SummarySugSigSNP.csv}: This table contains the SNP passing the
#' suggestive significant cut-off of either C-GWAS or MinGWAS. For each SNP, its
#' position (CHR and BP), rsID (SNP), adjusted C-GWAS p-value (C-GWAS-P) and
#' MinGWAS p-value (MinGWAS-P) is always provided. Minor allele frequency (MAF)
#' will be provided if user define available frequency data path \code{mafFilePath}.
#' All SNPs were clustered to different loci base on distance parameter
#' \code{lociInr}. Loci are distinguished using the descending order of significance
#' of lead SNP in each loci from C-GWAS (C-GWAS-Loci) and MinGWAS (MinGWAS-Loci)
#' respectively. NA indicates the SNP is not suggestive significant in C-GWAS or
#' MinGWAS.
#'
#' \code{C-GWAS.p}: This dataset contains position (CHR and BP), rsID (SNP),
#' adjusted C-GWAS p-value (C-GWASadjustedP), unadjusted C-GWAS p-value
#' (C-GWASrawP), adjusted MinGWAS p-value (MinGWASadjustedP), unadjusted MinGWAS
#' p-value (MinGWASrawP) of all SNPs. Minor allele frequency (MAF) will be also
#' included if user define available frequency data path \code{mafFilePath}.
#'
#' \code{SummaryGetI.txt}: This table contains mean chi-square (square of test
#' statistics) of all SNPs before adjustment (RawMeanX2), Genomic control value
#' (GClambda), getI result (EstInf) and mean chi-square of all SNPs after getI
#' based adjustment (AdjMeanX2) for each GWAS (traitName).
#'
#' \code{SummaryGetPsi.txt}: This table contains Pearson correlation coefficient
#' of test statistics from all SNPs (StatCor), getPsi result (Psi) and getPi
#' result (considering all SNPs, allPi; considering only significant SNPs,
#' sigPi) for all GWAS pairs (GWAS1 and GWAS2).
#'
#' \code{EbICoW.txt}: This table contains effect of EbICoW GWAS (define as mean
#' chi-square of all SNPs minus one, AllEff), Equivalent sample size (Ess),
#' EbICoW GWAS combination members (GWAS*), corresponding weight of effect size
#' in combination (bw*), corresponding weight of test statistics in combination
#' (tw*) and rounds of all iterations involved in combine (r*) for each EbICoW
#' GWAS (traitName).
#'
#' \code{SummaryEbICoWGetPsi.txt}: This table contains Pearson correlation
#' coefficient of test statistics from all SNPs (StatCor), getPsi result (Psi)
#' and getPi result (considering all SNPs, allPi; considering only significant
#' SNPs, sigPi) for all EbICoW GWAS pairs (GWAS1 and GWAS2).
#'
#' \code{Summaryi-EbICoW.txt}:  This table contains the details of each i-EbICoW
#' iteration round.
#' 1) Iteration round count (Round).
#' 2) Two GWASs used in the iteration (GWAS1 and GWAS2).
#' 3) Name of combined GWAS (NewGWAS).
#' 4) Type of Pi and h chosen by optimize function (Piall and hall, AEC; Pisig
#' and hsig, SEC; Pistb and hstb, STB).
#' 5) Whether pass evaluation of evaluate function (Evaluate).
#' 6) Effect of GWAS1, GWAS2 and NewGWAS (define as mean chi-square of all SNPs
#' minus one, Eff1X2, Eff2X2	and NewEffX2).
#' 7) getPsi and getPi result (Psi; considering all SNPs, allPi; considering
#' only significant SNPs, sigPi).
#' 8) Significant SNP number using parameter threshold \code{pMain}, including
#' in GWAS1 (sigSNP1N), GWAS2 (sigSNP2N), MinGWAS of GWAS1 and GWAS2
#' (sigSNPminN), combined GWAS using AEC (sigSNPAECN), combined GWAS using SEC
#' (sigSNPSECN), combined GWAS using STB (sigSNPSTBN), combined GWAS using Wald
#' test (sigSNPWDN), intersection of AEC combined GWAS and GWAS1-GWAS2 union
#' (sigSNPAECuniN), intersection of AEC combined GWAS and GWAS1-GWAS2 union
#' (sigSNPSECuniN), intersection of AEC combined GWAS and GWAS1-GWAS2 union
#' (sigSNPSTBuniN), intersection of AEC combined GWAS and GWAS1-GWAS2
#' intersection (sigSNPAECinterN), intersection of AEC combined GWAS and
#' GWAS1-GWAS2 intersection (sigSNPSECinterN), intersection of AEC combined GWAS
#' and GWAS1-GWAS2 intersection (sigSNPSTBinterN). 9) Equivalent sample size of
#' GWAS1, GWAS2 and NewGWAS (Ess1, Ess2 and EssNew), 10) Weight of effect size
#' used in combining GWAS1 and GWAS2 chosen by optimize function (bw1 and bw2),
#' 11) Weight of test statistics used in combining GWAS1 and GWAS2 chosen by
#' optimize function (tw1 and tw2).
#'
#' \code{NullCGWAScorrection.txt}: This dataset contains LOESS model samples
#' (ExpectedQuantile and ObservedP) simulated in getCoef function using in
#' adjusting C-GWAS result.
#'
#' \code{NullCGWASdistribution.jpg}: Simulated null distributions of the
#' unadjusted p-values (gray) and p-values corrected using the getCoef function
#' (orange) are compared with the uniform distribution (black) for C-GWAS
#' simulation. 95% confidence interval is indicated using dashed lines.
#'
#' \code{NullMinpcorrection.txt}: This dataset contains LOESS model samples
#' (ExpectedQuantile and ObservedP) simulated in getCoef function using in
#' adjusting MinGWAS result.
#'
#' \code{NullMinpdistribution.jpg}: Simulated null distributions of the
#' unadjusted p-values (gray) and p-values corrected using the getCoef function
#' (orange) are compared with the uniform distribution (black) for MinGWAS
#' simulation. 95% confidence interval is indicated using dashed.
#'
#' @examples
#'
#' ExDataDir <- system.file("extdata", package = "CGWAS")
#' gwasFileName <- c("AlL-ChL.assoc", "AlL-Sn.assoc", "EnL-AlL.assoc",
#'                   "EnL-ChL.assoc", "Ls-ChL.assoc", "N-EnL.assoc",
#'                   "N-Prn.assoc", "Prn-AlL.assoc", "Prn-EnL.assoc",
#'                   "Prn-Sn.assoc", "Sn-ChL.assoc", "Sn-Ls.assoc")
#' gwasFilePath <- file.path(ExDataDir, gwasFileName)
#' snpFilePath <- file.path(ExDataDir, 'SnpInfo')
#' outputPath <- getwd()
#' traitName <- c("AlL-ChL", "AlL-Sn", "EnL-AlL",
#'                "EnL-ChL", "Ls-ChL", "N-EnL",
#'                "N-Prn", "Prn-AlL", "Prn-EnL",
#'                "Prn-Sn", "Sn-ChL", "Sn-Ls")
#' mafFilePath <- file.path(ExDataDir, 'MAF')
#' indSNPN = 1e5
#'
#' cgwas(gwasFilePath, snpFilePath, outputPath,
#'       traitName = traitName, mafFilePath = mafFilePath, indSNPN = indSNPN)
#'
#' @export
cgwas <- function(gwasFilePath,
                  snpFilePath,
                  outputPath,
                  traitName = basename(gwasFilePath[1]),
                  exNa = TRUE,
                  mafFilePath = NULL,
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
\\____/       \\____/   |__/|__//_/  |_|/____/      version 1.9.1\n\n")

  cgwasenv <- CGWAS_ENV$new()

  # A string list of path to GWASs summary data files.
  # expose
  cgwasenv$.GWAS_FILE_PATH <- gwasFilePath

  # A string list of path to GWASs summary data files.
  # expose
  cgwasenv$.SNP_FILE_PATH <- snpFilePath
  # cgwasenv$.ASSOC_COLUMN_INDEX <- assocColInx

  # A list of trait names.
  # expose
  cgwasenv$.TRAIT_NAME <- traitName

  # Trait numbers.
  # hide
  cgwasenv$.TRAIT_NUM <- length(cgwasenv$.TRAIT_NAME)

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
  if (pMain <= 0 | pMain > 3/1e4) {
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

  # create directory of output
  cgwasenv$.CGWAS_DETAIL_PATH <- file.path(cgwasenv$.CGWAS_DIR, "Details")
  cgwasenv$.CGWAS_iEbICoW_PATH <- file.path(cgwasenv$.CGWAS_DIR, "Details", "i-EbICoW")
  cgwasenv$.CGWAS_RESULT_PATH <- file.path(cgwasenv$.CGWAS_DIR, "Result")
  dir.create(cgwasenv$.CGWAS_DETAIL_PATH, recursive = T, showWarnings = F)
  dir.create(cgwasenv$.CGWAS_iEbICoW_PATH, recursive = T, showWarnings = F)
  dir.create(cgwasenv$.CGWAS_RESULT_PATH, recursive = T, showWarnings = F)
  if (file.exists(file.path(cgwasenv$.CGWAS_RESULT_PATH, 'LogFile'))) {
    file.remove(file.path(cgwasenv$.CGWAS_RESULT_PATH, 'LogFile'))
  }

  # output basic paramenters
  paramOutput(cgwasenv)

  # Step1 Data formating
  ts <- Sys.time()
  step1(cgwasenv)
  logOutput('Step 1 takes ',
            round(as.numeric(difftime(Sys.time(), ts, units = 'secs')), 1),
            ' Secs.\n\n\n', cgwasenv = cgwasenv)

    # Step2 GetI
  ts <- Sys.time()
  step2(cgwasenv)
  logOutput('Step 2 takes ',
            round(as.numeric(difftime(Sys.time(), ts, units = 'secs')), 1),
            ' Secs.\n\n\n', cgwasenv = cgwasenv)

  # Step3 GetPsi
  ts <- Sys.time()
  step3(cgwasenv)
  logOutput('Step 3 takes ',
            round(as.numeric(difftime(Sys.time(), ts, units = 'secs')), 1),
            ' Secs.\n\n\n', cgwasenv = cgwasenv)

  # Step4 i-EbICoW
  ts <- Sys.time()
  step4(cgwasenv)
  logOutput('Step 4 takes ',
            round(as.numeric(difftime(Sys.time(), ts, units = 'secs')), 1),
            ' Secs.\n\n\n', cgwasenv = cgwasenv)

  # Step5 Truncated Wald Test
  ts <- Sys.time()
  step5(cgwasenv)
  logOutput('Step 5 takes ',
            round(as.numeric(difftime(Sys.time(), ts, units = 'secs')), 1),
            ' Secs.\n\n\n', cgwasenv = cgwasenv)

  # Step6 Summary of C-GWAS
  ts <- Sys.time()
  step6(cgwasenv)
  logOutput('Step 6 takes ',
            round(as.numeric(difftime(Sys.time(), ts, units = 'secs')), 1),
            ' Secs.\n\n\n', cgwasenv = cgwasenv)
}
