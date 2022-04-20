 #' A Reference Class to save global variables
 #' @noRd
CGWAS_ENV <- setRefClass("CGWAS_ENV",
                         fields = c(".GWAS_FILE_PATH",".SNP_FILE_PATH",
                                    ".TRAIT_NAME",".TRAIT_NUM",".SNP_N",
                                    ".EXCLUDE_NA",".MRAF_FILE_PATH",
                                    ".MRAF_FILE_EXIST",".CGWAS_DIR",
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
#' Combining GWAS summary statistics of multiple potentially related traits.
#'
#' C-GWAS is a powerful method for combining GWAS summary statistics of multiple
#' potentially related traits and detect SNPs with multi-trait effects. C-GWAS
#' begins with GWASs summary as inputs and outputs a single vector of combined
#' p-values testing if the null is deviated. For each SNP, the null is the
#' absence of any effect on all traits, and the alternative is that its effect
#' deviates from 0 for at least one trait. C-GWAS integrates two different
#' statistical methods with complementary statistical features to ensure the
#' optimal power under various and complex scenarios while keeping a stable
#' study-wide type-I error rate. The first method is called iterative effect
#' based inverse covariance weighting (\code{i-EbICoW}) and the second method
#' is called truncated Wald test (\code{TWT}). C-GWAS controls the study-wide
#' type-I error rate in an empirical manner via simulations and adjust the
#' resultant p-values in such a way that they are directly comparable with those
#' from traditional GWAS of a single trait.
#'
#' The final and intermediate results of C-GWAS analysis are saved in the output
#' directory \code{outputPath}. Two new folders will be created in this
#' directory: \code{Results/} and \code{Details/}.
#'
#' @param gwasFilePath a vector of strings containing the paths to input GWASs
#' summary statistics files. These files should be in the space or tab delimited
#' format. Each file contains two and only two columns. The first column is the
#' regression betas of all SNPs and the second column is the P-values of all
#' SNPs. These files require a header line with two items: BETA and P. Note that
#' for each SNP, all betas must be based on the same reference allele.
#' @param snpFilePath string indicating the path to input SNP information file,
#' a space delimited file consisting of three columns representing chromosome,
#' base pair and SNP identifier. The file requires a header line with three
#' items CHR, BP and SNP in the specified order.
#' @param outputPath a string indicating the path to output directory. In this
#' directory, two new folders will be created, the \code{Results/} contains main
#' results and the \code{Details/} contains all intermediate results.
#' @param traitName a vector of strings of trait names. The order must be the
#' same as input GWASs. The basename of GWAS files will be set to trait names by
#' default.
#' @param exNa logical. If \code{TRUE}, SNPs with missing values in at least one
#' input GWAS will be removed and will not appear in the output files. All
#' removed SNPs will be listed in a file named \code{ExcludedSNP.txt}; if
#' \code{FALSE}, the NA is replaced with BETA=0 and P=1. Default is \code{TRUE}.
#' @param mrafFilePath a string indicating the path to the mean reference allele
#' frequency (MRAF) file, which requires a header MRAF, and contains one column
#' of mean frequency of the reference allele of each SNP. Note that the
#' reference allele must correspond to the beta in the input GWAS. Reference
#' allele frequencies are used to estimate the weights of input GWASs and all
#' intermediate GWASs. This MRAF file is recommended but not obligatory with
#' default \code{NULL}. If \code{NULL}, weights will still be estimated by
#' approximations.
#' @param keepIEb logical. If \code{TRUE}, all intermediate results from
#' \code{i-EbICoW} will be saved in \code{Details/i-EbICoW}. Default is
#' \code{FALSE}.
#' @param threadN number of threads/cores to be used for parallel computing. The
#' default is to automatically detect the number of CPU cores on the current
#' host and uses half of the available cores for parallel computing.
#' @param indSNPN the total number of unlinked SNPs of the input GWASs, the
#' minimal \code{indSNPN} is 1e4, e.g. for exon-wide association studies. The
#' default value is 1e6, which is fine for typical GWASs. \code{indSNPN} will be
#' automatically ceiling to nearest 500 in preprocess. This parameter is to
#' specify the number of SNPs for simulation analysis to adjust for extra
#' multiple testing induced by C-GWAS procedures.
#' @param simulDep the number of simulations, the minimal \code{simulDep} is
#' 50. \code{simulDep} will be automatically ceiling to nearest 50 in preprocess.
#' The default value is 100.
#' @param pStudy study-wide significant threshold after adjusting for extra
#' multiple testing due to C-GWAS procedures, which should take values in
#' \eqn{(0,5e-6]}. The default value is \code{0.05/indSNPN}. When the
#' \code{indSNPN} is default 1e6, the study-wide significance threshold of
#' C-GWAS equals to 5e-8 and corresponds to the genome-wide significance
#' threshold of traditional GWASs, meaning that C-GWAS results are directly
#' comparable with any standard GWAS.
#' @param pSugst study-wide suggestive significant threshold, which should take
#' values in \eqn{(0,1e-4]}. The default value is \code{1/indSNPN}. When
#' replication samples are available, this relaxed threshold is recommended for
#' the initial screening.
#' @param pMain an internal cut-off threshold of \code{i-EbICoW} to filter SNPs
#' with significant effect during iterations, range \eqn{(0,3/1e4]}. The default
#' value is \code{3/indSNPN}.
#' @param ebicoPwrInc an internal threshold of \code{i-EbICoW}, requesting that
#' the number of significant signals from EbICoW is at least \code{ebicoPwrInc}
#' times of that from the Wald test and/or that from the adjusted minimal P from
#' individual GWAS, range \eqn{(0, inf)}. This value is typically equal to or
#' slightly greater than 1. The default value is 1.1.
#' @param minCorDiff an internal threshold of \code{i-EbICoW}, specifying the
#' minimal difference between the background correlation (\code{Psi}) and the
#' effect correlation (\code{Pi}) for combining two GWASs in \code{i-EbICoW},
#' range \eqn{[0, 2)}. The default value is 0.05.
#' @param hiCorRestc an internal threshold of \code{i-EbICoW} to avoid abnormal
#' values when solving large \code{Psi} in the presence of strong collinearity,
#' range \eqn{(0,1)}. The default value is 0.5. GWASs with
#' \code{psi > hiCorRestc} in \code{i-EbICoW} iterations are forced combined.
#' @param sampleNInc an internal threshold of \code{i-EbICoW}, requesting that
#' the effective sample size (Ess) of EbICoW after combination should be larger
#' than the weighted sum of Ess of the two GWASs before combination. Typically,
#' \code{0 <= sampleNInc <= 1}, default is set to 0.5.
#' @param twtStCo a vector of internal thresholds of \code{TWT}, a vector of
#' \code{TWT} stratification cutoffs, range \eqn{(0,1)}. The default thresholds
#' are \code{10^(seq(0,1-ceiling(log10(indSNPN)),-1/3))[-1]}.
#' @param LoInrP a vector of internal interval cut-offs for fitting the LOESS
#' models, which should be in the descending order in the range
#' \code{[10/indSNPN,0.1]}. The default cut-offs are \code{c(0.05,0.001)}.
#' @param LoSpanV a vector of internal parameters for smoothing the LOESS curves
#' in the intervals specified by \code{LoInrP}, range \eqn{(0,1)}. The default
#' values are \code{c(0.03,0.1,0.75)}. Note that the length of
#' \code{LoSpanV = length(LoInrP)+1}.
#' @param lociInr an integer specifying the physical distance in base pair
#' between two SNPs that should be considered as in distinct loci, value range
#' \eqn{[1e5,1e6]}. The default value is \eqn{2.5e5}, i.e., \code{250 kbp}.
#'
#' @section Main Results:
#'
#' The \code{Results/} folder contains main results including:
#'
#' \describe{
#' \item{\code{CGWAS-GWAS.jpg}}{Combined Manhattan plots displays adjusted
#' C-GWAS p-value (upper parts) and adjusted MinGWAS p-value (lower parts) of
#' all SNPs at the \eqn{–log10} scale. The study-wide significance is indicated
#' using dashed lines and the suggestive significance is indicated using solid
#' lines. Loci with the lead SNP passing the study-wide significant threshold on
#' both sides are highlighted using different colors and shapes.}
#' \item{\code{CGWASminpQQ.jpg}}{Adjusted p-values at the \eqn{–log10} scale
#' (y-axis) of MinGWAS (blue) and C-GWAS (orange) for all SNPs are plotted
#' against their expected uniformly distributed quantiles at the \eqn{–log10}
#' scale (x-axis) under the null.}
#' \item{\code{SummarySugSigSNP.csv}}{This table contains all SNP passing the
#' suggestive significance threshold in either C-GWAS or MinGWAS. For each SNP,
#' its position (CHR and BP), rsID (SNP), adjusted C-GWAS p-value (C-GWAS-P) and
#' MinGWAS p-value (MinGWAS-P) are provided. Mean reference allele frequency
#' (MRAF) will be provided if user define available frequency data path
#' \code{mrafFilePath}. Distinct loci are distinguished by loci index. NA
#' indicates that the SNP did not pass the suggestive significance threshold.}
#' \item{\code{C-GWAS.p}}{This table contains position (CHR and BP), rsID (SNP),
#' adjusted C-GWAS p-value (C-GWASadjustedP), unadjusted C-GWAS p-value
#' (C-GWASrawP), adjusted MinGWAS p-value (MinGWASadjustedP), unadjusted MinGWAS
#' p-value (MinGWASrawP) of all SNPs. Mean reference allele frequency (MRAF)
#' will be also included if user provided the reference allele frequencies in
#' \code{mrafFilePath}.}
#' }
#'
#' @section Intermediate Results:
#'
#' The \code{Details/} folder contains all intermediate results including
#'
#' \describe{
#' \item{\code{SummaryGetI.txt}}{This table contains mean chi-square (square of
#' test statistics) of all SNPs before adjustment (RawMeanX2), Genomic control
#' value (GClambda), \code{getI} result (EstInf) and mean chi-square of all SNPs
#' after \code{getI} based adjustment (AdjMeanX2) for each GWAS
#' (\code{traitName}).}
#' \item{\code{SummaryGetPsi.txt}}{This table contains Pearson correlation
#' coefficient of test statistics from all SNPs (StatCor), \code{getPsi} result
#' (\code{Psi}) and \code{getPi} result (considering all SNPs, allPi; considering only
#' significant SNPs, sigPi) for all GWAS pairs (GWAS1 and GWAS2).}
#' \item{\code{EbICoW.txt}}{This table contains effect of EbICoW GWAS (define as
#' mean chi-square of all SNPs minus one, AllEff), Equivalent sample size (Ess),
#' EbICoW GWAS combination members (GWAS*), corresponding weight of effect size
#' in combination (bw*), corresponding weight of test statistics in combination
#' (tw*) and rounds of all iterations involved in combine (r*) for each EbICoW
#' GWAS (\code{traitName}).}
#' \item{\code{SummaryEbICoWGetPsi.txt}}{This table contains Pearson correlation
#' coefficient of test statistics from all SNPs (StatCor), \code{getPsi} result
#' (\code{Psi}) and \code{getPi} result (considering all SNPs, allPi; considering only
#' significant SNPs, sigPi) for all EbICoW GWAS pairs (GWAS1 and GWAS2).}
#' \item{\code{Summaryi-EbICoW.txt}}{This table contains the details of each
#' \code{i-EbICoW} iteration round.
#' 1) Iteration round count (Round);
#' 2) Two GWASs used in the iteration (GWAS1 and GWAS2);
#' 3) Name of combined GWAS (NewGWAS);
#' 4) Type of \code{Pi} and h chosen by optimize function (Piall and hall, AEC; Pisig
#' and hsig, SEC; Pistb and hstb, STB);
#' 5) Whether pass evaluation of evaluate function (Evaluate);
#' 6) Effect of GWAS1, GWAS2 and NewGWAS (define as mean chi-square of all SNPs
#' minus one, Eff1X2, Eff2X2 and NewEffX2);
#' 7) \code{getPsi} and \code{getPi} result (\code{Psi}; considering all SNPs, allPi;
#' considering only significant SNPs, sigPi);
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
#' and GWAS1-GWAS2 intersection (sigSNPSTBinterN);
#' 9) Equivalent sample size of GWAS1, GWAS2 and NewGWAS (Ess1, Ess2 and EssNew);
#' 10) Weight of effect size used in combining GWAS1 and GWAS2 chosen by
#' optimize function (bw1 and bw2);
#' 11) Weight of test statistics used in combining GWAS1 and GWAS2 chosen by
#' optimize function (tw1 and tw2).}
#' \item{\code{NullCGWAScorrection.txt}}{This dataset contains LOESS model
#' samples (ExpectedQuantile and ObservedP) simulated in \code{getCoef} function
#' using in adjusting C-GWAS result.}
#' \item{\code{NullCGWASdistribution.jpg}}{Simulated null distributions of the
#' unadjusted p-values (gray) and p-values corrected using the \code{getCoef}
#' function (orange) are compared with the uniform distribution (black) for
#' C-GWAS simulation.}
#' \item{\code{NullMinpcorrection.txt}}{This dataset contains LOESS model
#' samples (ExpectedQuantile and ObservedP) simulated in \code{getCoef} function
#' using in adjusting MinGWAS result.}
#' \item{\code{NullMinpdistribution.jpg}}{Simulated null distributions of the
#' unadjusted p-values (gray) and p-values corrected using the \code{getCoef}
#' function (orange) are compared with the uniform distribution (black) for
#' MinGWAS simulation.}
#' \item{\code{Details/i-EbICoW}}{If user choose to keep \code{i-EbICoW} output
#' (\code{keepIEb = TRUE}), all summary statistics (effect sizes, *.beta and
#' test statistics, *.stat) of EbICoW GWASs are saved in \code{Details/i-EbICoW}.}
#' }
#'
#' @examples
#' # example of the input GWAS file
#' f1 <- read.table(file.path(system.file("extdata", package = "CGWAS"),
#' 'Y1.assoc'), header=TRUE)
#' head(f1)
#'
#' # example of the input mean reference allele frequency file
#' f2 <- read.table(file.path(system.file("extdata", package = "CGWAS"),
#' 'MRAF'), header=TRUE)
#' head(f2)
#'
#' # example of the input SNP information file
#' f3 <- read.table(file.path(system.file("extdata", package = "CGWAS"),
#' 'SnpInfo'), header=TRUE)
#' head(f3)
#'
#' # demo of whole C-GWAS procedure implementing
#' # the output files are in the current working directory
#' ExDataDir <- system.file("extdata", package = "CGWAS")
#' gwasFileName <- c("Y1.assoc", "Y2.assoc", "Y3.assoc",
#'                   "Y4.assoc", "Y5.assoc", "Y6.assoc",
#'                   "Y7.assoc", "Y8.assoc", "Y9.assoc",
#'                   "Y10.assoc", "Y11.assoc", "Y12.assoc")
#' gwasFilePath <- file.path(ExDataDir, gwasFileName)
#' snpFilePath <- file.path(ExDataDir, 'SnpInfo')
#' outputPath <- getwd()
#' traitName <- c("Y1", "Y2", "Y3",
#'                "Y4", "Y5", "Y6",
#'                "Y7", "Y8", "Y9",
#'                "Y10", "Y11L", "Y12")
#' mrafFilePath <- file.path(ExDataDir, 'MRAF')
#' indSNPN = 1e5
#'
#' cgwas(gwasFilePath, snpFilePath, outputPath,
#'       traitName = traitName, mrafFilePath = mrafFilePath, indSNPN = indSNPN)
#'
#' @section Future Development Plan:
#' The current version is 0.9.3, which implementing the whole C-GWAS procedure
#' into a single function \code{cgwas}. Standing along functions facilitating
#' the C-GWAS analysis will be implemented in a later version, e.g., \code{getI},
#' \code{getPsi}, \code{getPi}.
#'
#' @seealso
#' For more instruction of running CGWAS, please see
#' \url{https://github.com/FanLiuLab/CGWAS}
#'
#' @export
cgwas <- function(gwasFilePath,
                  snpFilePath,
                  outputPath,
                  traitName = basename(gwasFilePath),
                  exNa = TRUE,
                  mrafFilePath = NULL,
                  keepIEb = FALSE,
                  threadN = ceiling(parallel::detectCores()/2),
                  indSNPN = 1e6,
                  simulDep = 100,
                  pStudy = 0.05/indSNPN,
                  pSugst = 1/indSNPN,
                  pMain = 3/indSNPN,
                  ebicoPwrInc = 1.1,
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
\\____/       \\____/   |__/|__//_/  |_|/____/      version 0.9.3\n\n")

  cgwasenv <- CGWAS_ENV$new()

  # A string list of path to GWASs summary data files.
  # expose
  cgwasenv$.GWAS_FILE_PATH <- gwasFilePath

  # A string list of path to GWASs summary data files.
  # expose
  cgwasenv$.SNP_FILE_PATH <- snpFilePath

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
  cgwasenv$.MRAF_FILE_PATH <- mrafFilePath

  # whether MRAF file exist.
  # hide
  cgwasenv$.MRAF_FILE_EXIST <- !is.null(mrafFilePath)

  # Results output directory.
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
  # integer, >=1e4, default 1e6
  # roundup to nearest 500
  # expose
  if (indSNPN < 1e4) {
    return(message("ERROR: parament indSNP invalid!"))
  } else {
    cgwasenv$.IND_SNP_N <- ceiling(indSNPN / 500) * 500
  }

  # Quantile simulation depth
  # integer, >=50, default 100
  # roundup to nearest 50
  # expose
  if (simulDep < 50) {
    return(message("ERROR: parament simulDep invalid!"))
  } else {
    cgwasenv$.SIMUL_DEP <- ceiling(simulDep / 50) * 50
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
  # >0, default 1.1
  # expose
  if (ebicoPwrInc <= 0) {
    return(message("ERROR: parament ebicoPwrInc invalid!"))
  } else {
    cgwasenv$.MIN_EbICo_POWER_INC <- ebicoPwrInc
  }

  # min correlation difference
  # >=0 & <2, default 0.05
  # expose
  if (minCorDiff < 0 & minCorDiff >= 2) {
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
  cgwasenv$.CGWAS_RESULT_PATH <- file.path(cgwasenv$.CGWAS_DIR, "Results")
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
