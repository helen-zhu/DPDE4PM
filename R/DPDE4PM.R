#' DPDE4PM merged RNA-seq peaks called on different samples
#'
#' @param GENE A 'character' gene id corresponding to the gene_id's found in the GTF files and the 'name' column in the peak files
#' @param PEAKS A data frame containing the following columns, and potentially extras, usually found in a BED12 file, base 0 system
#' \describe{
#'   \item{chr}{chromosomes, same as in GTF file}
#'   \item{start}{starting position of the peak, base 0}
#'   \item{end}{end position of the peak, base 0}
#'   \item{name}{gene id}
#'   \item{score}{p-value associated with the peak}
#'   \item{strand}{strand of the gene}
#'   \item{blockCount}{number of segments in the peak}
#'   \item{blockSizes}{size of segments in the peak, BED12 notation}
#'   \item{blockStarts}{starting positions of segments, BED12 notation}
#'   \item{sample}{sample_id of samples}
#' }
#' @param GTF The GTF file used to generate the peaks. This is used to determine the genomic coordinates of the gene.
#' @param ANNOTATION A object created by the read.gtf function. This allows the user to provide an annotation object to save compute time during parallelization.
#' @param RESOLUTION The width (bps) used to sample points from the peaks. This is likely optimized by choosing the window size used to generate the peaks.
#' @param DP.ITERATIONS Number of iterations used to fit the Dirichlet Process
#' @param OUTPUTDIR Output directory
#' @param WEIGHT.THRESHOLD A proportion (out of 1) used to determine the percentage of data accounted for by the GMM to be called a peak
#' @param N.SD Number of standard deviations from the mean for each fitted Gaussian that should considered part of the joint peak
#' @param PLOT.RESULT TRUE or FALSE, whether plots should be generated
#' @param ALPHA.PRIORS A length 2 numeric vector indicating alpha, beta of the Gamma distribution from which the concentration weight parameter alpha is drawn, see the R package dirichletprocess
#' @param SEED a seed for reproducibility
#'
#' @return A dataframe with BED12 columns and additional columns for each sample in the PEAKS dataframe and associated p-value with that peak.
#'
#' @importFrom grDevices dev.off pdf rainbow
#' @importFrom graphics layout mtext par
#' @importFrom stats end sd start
#' @importFrom utils write.table
#'
#' @export DPDE4PM
#'
#' @examples
DPDE4PM = function(
  GENE,
  PEAKS,
  GTF = NULL,
  ANNOTATION = NULL,
  RESOLUTION = 50,
  DP.ITERATIONS = 1000,
  WEIGHT.THRESHOLD = 0.2,
  N.SD = 1,
  OUTPUTDIR =".",
  PLOT.RESULT=F,
  ALPHA.PRIORS = c(1,2),
  SEED = 123
){

  # Making a list of parameters to pass back and forth
  PARAMETERS = list()
  PARAMETERS$GENE = GENE
  PARAMETERS$GTF = GTF
  PARAMETERS$RESOLUTION = RESOLUTION
  PARAMETERS$DP.ITERATIONS = DP.ITERATIONS
  PARAMETERS$WEIGHT.THRESHOLD = WEIGHT.THRESHOLD
  PARAMETERS$N.SD = N.SD
  PARAMETERS$OUTPUTDIR = OUTPUTDIR
  PARAMETERS$PLOT.RESULT = PLOT.RESULT
  PARAMETERS$ALPHA.PRIORS = ALPHA.PRIORS
  PARAMETERS$SEED = SEED

  # Check if peaks are in the right format

  # Creating a list of samples
  ALL.SAMPLES = sort(unique(PEAKS$sample))
  PARAMETERS$ALL.SAMPLES = ALL.SAMPLES

  # Check if OUTPUTDIR exists
  if(!dir.exists(PARAMETERS$OUTPUTDIR)){
    dir.create(PARAMETERS$OUTPUTDIR, recursive = T)
  }

  # If the gene doesn't have peaks
  if(PARAMETERS$GENE %in% PEAKS$name){
    warning("No Peaks are Found for This Gene in PEAKS!", call. = TRUE, domain = NULL)
    return(.generate.null.result(PARAMETERS))
  }

  # Turning PEAKS into a GRanges Object
  PEAKSGR = .retrieve.peaks.as.granges(PEAKS = PEAKS, GENE = PARAMETERS$GENE, DF = F)

  # Import GTF as a GRanges Object
  annot.format = F
  if(!is.null(ANNOTATION)){
    annot.format = .check.annotation(ANNOTATION, PEAKSGR, GENE = PARAMETERS$GENE)
  }
  if(!annot.format){
    ANNOTATION = read.gtf(PARAMETERS)
  }

  # Get Gene Information
  GENEINFO = .get.gene.anno(PARAMETERS, ANNOTATION)

  # Converting to RNA
  GENEPEAKSGR = GenomicRanges::shift(PEAKSGR, -1*GENEINFO$left+1)
  GenomicRanges::start(GENEPEAKSGR) = GENEINFO$DNA2RNA[GenomicRanges::start(GENEPEAKSGR)+1]
  GenomicRanges::end(GENEPEAKSGR) = GENEINFO$DNA2RNA[GenomicRanges::end(GENEPEAKSGR)]

  # Reduce Overlapping Peaks in the Same Sample
  GENEPEAKSGR = S4Vectors::split(GENEPEAKSGR, GENEPEAKSGR$sample)
  GENEPEAKSGR = unlist(GenomicRanges::reduce(GENEPEAKSGR))

  # Examining Weights & Distributions of the GMM
  dp = .dpc.peaks(GENEPEAKSGR, PARAMETERS)

  # Generating Peaks
  merged.peaks.rna = .generate.peaks.from.gmm(dp = dp, PARAMETERS, GENEINFO)
  merged.peaks.genome = .rna.peaks.to.genome(merged.peaks.rna, GENEINFO)

  # Plotting
  if(PARAMETERS$PLOT.RESULT){
    plotting.data = .generate.merged.peaks.plotting(dp, GENEINFO, GENEPEAKSGR)
    .plot.merged.peaks(plotting.data, merged.peaks.rna, PARAMETERS)
  }

  # Return a Data Frame of Merged Peaks
  if(length(merged.peaks.genome) == 0){
    warning("No Peaks are Found for This Gene After Fitting!", call. = TRUE, domain = NULL)
    OUTPUT.TABLE = .generate.null.result(PARAMETERS)
  } else {

    # Creating a BED12 File
    GenomicRanges::start(merged.peaks.genome) = GenomicRanges::start(merged.peaks.genome)-1
    PEAKS.FINAL = .bed6tobed12(MERGED.PEAKS = merged.peaks.genome, ID.COLS = c("name", "i"))

    # Merging P-Values
    SAMPLE.PVAL = .merge.p(PEAKSGR, MERGED.PEAKS = merged.peaks.genome, ANNOTATION, PARAMETERS, ID.COLS = c("name", "i"))

    # Write Output Tables & Return Files
    OUTPUT.TABLE = merge(PEAKS.FINAL, SAMPLE.PVAL, by = "peak", all = T)
    OUTPUT.TABLE = OUTPUT.TABLE[,colnames(OUTPUT.TABLE) != "peak"]
  }

  return(OUTPUT.TABLE)

}
