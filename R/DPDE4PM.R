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
#' @param WRITE.OUTPUT TRUE or FALSE, whether an output file should be saved
#' @param OUTPUT.TAG A character string indicating a tag to track the generated files
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
  GTF,
  ANNOTATION,
  RESOLUTION = 50,
  DP.ITERATIONS = 1000,
  WEIGHT.THRESHOLD = 0.2,
  N.SD = 1,
  OUTPUTDIR =".",
  PLOT.RESULT=F,
  WRITE.OUTPUT=T,
  OUTPUT.TAG="",
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
  PARAMETERS$WRITE.OUTPUT = WRITE.OUTPUT
  PARAMETERS$OUTPUT.TAG = OUTPUT.TAG
  PARAMETERS$ALPHA.PRIORS = ALPHA.PRIORS
  PARAMETERS$SEED = SEED

  # Error messages
  # Check if OUTPUTDIR exists
  if(!dir.exists(PARAMETERS$OUTPUTDIR)){
    dir.create(PARAMETERS$OUTPUTDIR, recursive = T)
  }
  # 2. Check if peaks are in the right format
  # 3. Check if gene is in peaks & gtf

  # Making a vector of all samples so that the output is consistent
  ALL.SAMPLES = sort(unique(PEAKS$sample))
  PARAMETERS$ALL.SAMPLES = ALL.SAMPLES

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

  # Creating some big peaks
  REDUCED.GENE.PEAKS.GR = reduce(GENEPEAKSGR)

  # Generating Data
  TILED.PEAKS.GR = GenomicRanges::tile(GENEPEAKSGR, width = PARAMETERS$RESOLUTION)
  startvec = data.frame(TILED.PEAKS.GR, stringsAsFactors = F)

  # Dirichlet Process
  startvec.mean = mean(as.vector(startvec$start))
  startvec.sd = sd(as.vector(startvec$start))
  startvec.sd = ifelse( is.na(startvec.sd) | startvec.sd == 0, 1, startvec.sd) # This is for when there's 1 short peak or all samples have the exact same peak
  startvec.scaled = (as.vector(startvec$start) - startvec.mean)/startvec.sd

  set.seed(PARAMETERS$SEED)
  dp = dirichletprocess::DirichletProcessGaussian(
    y = startvec.scaled,
    g0Priors = c(0, 1, 1, 1),
    alphaPriors = PARAMETERS$ALPHA.PRIORS)

  dp = dirichletprocess::Fit(
    dpObj = dp,
    its = PARAMETERS$DP.ITERATIONS,
    updatePrior = F,
    progressBar = TRUE)

  # Examining Weights & Distributions of the GMM
  dp_data = data.frame(
    "i" = 1,
    "Weights" = dp$weights,
    "Mu" = c(dp$clusterParameters[[1]]),
    "Sigma" = c(dp$clusterParameters[[2]]),
    stringsAsFactors = F)
  dp_data$Mu = (dp_data$Mu*startvec.sd)+startvec.mean
  dp_data$Sigma = (dp_data$Sigma*startvec.sd)

  # Generating Peaks
  merged.peaks = .generate.peaks.from.gmm(dp_data, PARAMETERS, GENEINFO)
  merged.peaks.rna = merged.peaks[[1]]
  merged.peaks.genome = merged.peaks[[2]]

  # Peak Coverage
  PEAK.COVERAGE = GenomicRanges::coverage(GENEPEAKSGR)
  BINS = unlist(GenomicRanges::tile(REDUCED.GENE.PEAKS.GR, width = 1))
  BIN.COUNTS = data.frame(GenomicRanges::binnedAverage(BINS, PEAK.COVERAGE, "Coverage"), stringsAsFactors = F)

  # Plotting Data
  x.norm <- seq(min(startvec.scaled), max(startvec.scaled), by=0.01)
  y.fit <- data.frame(replicate(100, dirichletprocess::PosteriorFunction(dp)(x.norm)))
  fit.frame <- data.frame(x=x.norm, y=rowMeans(y.fit))
  fit.frame$x = (fit.frame$x*startvec.sd)+startvec.mean
  fit.frame$y = fit.frame$y*max(BIN.COUNTS$Coverage)

  # Plotting DP data
  sample.points = seq(1, GENEINFO$exome_length, 10)
  scaling.factor = length(unique(PEAKSGR$sample))*100
  plot.dp.data = data.frame("sample.points" = sample.points, stringsAsFactors = F)
  for(i in 1:nrow(dp_data)){
    norm.tmp = dnorm(sample.points, mean = dp_data$Mu[i], sd = dp_data$Sigma[i])*dp_data$Weights[i]*scaling.factor
    plot.dp.data = cbind(plot.dp.data, norm.tmp)
  }
  colnames(plot.dp.data) = c("sample.points", paste0("V", 1:nrow(dp_data)))

  # Plotting
  if(PARAMETERS$PLOT.RESULT){
    .plot.merged.peaks(startvec, fit.frame, BIN.COUNTS, merged.peaks.rna, plot.dp.data, PARAMETERS)
  }

  # Return a Data Frame of Merged Peaks
  if(length(merged.peaks.genome) == 0){
    output.ncol = length(PARAMETERS$ALL.SAMPLES) + 12
    OUTPUT_TABLE = data.frame(matrix(ncol = output.ncol, nrow = 0))
    names(OUTPUT_TABLE) = c("chr", "start", "end", "name", "score", "strand", "thickStart", "thickEnd",
                            "itemRgb", "blockCount", "blockSizes","blockStarts", PARAMETERS$ALL.SAMPLES)
  } else {

    # Creating a BED12 File
    GenomicRanges::start(merged.peaks.genome) = GenomicRanges::start(merged.peaks.genome)-1
    PEAKS.FINAL = .bed6tobed12(MERGED.PEAKS = merged.peaks.genome, ID.COLS = c("name", "i", "j"))

    # Merging P-Values
    SAMPLE.PVAL = .merge.p(PEAKSGR, MERGED.PEAKS = merged.peaks.genome, ANNOTATION, PARAMETERS, ID.COLS = c("name", "i", "j"))

    # Write Output Tables & Return Files
    OUTPUT.TABLE = merge(PEAKS.FINAL, SAMPLE.PVAL, by = "peak", all = T)
    OUTPUT.TABLE = OUTPUT.TABLE[,colnames(OUTPUT.TABLE) != "peak"]
  }

  if(PARAMETERS$WRITE.OUTPUT){
    filename = paste0(PARAMETERS$OUTPUTDIR, "/", PARAMETERS$GENE, ".", PARAMETERS$OUTPUT.TAG, ".MergedPeaks.tsv")
    write.table(
      OUTPUT.TABLE,
      file = filename,
      sep = "\t",
      col.names = T,
      row.names = F,
      quote = F
    )
  }

  return(OUTPUT.TABLE)

}
