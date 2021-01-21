#' Title
#'
#' @param GENE
#' @param PEAKS
#' @param GTF
#' @param RESOLUTION
#' @param DP.ITERATIONS
#' @param OUTPUTDIR
#'
#' @return
#' @export
#'
#' @examples
DPDE4PM = function(
  GENE,
  PEAKS,
  GTF,
  RESOLUTION = 50,
  DP.ITERATIONS = 1000,
  WEIGHT.THRESHOLD = 0.2,
  N.SD = 1,
  OUTPUTDIR =".",
  PLOT.RESULT=F
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

  # Error messages
  # 1. Check if OUTPUTDIR exists
  # 2. Check if peaks are in the right format
  # 3. Check if gene is in peaks & gtf

  # Import GTF as a GRanges Object
  ANNOTATION = .read.gtf(PARAMETERS)

  # Get Gene Information
  GENEINFO = .get.gene.anno(PARAMETERS, ANNOTATION)

  # Turning PEAKS into a GRanges Object
  PEAKSGR = makeGRangesFromDataFrame(PEAKS, keep.extra.columns = T)

  # Split into peaks
  GENEPEAKSGR = PEAKSGR[PEAKSGR$name == PARAMETERS$GENE]

  # Converting to RNA
  GENEPEAKSGR = shift(GENEPEAKSGR, -1*GENEINFO$left+1)
  start(GENEPEAKSGR) = GENEINFO$DNA2RNA[start(GENEPEAKSGR)+1]
  end(GENEPEAKSGR) = GENEINFO$DNA2RNA[end(GENEPEAKSGR)]

  # Reduce Overlapping Peaks in the Same Sample
  GENEPEAKSGR = split(GENEPEAKSGR, GENEPEAKSGR$sample)
  GENEPEAKSGR = unlist(reduce(GENEPEAKSGR))

  # Creating some big peaks
  REDUCED.GENE.PEAKS.GR = reduce(GENEPEAKSGR)

  # Initializing Data
  plot.startvec = data.frame(stringsAsFactors = F)
  plot.fit.frame = data.frame(stringsAsFactors = F)
  plot.bin.counts = data.frame(stringsAsFactors = F)
  plot.dp = data.frame(stringsAsFactors = F)
  plot.merged.peaks = GRanges()
  merged.peaks.genome = GRanges()

  for(i in 1:length(REDUCED.GENE.PEAKS.GR)){

    # Generating Data
    OVERLAP.PEAKS.GR = GENEPEAKSGR[queryHits(findOverlaps(GENEPEAKSGR, REDUCED.GENE.PEAKS.GR[i]))]
    TILED.PEAKS.GR = tile(OVERLAP.PEAKS.GR, width = PARAMETERS$RESOLUTION)
    startvec = data.frame(TILED.PEAKS.GR, stringsAsFactors = F)

    # Dirichlet Process
    startvec.mean = mean(as.vector(startvec$start))
    startvec.sd = sd(as.vector(startvec$start))
    startvec.scaled = (as.vector(startvec$start) - startvec.mean)/startvec.sd

    dp = DirichletProcessGaussian(
      y = startvec.scaled,
      g0Priors = c(0, 1, 1, 1),
      alphaPriors = c(2, 4))

    dp = Fit(dpObj = dp,
             its = PARAMETERS$DP.ITERATIONS,
             updatePrior = F,
             progressBar = TRUE)

    # Examining Weights & Distributions of the GMM
    dp_data = data.frame(
      "i" = i,
      "Weights" = dp$weights,
      "Mu" = c(dp$clusterParameters[[1]]),
      "Sigma" = c(dp$clusterParameters[[2]]),
      stringsAsFactors = F
    )
    dp_data$Mu = (dp_data$Mu*startvec.sd)+startvec.mean
    dp_data$Sigma = (dp_data$Sigma*startvec.sd)

    # Generating Peaks
    merged.peaks = .generate.peaks.from.gmm(dp_data, PARAMETERS, GENEINFO)
    merged.peaks.rna = merged.peaks[[1]]
    merged.peaks.gen = merged.peaks[[2]]

    # Peak Coverage
    PEAK.COVERAGE = coverage(OVERLAP.PEAKS.GR)
    BINS = tile(REDUCED.GENE.PEAKS.GR[i], width = 1)[[1]]
    BIN.COUNTS = data.frame(GenomicRanges::binnedAverage(BINS, PEAK.COVERAGE, "Coverage"), stringsAsFactors = F)

    # Plotting Data
    x.norm <- seq(min(startvec.scaled), max(startvec.scaled), by=0.01)
    y.fit <- data.frame(replicate(100, PosteriorFunction(dp)(x.norm)))
    fit.frame <- data.frame(x=x.norm, y=rowMeans(y.fit))
    fit.frame$x = (fit.frame$x*startvec.sd)+startvec.mean
    fit.frame$y = fit.frame$y*max(BIN.COUNTS$Coverage)

    # Updating Data
    plot.startvec = rbind(plot.startvec, startvec)
    plot.fit.frame = rbind(plot.fit.frame, fit.frame)
    plot.bin.counts = rbind(plot.bin.counts, BIN.COUNTS)
    plot.dp = rbind(plot.dp, dp_data)
    plot.merged.peaks = c(plot.merged.peaks, merged.peaks.rna)
    merged.peaks.genome = c(merged.peaks.genome, merged.peaks.gen)
  }

  # Plotting
  if(PARAMETERS$PLOT.RESULT){
    .plot.merged.peaks(plot.startvec, plot.fit.frame, plot.bin.counts, plot.dp, plot.merged.peaks, PARAMETERS)
  }

  # Merging P-Values
  # Find overlapping peaks with each peak and use Fisher's Method

  # Return a Data Frame of Merged Peaks
  PEAKS.FINAL = .bed12tobed6(merged.peaks.genome)

}
