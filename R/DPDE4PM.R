source("/cluster/home/helenzhu/code/snakemake_M6A_16_PeakMatTest/Rscripts/000_HEADER.R")

# Genomic Data Manipulation
library(GenomicFeatures)
library(GenomicRanges)

# Dirichlet Process
library(dirichletprocess)

# Plotting
library(ggplot2)
library(Sushi)

# Loading Test Data -------------------------------------------------------

all_samples = setdiff(all_samples, exclude_samples)

PEAKS = do.call(rbind, lapply(all_samples, function(i){
  cat(i, "\n")
  tmp = read.table(paste0("/cluster/home/helenzhu/Cluster_Helen/Snakemake_M6A_hg38/G_exomePeak_Peaks_SS/", i, "/peak.bed6"), header = F, sep = "\t", stringsAsFactors = F)
  tmp$V7 = i
  tmp
}))
colnames(PEAKS) = c("chr", "start", "end", "name", "score", "strand", "sample")

GTF = "/cluster/home/helenzhu/Cluster_Helen/Snakemake_M6A_hg38/GRCh38.p13_Build34/gencode.v34.chr_patch_hapl_scaff.annotation.gtf"

GENE = "ENSG00000129514.8"

RESOLUTION = 50
DP.ITERATIONS = 1000
OUTPUTDIR ="."

# Writing Functions -------------------------------------------------------

DPDE4PM = function(
  GENE,
  PEAKS,
  GTF,
  RESOLUTION = 50,
  DP.ITERATIONS = 1000,
  OUTPUTDIR ="."
){
  
  # Making a list of parameters to pass back and forth
  PARAMETERS = list()
  PARAMETERS$GENE = GENE
  # PARAMETERS$PEAKS = PEAKS
  PARAMETERS$GTF = GTF
  PARAMETERS$RESOLUTION = RESOLUTION
  PARAMETERS$DP.ITERATIONS = DP.ITERATIONS
  PARAMETERS$OUTPUTDIR = OUTPUTDIR

  # Import GTF as a GRanges Object
  ANNOTATION = .read.gtf(PARAMETERS)

  # Get DNA2RNA
  GENEINFO = .get.gene.anno(PARAMETERS, ANNOTATION)
  GENE.EXON.GR = GRanges(seqnames = GENEINFO$chr, IRanges(0, GENEINFO$exome_length), strand = GENEINFO$strand)
  
  # Turning PEAKS into a GRanges Object
  PEAKSGR = makeGRangesFromDataFrame(PEAKS, keep.extra.columns = T)
  
  # Split into peaks
  GENEPEAKSGR = PEAKSGR[PEAKSGR$name == PARAMETERS$GENE]
  # GENEPEAKSGR = GENEPEAKSGR[GENEPEAKSGR$sample == "CPCG0100"]
  
  # Converting to RNA
  GENEPEAKSGR = shift(GENEPEAKSGR, -1*GENEINFO$left+1)
  start(GENEPEAKSGR) = GENEINFO$DNA2RNA[start(GENEPEAKSGR)+1]
  end(GENEPEAKSGR) = GENEINFO$DNA2RNA[end(GENEPEAKSGR)]
  
  # Reduce Overlapping Peaks in the Same Sample
  GENEPEAKSGR = split(GENEPEAKSGR, GENEPEAKSGR$sample)
  GENEPEAKSGR = unlist(reduce(GENEPEAKSGR))

  # Creating some big peaks
  # REDUCED.GENE.PEAKS.GR = reduce(GENEPEAKSGR)
  
  # Initializing Data
  plot.startvec = data.frame(stringsAsFactors = F)
  plot.fit.frame = data.frame(stringsAsFactors = F)
  plot.bin.counts = data.frame(stringsAsFactors = F)
  plot.dp = data.frame(stringsAsFactors = F) # dp parameters
  
  # for(i in 1:length(REDUCED.GENE.PEAKS.GR)){
    
    # Generating Data
    # OVERLAP.PEAKS.GR = GENEPEAKSGR[queryHits(findOverlaps(GENEPEAKSGR, REDUCED.GENE.PEAKS.GR[i]))]
    # TILED.PEAKS.GR = tile(OVERLAP.PEAKS.GR, width = PARAMETERS$RESOLUTION)
    # startvec = data.frame(TILED.PEAKS.GR, stringsAsFactors = F)
    
    TILED.PEAKS.GR = tile(GENEPEAKSGR, width = PARAMETERS$RESOLUTION)
    startvec = data.frame(TILED.PEAKS.GR, stringsAsFactors = F)
    
    # Dirichlet Process
    startvec.mean = mean(as.vector(startvec$start))
    startvec.sd = sd(as.vector(startvec$start))
    startvec.scaled = (as.vector(startvec$start) - startvec.mean)/startvec.sd
    dp = DirichletProcessGaussian(startvec.scaled)
    dp = Fit(dp, PARAMETERS$DP.ITERATIONS) # progressBar = FALSE)
    
    # Peak Coverage
    # PEAK.COVERAGE = coverage(OVERLAP.PEAKS.GR)
    # BINS = tile(REDUCED.GENE.PEAKS.GR[i], width = 1)[[1]]
    # BIN.COUNTS = data.frame(GenomicRanges::binnedAverage(BINS, PEAK.COVERAGE, "Coverage"), stringsAsFactors = F)
    PEAK.COVERAGE = coverage(GENEPEAKSGR)
    BINS = unlist(tile(REDUCED.GENE.PEAKS.GR, width = 1))
    BIN.COUNTS = data.frame(GenomicRanges::binnedAverage(BINS, PEAK.COVERAGE, "Coverage"), stringsAsFactors = F)
    
    # Plotting Data
    x.norm <- seq(min(startvec.scaled)-1, max(startvec.scaled)+1, by=0.01)
    y.fit <- data.frame(replicate(100, PosteriorFunction(dp)(x.norm)))
    fit.frame <- data.frame(x=x.norm, y=rowMeans(y.fit))
    fit.frame$x = (fit.frame$x*startvec.sd)+startvec.mean
    fit.frame$y = fit.frame$y*max(BIN.COUNTS$Coverage)
    
    # Updating Data
    plot.startvec = rbind(plot.startvec, startvec)
    plot.fit.frame = rbind(plot.fit.frame, fit.frame)
    plot.bin.counts = rbind(plot.bin.counts, BIN.COUNTS)
    # plot.dp
  # }

  
  filename = "~/figures/testing.DP.pdf"
  pdf(filename)
  
  ggplot() + 
    geom_histogram(data=plot.startvec, aes(x=start), binwidth = 50, colour = "black", fill="white") + 
    geom_line(data=plot.fit.frame, aes(x=x,y=y), colour='red') +
    theme_bw() +
    ggtitle(PARAMETERS$GENE)
  
  ggplot() + 
    geom_line(data=plot.bin.counts, aes(x=start,y=Coverage), colour='blue') +
    geom_line(data=plot.fit.frame, aes(x=x,y=y), colour='red') +
    theme_bw() +
    ggtitle(PARAMETERS$GENE)
  
  dev.off()

  
  # Dirichlet Density Estimation
  
  # Plotting
  
  # Saving Result
  
  # Merging P-Values
  
  
}