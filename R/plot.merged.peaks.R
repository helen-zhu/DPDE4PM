
#' Makes a plot of the merged peaks and the fitted curve
#'
#' @param plot.startvec A data frame of the points surveyed
#' @param plot.fit.frame A data frame with the fitted Dirichlet process model on the points in the startvec
#' @param plot.bin.counts A data frame of the coverage at a base pair resolution
#' @param plot.merged.peaks A GRanges object containing the RNA coordinates of the merged peaks in BED6 format
#' @param PARAMETERS A PARAMETERS list with the parameters indicated in the DPDE4PM function
.plot.merged.peaks = function(
  plot.startvec,
  plot.fit.frame,
  plot.bin.counts,
  plot.merged.peaks,
  PARAMETERS
){

  merged.peaks.df = data.frame(plot.merged.peaks, stringsAsFactors = F)

  filename = paste0(PARAMETERS$OUTPUTDIR, "/", PARAMETERS$GENE, ".", PARAMETERS$OUTPUT.TAG, ".MergedPeaks.pdf")
  pdf(filename)

  p1 = ggplot() +
    geom_histogram(data=plot.startvec, aes(x=start), binwidth = 50, colour = "black", fill="white") +
    geom_line(data=plot.fit.frame, aes(x=x,y=y), colour='red') +
    theme_bw() +
    ggtitle(PARAMETERS$GENE) +
    ylab("Binned Counts From Peaks (Used to Fit GMM)") +
    xlab("Transcript Coordinate") +
    annotate("rect", xmin=merged.peaks.df$start, xmax=merged.peaks.df$end, ymin=-1 , ymax=-0.1, alpha=0.5, color="black", fill=rainbow(nrow(merged.peaks.df)))
  print(p1)

  p2 = ggplot() +
    geom_line(data=plot.bin.counts, aes(x=start,y=Coverage), colour='blue') +
    geom_line(data=plot.fit.frame, aes(x=x,y=y), colour='red') +
    theme_bw() +
    ggtitle(PARAMETERS$GENE) +
    ylab("Coverage (at BP resolution)") +
    xlab("Transcript Coordinate") +
    annotate("rect", xmin=merged.peaks.df$start, xmax=merged.peaks.df$end, ymin=-1 , ymax=-0.1, alpha=0.5, color="black", fill=rainbow(nrow(merged.peaks.df)))
  print(p2)

  dev.off()

}
