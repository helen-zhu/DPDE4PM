
#' Makes a plot of the merged peaks and the fitted curve
#'
#' @param plot.startvec A data frame of the points surveyed
#' @param plot.fit.frame A data frame with the fitted Dirichlet process model on the points in the startvec
#' @param plot.bin.counts A data frame of the coverage at a base pair resolution
#' @param plot.merged.peaks A GRanges object containing the RNA coordinates of the merged peaks in BED6 format
#' @param plot.dp.data A data frame of individual normal distributions fitted
#' @param PARAMETERS A PARAMETERS list with the parameters indicated in the DPDE4PM function
.plot.merged.peaks = function(
  plot.startvec,
  plot.fit.frame,
  plot.bin.counts,
  plot.merged.peaks,
  plot.dp.data,
  PARAMETERS
){

  merged.peaks.df = data.frame(plot.merged.peaks, stringsAsFactors = F)
  melted.plot.dp.data = reshape2::melt(plot.dp.data, id = "sample.points")

  filename = paste0(PARAMETERS$OUTPUTDIR, "/", PARAMETERS$GENE, ".", PARAMETERS$OUTPUT.TAG, ".MergedPeaks.pdf")
  pdf(filename)

  p1 = ggplot2::ggplot() +
    ggplot2::geom_histogram(data=plot.startvec, ggplot2::aes(x=start), binwidth = 50, colour = "black", fill="white") +
    ggplot2::geom_line(data=plot.fit.frame, ggplot2::aes(x=x, y=y), colour='red') +
    ggplot2::theme_bw() +
    ggplot2::ggtitle(PARAMETERS$GENE) +
    ggplot2::ylab("Binned Counts From Peaks (Used to Fit GMM)") +
    ggplot2::xlab("Transcript Coordinate") +
    ggplot2::annotate("rect", xmin=merged.peaks.df$start, xmax=merged.peaks.df$end, ymin=-1 , ymax=-0.1, alpha=0.5, color="black", fill=rainbow(nrow(merged.peaks.df)))

  p1 = p1 + geom_line(data = melted.plot.dp.data, ggplot2::aes(x=sample.points, y=value, color = variable)) +
    theme(legend.position = "none")

  print(p1)

  p2 = ggplot2::ggplot() +
    ggplot2::geom_line(data=plot.bin.counts, ggplot2::aes(x=start, y=Coverage), colour='blue') +
    ggplot2::geom_line(data=plot.fit.frame, ggplot2::aes(x=x,y=y), colour='red') +
    ggplot2::theme_bw() +
    ggplot2::ggtitle(PARAMETERS$GENE) +
    ggplot2::ylab("Coverage (at BP resolution)") +
    ggplot2::xlab("Transcript Coordinate") +
    ggplot2::annotate("rect", xmin=merged.peaks.df$start, xmax=merged.peaks.df$end, ymin=-1 , ymax=-0.1, alpha=0.5, color="black", fill=rainbow(nrow(merged.peaks.df)))

  p2 = p2 + geom_line(data = melted.plot.dp.data, ggplot2::aes(x=sample.points, y=value, color = variable)) +
    theme(legend.position = "none")

  print(p2)

  dev.off()

}
