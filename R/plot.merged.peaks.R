
#' Title
#'
#' @param plot.startvec 
#' @param plot.fit.frame 
#' @param plot.bin.counts 
#' @param plot.dp 
#' @param PARAMETERS 
#'
#' @return
#' @export
#'
#' @examples
.plot.merged.peaks = function(
  plot.startvec, 
  plot.fit.frame, 
  plot.bin.counts, 
  plot.dp,
  plot.merged.peaks,
  PARAMETERS
){
  
  merged.peaks.df = data.frame(plot.merged.peaks, stringsAsFactors = F)
  
  filename = paste0(PARAMETERS$OUTPUTDIR, "/", PARAMETERS$GENE, ".MergedPeaks.pdf")
  pdf(filename)
  
  p1 = ggplot() +
    geom_histogram(data=plot.startvec, aes(x=start), binwidth = 50, colour = "black", fill="white") +
    geom_line(data=plot.fit.frame, aes(x=x,y=y), colour='red') +
    theme_bw() +
    ggtitle(PARAMETERS$GENE) +
    annotate("rect", xmin=merged.peaks.df$start, xmax=merged.peaks.df$end, ymin=-1 , ymax=-0.1, alpha=0.2, color="black", fill=rainbow(nrow(merged.peaks.df)))
  print(p1)
  
  p2 = ggplot() +
    geom_line(data=plot.bin.counts, aes(x=start,y=Coverage), colour='blue') +
    geom_line(data=plot.fit.frame, aes(x=x,y=y), colour='red') +
    theme_bw() +
    ggtitle(PARAMETERS$GENE) +
    annotate("rect", xmin=merged.peaks.df$start, xmax=merged.peaks.df$end, ymin=-1 , ymax=0.1, alpha=0.2, color="black", fill=rainbow(nrow(merged.peaks.df)))
  print(p2)
  
  dev.off()
  
}