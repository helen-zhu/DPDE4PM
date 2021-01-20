
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
  PARAMETERS
){
  
  filename = paste0(PARAMETERS$OUTPUTDIR, "/", PARAMETERS$GENE, ".MergedPeaks.pdf")
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
  
}