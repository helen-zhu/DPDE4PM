
#' Title
#'
#' @param PEAKS 
#' @param GENEINFO 
#'
#' @return
#' @export
#'
#' @examples
.retrieve.peaks.as.granges = function(PEAKS, GENEINFO){
  
  PEAKS = PEAKS[PEAKS$name == GENEINFO$gene,]

  # Verify that all the columns have the
  BED6 = data.frame(stringsAsFactors = F)
  for(i in 1:nrow(PEAKS)){
    tmp = data.frame(
      "chr" = PEAKS$chr[i],
      "start" = PEAKS$start[i] + strtoi(strsplit(PEAKS$blockStarts[i], split = ",")[[1]]),
      "end" = PEAKS$start[i] + strtoi(strsplit(PEAKS$blockStarts[i], split = ",")[[1]]) + strtoi(strsplit(PEAKS$blockSizes[i], split = ",")[[1]]),
      "name" = PEAKS$name[i],
      "score" = PEAKS$score[i],
      "strand" = PEAKS$strand[i],
      "peak" = i,
      "sample" = PEAKS$sample[i],
      stringsAsFactors = F
    )
    BED6 = rbind(BED6, tmp)
  }
  
  makeGRangesFromDataFrame(BED6, keep.extra.columns = T)
}