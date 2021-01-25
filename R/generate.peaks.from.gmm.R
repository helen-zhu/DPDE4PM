
#' Title
#'
#' @param dp_data 
#' @param PARAMETERS 
#' @param GENEINFO 
#'
#' @return
#' @export
#'
#' @examples
.generate.peaks.from.gmm = function(
  dp_data,
  PARAMETERS,
  GENEINFO
){
  
  # Filtering by Threshold
  dp_data = dp_data[dp_data$Weights > PARAMETERS$WEIGHT.THRESHOLD,]
  if(nrow(dp_data) == 0){
    warning("No Peaks Survive Past The Threshold. Consider Lowering Requirements or Tuning Parameters.",
            call. = TRUE, domain = NULL)
    return(list(data.frame(), data.frame()))
  }
  dp_data = dp_data[order(-dp_data$Weights),]
  
  # Creating Peaks
  merged.peaks = data.frame(
    "chr" = GENEINFO$chr,
    "start" = round(dp_data$Mu - PARAMETERS$N.SD*dp_data$Sigma),
    "end" = round(dp_data$Mu + PARAMETERS$N.SD*dp_data$Sigma),
    "name" = GENEINFO$gene,
    "strand" = GENEINFO$strand,
    "weights" = dp_data$Weights,
    "i" = dp_data$i,
    "j" = seq(1, nrow(dp_data), 1),
    stringsAsFactors = F
  )
  
  # Filtering Overlapping Peaks
  merged.peaks.gr =  makeGRangesFromDataFrame(merged.peaks, keep.extra.columns = T)
  merged.peaks.filtered.rna = do.call(c, lapply(1:length(merged.peaks.gr), function(i){
    if(i == 1){
      tmp = merged.peaks.gr[1]
    }else{
      tmp = setdiff(merged.peaks.gr[i], merged.peaks.gr[c(1:(i-1))])
      if(length(tmp) > 0){mcols(tmp) = mcols(merged.peaks.gr[i])}
    }
    tmp
  }))
  
  # Transferring to Genomic Coordinates
  merged.peaks.genome = merged.peaks.filtered.rna
  end(merged.peaks.genome) = GENEINFO$RNA2DNA[end(merged.peaks.genome)]
  start(merged.peaks.genome) = GENEINFO$RNA2DNA[start(merged.peaks.genome)]
  anno_gr = makeGRangesFromDataFrame(GENEINFO$anno)
  
  # Filtering Out Introns
  merged.peaks.filtered.genome = do.call(c, lapply(1:length(merged.peaks.genome), function(i){
      tmp = intersect(merged.peaks.genome[c(i)], anno_gr)
      if(length(tmp) > 0){mcols(tmp) = mcols(merged.peaks.genome[i])}
      tmp
  }))

  return(list(merged.peaks.filtered.rna, merged.peaks.filtered.genome))
}