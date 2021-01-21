
#' Title
#'
#' @param PEAKS 
#' @param ID.COL 
#'
#' @return
#' @export
#'
#' @examples
.bed6tobed12 = function(
  MERGED.PEAKS,
  ID.COLS 
){
  
  MERGED.PEAKS$tag = apply(MERGED.PEAKS[,ID.COLS], 1, function(x) paste(x, collapse = ":"))
  
  BED12 = do.call(rbind, lapply(unique(MERGED.PEAKS$tag), function(itag) {
    TMP.PEAK = MERGED.PEAKS[MERGED.PEAKS$tag == itag,]
    TMP.PEAK = TMP.PEAK[order(TMP.PEAK$start),]
    bed12 = data.frame(
      "chr" = TMP.PEAK$chr[1],
      "start" = min(TMP.PEAK$start),
      "end" = max(TMP.PEAK$end),
      "name" = TMP.PEAK$name[1],
      "score" = 0, # This might be changed to a p-value if there is one
      "strand" = TMP.PEAK$strand[1],
      "thickStart" = min(TMP.PEAK$start),
      "thickEnd" = max(TMP.PEAK$end),
      "itemRgb" = 0,
      "blockCount" = nrow(TMP.PEAK),
      "blockSizes" = paste0(paste(TMP.PEAK$end - TMP.PEAK$start, collapse = ","), ","),
      "blockStarts" = paste0(paste(TMP.PEAK$start - min(TMP.PEAK$start), collapse = ","), ","),
      stringsAsFactors = F
      )
    
  }))
  return(BED12)
}
