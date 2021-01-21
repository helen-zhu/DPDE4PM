
#' Title
#'
#' @param PEAKSGR
#' @param MERGED.PEAKS
#' @param ANNOTATION
#' @param PARAMETERS
#' @param ID.COLS
#'
#' @return
#' @export
#'
#' @examples
merge.p = function(PEAKSGR, MERGED.PEAKS, ANNOTATION, PARAMETERS, ID.COLS){

  # Numeric Stability
  PEAKSGR$score = ifelse(PEAKSGR$score == 0, 2e-16, PEAKSGR$score)

  # Tag
  MERGED.PEAKS$tag = apply(mcols(MERGED.PEAKS)[,ID.COLS], 1, function(x) paste(x, collapse = ":"))

  # Finding Overlaps
  overlaps = findOverlaps(MERGED.PEAKS, PEAKSGR)
  overlapping_peaks = PEAKSGR[subjectHits(overlaps)]
  overlapping_peaks$peak = MERGED.PEAKS$tag[queryHits(overlaps)]
  overlapping_peaks$tag = paste0(overlapping_peaks$sample, ":", MERGED.PEAKS$tag[queryHits(overlaps)])

  overlapping_peaks = split(overlapping_peaks, overlapping_peaks$tag)
  results = do.call(rbind, lapply(1:length(overlapping_peaks), function(i){
    tmp = overlapping_peaks[[i]]
    p = ifelse(length(tmp) > 1, metap::sumlog(tmp$score)$p, tmp$score)
    data.frame(
      "peak" = tmp$peak[1],
      "sample" = tmp$sample[1],
      "p" = p,
      stringsAsFactors = F
    )
  }))

  # Format this into something usable
  pmat = reshape2::dcast(results, peak ~ sample, value.var = "p")
  pmat[is.na(pmat)] <- 1
  return(pmat)

}
