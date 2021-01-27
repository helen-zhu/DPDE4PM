
#' Creates a BED12 file from a BED6 file using unique columns of ID.COLS
#'
#' @param MERGED.PEAKS A GRanges Object, in BED6 format for the merged peaks
#' @param ID.COLS Column names of columns whose unique values are used to merge peaks
#'
#' @return A BED12 file of merged peaks
#' \describe{
#'   \item{chr}{chromosomes, same as in GTF file}
#'   \item{start}{starting position of the peak, base 0}
#'   \item{end}{end position of the peak, base 0}
#'   \item{name}{gene id}
#'   \item{score}{p-value associated with the peak}
#'   \item{strand}{strand of the gene}
#'   \item{blockCount}{number of segments in the peak}
#'   \item{blockSizes}{size of segments in the peak, BED12 notation}
#'   \item{blockStarts}{starting positions of segments, BED12 notation}
#' }
#'
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
      "peak" = itag,
      stringsAsFactors = F
      )

  }))
  return(BED12)
}
