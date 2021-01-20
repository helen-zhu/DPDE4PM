#' Title
#'
#' @param PARAMETERS
#'
#' @return
#'
#' @examples
.read.gtf <- function(PARAMETERS){

  # Creating a TXDB
  op <- options(warn = (-1))
  txdb=makeTxDbFromGFF(PARAMETERS$GTF,format="gtf")
  options(op)

  # Filtering the TXDB
  colkey <- columns(txdb)
  select_col <- match(c("EXONCHROM","TXID","EXONSTART","EXONEND","EXONSTRAND","GENEID","TXNAME"),colkey)
  op <- options(warn = (-1))
  ID = keys(txdb, "TXID")
  temp = select(txdb, ID , c(columns(txdb))[select_col], "TXID")
  select_col2 <- match(c("EXONCHROM","TXID","EXONSTART","EXONEND","EXONSTRAND","GENEID","TXNAME"),names(temp))
  temp <- temp[,select_col2]
  colnames(temp)=c("chr","feature","start","stop","strand","gene","transcript")
  options(op)
  temp$"feature" <- "exon";
  gtf <- temp

  # return data
  return(gtf)
}
