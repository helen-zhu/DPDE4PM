#' Title
#'
#' @param PARAMETERS
#' @param ANNOTATION
#'
#' @return
#'
#' @examples
.get.gene.anno <- function(PARAMETERS,ANNOTATION) {

  # extract batch annotation
  anno=ANNOTATION[ANNOTATION$gene == PARAMETERS$GENE,c(1,3:6)]
  anno_unique=unique(anno)

  # extract information
  strand=as.character(anno_unique[1,4])
  chr=as.character(anno_unique[1,1])
  left=min(anno_unique$start)
  right=max(anno_unique$stop)
  intervals=anno_unique[,2:3]-left+1
  gene=as.character(anno_unique[1,5])
  dna_length=right-left+1

  # prepare DNA2RNA
  DNA2RNA=rep(0,dna_length)
  no_intervals=length(intervals[,1])
  for (i in 1:no_intervals) {DNA2RNA[intervals[i,1]:intervals[i,2]]=1}
  exome_length=sum(DNA2RNA) # this is actally exome length
  DNA2RNA=cumsum(DNA2RNA)*DNA2RNA


  # summarize result
  batch_anno=list(gene=gene,chr=chr,strand=strand,left=left,right=right,
                  DNA2RNA=DNA2RNA,dna_length=dna_length,
                  exome_length=exome_length)
}
