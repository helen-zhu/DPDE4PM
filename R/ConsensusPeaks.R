
ConsensusPeaks = function(
  GENES = "all",
  PEAKS,
  RNA.OR.DNA = "",
  METHOD = c("dpc", "union", "corces"),
  GTF,
  ANNOTATION,
  RESOLUTION = 50,
  DP.ITERATIONS = 1000,
  WEIGHT.THRESHOLD = 0.2,
  N.SD = 1,
  OUTPUTDIR =".",
  PLOT.RESULT=F,
  WRITE.OUTPUT=T,
  OUTPUT.TAG="",
  ALPHA.PRIORS = c(1,2),
  SEED = 123
) {


  if(PARAMETERS$WRITE.OUTPUT){
    filename = paste0(PARAMETERS$OUTPUTDIR, "/", PARAMETERS$GENE, ".", PARAMETERS$OUTPUT.TAG, ".MergedPeaks.tsv")
    write.table(
      OUTPUT.TABLE,
      file = filename,
      sep = "\t",
      col.names = T,
      row.names = F,
      quote = F
    )
  }
  # Return something
}
