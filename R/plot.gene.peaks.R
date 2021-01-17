
plot.gene.peaks = function(
  GENE,
  PEAKS,
  GTF,
  OUTPUTDIR = "~/figures"
){
  
  # Making a list of parameters to pass back and forth
  PARAMETERS = list()
  PARAMETERS$GENE = GENE
  # PARAMETERS$PEAKS = PEAKS
  PARAMETERS$GTF = GTF
  PARAMETERS$OUTPUTDIR = OUTPUTDIR
  
  # Import GTF as a GRanges Object
  ANNOTATION = .read.gtf(PARAMETERS)
  
  # Gene Bed
  gene.bed = ANNOTATION[ANNOTATION$gene == PARAMETERS$GENE, c("chr", "start", "stop", "gene", "strand")]
  gene.chr = unique(gene.bed$chr)
  gene.chromstart = min(gene.bed$start) - 100
  gene.chromend = max(gene.bed$stop) + 100
  
  # Peak Bed
  peak.bed = PEAKS[PEAKS$name == PARAMETERS$GENE,]
  unique.samples = unique(peak.bed$sample)
  sample_number = structure(1:length(unique.samples), names = unique.samples)
  peak.bed$row = sample_number[peak.bed$sample]
  sample_colour = structure(rainbow(length(unique.samples)), names = unique.samples)
  peak.bed$col = sample_colour[peak.bed$sample]
  
  # Plotting
  filename = paste0(PARAMETERS$OUTPUTDIR, "/", PARAMETERS$GENE, "_Peaks.pdf")
  pdf(filename)
  
  layout(matrix(c(1,2),2, 1, byrow = TRUE), heights = c(0.9,0.1))
  # Bottom, Left, Top, Right
  par(mar=c(0,4,0,4))
  par(oma=c(8,3,3,3))
  
  plotBed(
    peak.bed,
    chrom = gene.chr,
    chromstart = gene.chromstart,
    chromend = gene.chromend,
    # Type
    # type = "circles",
    # Colors
    color=peak.bed$col,
    # Rows
    row = "given",
    rownumber = peak.bed$row
    # rowlabels = unique.samples
  )
  mtext(
    PARAMETERS$GENE,
    side=3, 
    adj=-0.065,
    line=0.5,
    font=2
  )
  plotBed(
    gene.bed,
    chrom = gene.chr,
    chromstart = gene.chromstart,
    chromend = gene.chromend
  )
  labelgenome(
    gene.chr,
    gene.chromstart,
    gene.chromend,
    n=2,
    scale="Kb"
  )
  dev.off()
  
}