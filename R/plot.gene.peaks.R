
#' Creates a pile up of bed files using the Sushi R package
#'
#' @param GENE A 'character' gene id corresponding to the gene_id's found in the GTF files and the 'name' column in the peak files
#' @param PEAKS data frame containing the following columns, and potentially extras, usually found in a BED12 file, base 0 system
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
#'   \item{sample}{sample_id of samples}
#' }
#' @param GTF The GTF file used to generate the peaks. This is used to determine the genomic coordinates of the gene.
#' @param OUTPUTDIR Output directory
#'
#' @export
#'
#' @examples
plot.gene.peaks = function(
  GENE,
  PEAKS,
  GTF,
  OUTPUTDIR = "."
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

  Sushi::plotBed(
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
  Sushi::plotBed(
    gene.bed,
    chrom = gene.chr,
    chromstart = gene.chromstart,
    chromend = gene.chromend
  )
  Sushi::labelgenome(
    gene.chr,
    gene.chromstart,
    gene.chromend,
    n=2,
    scale="Kb"
  )
  dev.off()

}
