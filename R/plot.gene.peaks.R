
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
#' @param OUTPUT.TAG A character string indicating a tag to track the generated files
#'
#' @export plot.gene.peaks
#'
#' @examples
plot.gene.peaks = function(
  GENE,
  PEAKS,
  GTF,
  OUTPUTDIR = ".",
  OUTPUT.TAG = ""
){

  # Making a list of parameters to pass back and forth
  PARAMETERS = list()
  PARAMETERS$GENE = GENE
  # PARAMETERS$PEAKS = PEAKS
  PARAMETERS$GTF = GTF
  PARAMETERS$OUTPUTDIR = OUTPUTDIR
  PARAMETERS$OUTPUT.TAG = OUTPUT.TAG

  # Import GTF as a GRanges Object
  ANNOTATION = DPDE4PM:::.read.gtf(PARAMETERS)

  # Plotting Peaks
  plotting.peaks = .retrieve.peaks.as.granges(PEAKS = PEAKS, GENE = PARAMETERS$GENE, DF = T)

  # Gene Bed
  gene.bed = ANNOTATION[ANNOTATION$gene == PARAMETERS$GENE, c("chr", "start", "stop", "gene", "strand")]
  gene.chr = unique(gene.bed$chr)

  # Plotting
  filename = paste0(PARAMETERS$OUTPUTDIR, "/", PARAMETERS$GENE, ".", PARAMETERS$OUTPUT.TAG, ".Peaks.pdf")
  pdf(filename)

  # Code for ggplot
  p1 = ggplot(plotting.peaks, aes(y = sample, x = start, xend = end)) +
    geom_dumbbell() +
    theme_classic() +
    theme(axis.text.y = element_text(size = 0),
          axis.ticks.y = element_blank()) +
    ggtitle(PARAMETERS$GENE) +
    xlab(gene.chr) + ylab("Sample")

  p1 = p1 + annotate("rect",
                     xmin=gene.bed$start,
                     xmax=gene.bed$stop,
                     ymin=-2,
                     ymax=0,
                     alpha=0.2,
                     color="black",
                     fill=rainbow(nrow(gene.bed)))
  print(p1)

  dev.off()

}
