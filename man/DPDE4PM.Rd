% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DPDE4PM.R
\name{DPDE4PM}
\alias{DPDE4PM}
\title{DPDE4PM merged RNA-seq peaks called on different samples}
\usage{
DPDE4PM(
  GENE,
  PEAKS,
  GTF = NULL,
  ANNOTATION = NULL,
  RESOLUTION = 50,
  DP.ITERATIONS = 1000,
  WEIGHT.THRESHOLD = 0.2,
  N.SD = 1,
  OUTPUTDIR = ".",
  PLOT.RESULT = F,
  WRITE.OUTPUT = T,
  OUTPUT.TAG = "",
  ALPHA.PRIORS = c(1, 2),
  SEED = 123
)
}
\arguments{
\item{GENE}{A 'character' gene id corresponding to the gene_id's found in the GTF files and the 'name' column in the peak files}

\item{PEAKS}{A data frame containing the following columns, and potentially extras, usually found in a BED12 file, base 0 system
\describe{
  \item{chr}{chromosomes, same as in GTF file}
  \item{start}{starting position of the peak, base 0}
  \item{end}{end position of the peak, base 0}
  \item{name}{gene id}
  \item{score}{p-value associated with the peak}
  \item{strand}{strand of the gene}
  \item{blockCount}{number of segments in the peak}
  \item{blockSizes}{size of segments in the peak, BED12 notation}
  \item{blockStarts}{starting positions of segments, BED12 notation}
  \item{sample}{sample_id of samples}
}}

\item{GTF}{The GTF file used to generate the peaks. This is used to determine the genomic coordinates of the gene.}

\item{ANNOTATION}{A object created by the read.gtf function. This allows the user to provide an annotation object to save compute time during parallelization.}

\item{RESOLUTION}{The width (bps) used to sample points from the peaks. This is likely optimized by choosing the window size used to generate the peaks.}

\item{DP.ITERATIONS}{Number of iterations used to fit the Dirichlet Process}

\item{WEIGHT.THRESHOLD}{A proportion (out of 1) used to determine the percentage of data accounted for by the GMM to be called a peak}

\item{N.SD}{Number of standard deviations from the mean for each fitted Gaussian that should considered part of the joint peak}

\item{OUTPUTDIR}{Output directory}

\item{PLOT.RESULT}{TRUE or FALSE, whether plots should be generated}

\item{WRITE.OUTPUT}{TRUE or FALSE, whether an output file should be saved}

\item{OUTPUT.TAG}{A character string indicating a tag to track the generated files}

\item{ALPHA.PRIORS}{A length 2 numeric vector indicating alpha, beta of the Gamma distribution from which the concentration weight parameter alpha is drawn, see the R package dirichletprocess}

\item{SEED}{a seed for reproducibility}
}
\value{
A dataframe with BED12 columns and additional columns for each sample in the PEAKS dataframe and associated p-value with that peak.
}
\description{
DPDE4PM merged RNA-seq peaks called on different samples
}
