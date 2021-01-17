source("/cluster/home/helenzhu/code/snakemake_M6A_16_PeakMatTest/Rscripts/000_HEADER.R")

library(rtracklayer)
library(GenomicRanges)

# Importing Peaks ---------------------------------------------------------

all_samples = all_samples[!all_samples %in% exclude_samples]

# MeTPeak
metpeak = do.call(rbind, lapply(all_samples, function(i){
  cat(i, "\n")
  tmp = read.table(paste0("F_MeTPeak_Peaks/", i, "/peak.bed6"), header = F, sep = "\t", stringsAsFactors = F)
  tmp$V7 = i
  tmp
}))

# GRanges Object
metpeak_gr = GRanges(
  seqnames = metpeak$V1,
  IRanges(start = metpeak$V2, end = metpeak$V3),
  strand = metpeak$V6,
  mcols = metpeak[,c("V4", "V5", "V7")]
)

metpeak_gr_list = split(metpeak_gr, metpeak_gr$mcols.V4)

# Calculating Coverage ----------------------------------------------------

tmp = metpeak_gr_list[[1]]
reduced_tmp = reduce(tmp)
tmp2 = tmp[queryHits(findOverlaps(tmp, reduced_tmp[1]))]

# Getting a large vector of coordinates
tmp2 = tile(tmp2, width = 50)
df = data.frame(tmp2)["start"]

# Getting a data frame with the counts
# peak_coverage <- coverage(tmp2)
# bins = tile(reduced_tmp[1], width = 1)[[1]]
# gr1_bins <- GenomicRanges::binnedAverage(bins, peak_coverage, "Coverage")
# df = data.frame(gr1_bins)[,c( "Coverage", "start")]


# Dirichlet Process -------------------------------------------------------

df = scale(as.vector(df$start), center = T, scale = T)
dp = DirichletProcessGaussian(df)
dp = Fit(dp, 1000)

# faithfulTrans = (faithful$waiting - mean(faithful$waiting))/sd(faithful$waiting)
# dp = DirichletProcessGaussian(faithfulTrans)
# dp = Fit(dp, 1000)

filename = "~/figures/testing.DP.pdf"
pdf(filename)
plot(dp)
hist(df)
dev.off()