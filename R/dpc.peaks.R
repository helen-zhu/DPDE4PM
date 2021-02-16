
.dpc.peaks = function(GENEPEAKSGR, PARAMETERS){
  
  # Sampling from the set of peaks
  TILED.PEAKS.GR = GenomicRanges::tile(GENEPEAKSGR, width = PARAMETERS$RESOLUTION)
  startvec = data.frame(TILED.PEAKS.GR, stringsAsFactors = F)
  
  # Normalizing Data
  startvec.mean = mean(as.vector(startvec$start))
  startvec.sd = sd(as.vector(startvec$start))
  startvec.sd = ifelse( is.na(startvec.sd) | startvec.sd == 0, 1, startvec.sd) # This is for when there's 1 short peak or all samples have the exact same peak
  startvec.scaled = (as.vector(startvec$start) - startvec.mean)/startvec.sd
  
  # Running Dirichlet Process
  set.seed(PARAMETERS$SEED)
  dp = dirichletprocess::DirichletProcessGaussian(
    y = startvec.scaled,
    g0Priors = c(0, 1, 1, 1),
    alphaPriors = PARAMETERS$ALPHA.PRIORS)
  
  dp = dirichletprocess::Fit(
    dpObj = dp,
    its = PARAMETERS$DP.ITERATIONS,
    updatePrior = F,
    progressBar = TRUE)
  
  # Generating a data.frame from the dp object
  dp_data = data.frame(
    "Weights" = dp$weights,
    "Mu" = c(dp$clusterParameters[[1]]),
    "Sigma" = c(dp$clusterParameters[[2]]),
    stringsAsFactors = F)
  dp_data$Mu = (dp_data$Mu*startvec.sd)+startvec.mean
  dp_data$Sigma = (dp_data$Sigma*startvec.sd)
  
  return(list("dp" = dp, "dp_data" = dp_data, "startvec.scaled" = startvec.scaled, "startvec.mean" = startvec.mean, "startvec.sd" = startvec.sd))
  
}