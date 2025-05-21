# For the computation of coordinates
coord_pack = function() {
  yr_packages = installed.packages()
  rq_packages = c("devtools", "strawr", "Matrix", "prodlim", "mgcv", "nlme", "prodlim", "parallel", "data.table", "BiocManager")
  
  for (pack in rq_packages) {
    if (!pack %in% yr_packages[,1]){
      install.packages(pack)
    }
  }
  
  if(!"GenomicFeatures" %in% yr_packages[,1]){
  	BiocManager::install("GenomicFeatures")
  }
  
  if(!"FLAMINGOr" %in% yr_packages[,1]){
  	devtools::install("R/FLAMINGO-main/FLAMINGOr")

  }
}