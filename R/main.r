# load packages -----------------------------------------------------------
library(strawr)
library(GenomicFeatures)
library(HiCBricks)

  # FlamingoR functions slightly modified to make normalization optional
source("R/flamingo_coord.r")

# set working directory ---------------------------------------------------
#setwd()

# choose reference base -----------------------------------------------
all_size = getChromInfoFromUCSC("hg19") #hg38, hg19, mm10, etc.

# upload file -------------------------------------------------------------
  # hic file to use FlamingoR on
hic_file = "data/hic_bladder_results/GSE148079_SCABER_HiC.mcool" 

  # fragment resolutions to choose from
readHicBpResolutions(hic_file) #for .hic format only
Brick_list_mcool_resolutions(hic_file) #for .mcool format only

  # available chromosomes to map 
readHicChroms(hic_file) #for hic format only

# FlamingoR ----------------------------------------------------------------
  # FlamingoR predicts 3D coordinates one chromosome at a time (not in parallel)
  # This can be time-consuming depending on the chromosome and resolution

  # adapt the parameters to fit your hic file and desired output
chr_name = "chr2"
res = flamingo_rebuilt(hic_data_low = hic_file, 
                       #the input chromatin interaction data in .hic format, .mcool format or sparse matrix format.
                       file_format = 'mcool', 
                       #format of the input chromatin interaction data. Could be .hic, .mcool or sparse matrix format.
                       domain_res = 1e6, 
                       #size of the domains in bps.
                       frag_res = 1e4, 
                       #size of the fragment in bps (resolution). 
                       chr_name = chr_name, 
                       #name of the chromosome
                       chr_size = all_size[all_size$chrom == chr_name, "size"], 
                       #size of the chromosome in bps.
                       normalization = '', 
                       #normalization method.
                       downsampling_rates = 0.75, 
                       #fraction of contacts to be used for the reconstruction.
                       lambda = 10, 
                       #lagrangian coefficient.
                       max_dist = 0.01, 
                       #the maximum allowed distance between two consecutive points.
                       nThread = 7, 
                       #number of threads available for the reconstruction.
                       n_row = 45000 
                       #the number of rows to use the large matrix format ; reduce the number if the memory is limited.
                       )

# export ------------------------------------------------------------------
write.table(res,
            paste0('data/res_coordinates_chrom_', chr_name, '.txt'),
            sep = "\t",
            row.names = FALSE
            )
