# load packages ---------------------------------------------------------------
  # to download required packages
source("R/packages.r")
coord_pack()

library(strawr)
library(FLAMINGOr)
library(Matrix)
library(mgcv)
library(nlme)

# main function -----------------------------------------------
flamingo_rebuilt <- function(hic_data_low,
                             file_format,
                             domain_res, 
                             frag_res, 
                             chr_size,
                             chr_name,
                             normalization,
                             downsampling_rates,
                             lambda, 
                             max_dist,
                             nThread,
                             alpha = -0.25,
                             max_iter = 500,
                             hic_data_high = NULL,
                             norm_low = NULL,
                             norm_high = NULL,
                             n_row = 45000)
{
  # will delete directory Domain_data 
    unlink("Domain_data", recursive = TRUE) 
    unlink("Genomic_loc", recursive = TRUE)
    unlink("tmp_IF.txt.gz")
    unlink("tmp_PD.txt.gz")
    unlink("intra_domain.Rdata")

  
    # format hic
    if(file_format == 'hic'){
  	  flamingo_low_res_data_obj <- rebuilt.flamingo.construct_flamingo_from_hic(hic_file = hic_data_low,
  	                                                                            resolution = domain_res,
  	                                                                            chr_size = chr_size,
  	                                                                            chr_name = chr_name,
  	                                                                            normalization = normalization,
  	                                                                            alpha = alpha,
  	                                                                            n_row = n_row
  	                                                                            )
  	  
      flamingo_high_res_data_obj <- rebuilt.flamingo.construct_flamingo_from_hic(hic_file = hic_data_low,
                                                                                 resolution = frag_res,
                                                                                 chr_size = chr_size,
                                                                                 chr_name = chr_name,
                                                                                 normalization = normalization,
                                                                                 alpha = alpha,
                                                                                 n_row = n_row
                                                                                 )
      
    }

    
    # format mcool
    if(file_format == 'mcool'){
      flamingo_low_res_data_obj <- rebuilt.flamingo.construct_flamingo_from_mcool(mcool_file = hic_data_low,
                                                                        resolution = domain_res,
                                                                        chr_name=chr_name,
                                                                        normalization = normalization,
                                                                        alpha=alpha,
                                                                        n_row=n_row
                                                                        )
    flamingo_high_res_data_obj <- rebuilt.flamingo.construct_flamingo_from_mcool(mcool_file = hic_data_low,
                                                                        resolution = frag_res,
                                                                        chr_name=chr_name,
                                                                        normalization = normalization,
                                                                        alpha=alpha,
                                                                        n_row=n_row
                                                                        )
    }

    
  	print("Dividing domains...")
  	rebuilt.flamingo.divide_domain(flamingo_obj = flamingo_high_res_data_obj, 
    	domain_res = domain_res, frag_res = frag_res)
  	print("caching datasets...")
  	write_huge_mat(flamingo_high_res_data_obj@IF, "tmp_IF.txt.gz", 
    	nThread = nThread)
  	write_huge_mat(flamingo_high_res_data_obj@PD, "tmp_PD.txt.gz", 
    	nThread = nThread)
  	high_res_n_frag = flamingo_high_res_data_obj@n_frag
  	high_res_chr_name = flamingo_high_res_data_obj@chr_name
  	rm(flamingo_high_res_data_obj)
  	print("Reconstructing backbones...")
  	flamingo_backbone_prediction = flamingo.reconstruct_backbone_structure(flamingo_data_obj = flamingo_low_res_data_obj,
  	                                                                       sw = downsampling_rates,
  	                                                                       lambda = lambda,
  	                                                                       max_dist = max_dist,
  	                                                                       nThread = 1
  	                                                                       )

  	print("Reconstructing intra-domain structures...")
  	flamigo_intra_domain_prediction = rebuilt.flamingo.reconstruct_structure(sw = downsampling_rates, 
    	lambda = lambda, max_dist = max_dist, nThread = nThread)
  	save(flamigo_intra_domain_prediction, file = "intra_domain.Rdata")


 	print("Assembling structures...")
  	flamingo_high_res_IF = data.table::fread("tmp_IF.txt.gz", 
    	sep = "\t", header = F, nThread = nThread)
  	flamingo_high_res_IF = as.matrix(flamingo_high_res_IF)
  	flamingo_high_res_PD = data.table::fread("tmp_PD.txt.gz", 
    	sep = "\t", header = F, nThread = nThread)
 

  	flamingo_high_res_PD = as.matrix(flamingo_high_res_PD)
  	flamingo_high_res_data_obj = new("flamingo", IF = flamingo_high_res_IF, 
    	PD = flamingo_high_res_PD, n_frag = high_res_n_frag, 
    	chr_name = high_res_chr_name)


  	res = flamingo.assemble_structure(flamingo_backbone_prediction_obj = flamingo_backbone_prediction, 
    	flamingo_final_res_data_obj = flamingo_high_res_data_obj, 
    	list_of_flamingo_domain_prediction_obj = flamigo_intra_domain_prediction, 
    	max_iter = max_iter)

  	return(res)
}

# secondary functions -----------------------------------------------------
  # rebuilt.flamingo.construct_flamingo_from_hic()
rebuilt.flamingo.construct_flamingo_from_hic <- function(hic_file,
                                                         normalization,
                                                         resolution,
                                                         chr_name,
                                                         chr_size,
                                                         alpha = -0.25,
                                                         n_row = 45000){
  options(scipen = 999)

  normalized_data = strawr::straw(normalization,
                                  hic_file,
                                  chr_name, 
                                  chr_name, 
                                  unit='BP',
                                  binsize = resolution
                                  )
  
  normalized_data[normalized_data$counts == "NaN",] = 0
  
  n <- ceiling(chr_size/resolution)
  i_ind <- 1 + (normalized_data[,1]/resolution)
  j_ind <- 1 + (normalized_data[,2]/resolution)
  input_if = Matrix::sparseMatrix(i = i_ind,
                                  j = j_ind,
                                  x = normalized_data[,3],
                                  dims = c(n,n)
                                  )
  if(n < n_row){
    input_if <- as.matrix(input_if)
  }else{
    input_if <- convert_huge_mat(input_if)
  }
  
  input_if <- input_if + t(input_if)
  diag(input_if) <- diag(input_if)/2
  pd <- input_if^(alpha)
  res = new('flamingo',
            IF = input_if,
            PD = pd,
            n_frag = n,
            chr_name = chr_name
            )
  
  return(res)
}

  # rebuilt.flamingo.construct_flamingo_from_mcool()
rebuilt.flamingo.construct_flamingo_from_mcool <- function(mcool_file,
                                                           normalization,
                                                           resolution,
                                                           chr_name,
                                                           alpha=-0.25,
                                                           n_row=45000){
  options(scipen = 999)
  all_dir = rhdf5::h5ls(mcool_file)
  parent_dir = all_dir[1,2]
  target_dir = paste(c("", parent_dir, resolution),
                     collapse='/')
  mcool_dat = rhdf5::h5read(mcool_file,target_dir)
  available_normalization = setdiff(names(mcool_dat$bins),
                                    c('chrom','start','end','weight')
  )
  available = TRUE #
  if(!normalization %in% available_normalization){
    print('Normalization method not exist in the .mcool data!')
    available = FALSE
  }
  csr_rawcount = data.frame(bin_1 = mcool_dat$pixels$bin1_id,
                            bin_2 = mcool_dat$pixels$bin2_id,
                            value = mcool_dat$pixels$count)
  print("csr_rawcount")
  print(head(csr_rawcount))
  normalization_file = mcool_dat$bins[[normalization]]
  chr_id <- which(mcool_dat$bins$chrom == chr_name)
  offset <- mcool_dat$indexes$chrom_offset[which(mcool_dat$chroms$name==chr_name)]
  n <- length(chr_id)
  csr_rawcount <- subset(csr_rawcount,csr_rawcount[,1] %in% chr_id & csr_rawcount[,2] %in% chr_id)
  csr_rawcount[,1] <- csr_rawcount[,1]-offset+1
  csr_rawcount[,2] <- csr_rawcount[,2]-offset+1
  csr_rawcount <- as.matrix(csr_rawcount)
  normalization_file <- normalization_file[chr_id]
  if(available){
    for(i in 1:(dim(csr_rawcount)[1])){
      csr_rawcount[i,3] <- csr_rawcount[i,3]/(normalization_file[csr_rawcount[i,1]]*normalization_file[csr_rawcount[i,2]])
    }
    
    
  }
  
  input_if <- Matrix::sparseMatrix(i=csr_rawcount[,1],j=csr_rawcount[,2],x=csr_rawcount[,3])
  if(n<n_row){
    input_if <- as.matrix(input_if)
  }else{
    input_if <- convert_huge_mat(input_if)
  }
  
  ## matrice carre
  if(nrow(input_if)>ncol(input_if)){
    input_if = input_if[1:ncol(input_if),]
  }else if(nrow(input_if)<ncol(input_if)){
    input_if = input_if[,1:nrow(input_if)]
  }
  
  input_if <- input_if + t(input_if)
  diag(input_if) <- diag(input_if)/2
  pd <- input_if^(alpha)
  print("ici")
  print(dim(pd))
  res = new('flamingo',IF=input_if,PD=pd,n_frag=n,chr_name=chr_name)
  
  return(res)
  
}

  # write_huge_mat()
write_huge_mat = function (mat, file_path, nThread = 10){
  n = dim(mat)[1]
  block_size = 1000
  n_block = ceiling(n/block_size)
  pb <- progress::progress_bar$new(total = n_block)
  for (row in 1:n_block) {
    start_id = 1 + (row - 1) * 1000
    end_id = min(row * 1000, n)
    data.table::fwrite(mat[start_id:end_id, ], file_path, 
                       col.names = F, row.names = F, sep = "\t", quote = F, 
                       append = T, nThread = nThread)
    pb$tick()
  }
}

  # convert_huge_mat()
convert_huge_mat <- function(sparse_mat){
  print('Contact map is too large, large matrix mod is on')
  n <- dim(sparse_mat)[1]
  res <- matrix(0,n,n)
  n_bin <- ceiling(n/1000)
  bin_size = 1000
  pb <- progress::progress_bar$new(total = n_bin*n_bin)
  for(i in 1:n_bin){
    row_idx = (1+(i-1)*bin_size):min(n,i*bin_size)
    for(j in 1:n_bin){
      col_idx = (1+(j-1)*bin_size):min(n,j*bin_size)
      res[row_idx,col_idx] <- as.matrix(sparse_mat[row_idx,col_idx])
      pb$tick()
    }
  }
  return(res)
}

  # rebuilt.flamingo.divide_domain()
rebuilt.flamingo.divide_domain <- function(flamingo_obj,
                                           domain_res,
                                           frag_res){
  options(scipen = 999)
  dir.create("Domain_data")
  dir.create('Genomic_loc')
  n = flamingo_obj@n_frag
  bin_size = domain_res/frag_res
  if(bin_size != round(bin_size)){
    stop('The domain resolution is not an integer multiple of the fragment resolution!')
  }
  res = list()
  chr = flamingo_obj@chr_name
  pd = flamingo_obj@PD
  input_if = flamingo_obj@IF
  print('Processing Fragments...')
  n_domain = ceiling(n/bin_size) - 1 
  pb <- progress::progress_bar$new(total = n_domain)
  for(i in 1:n_domain){
    start_id <- (i-1)*bin_size+1
    end_id <- min(dim(pd)[1],i*bin_size)
    print(str(pd))#
    tmp_pd <- pd[start_id:end_id,start_id:end_id]
    tmp_input_if <- input_if[start_id:end_id,start_id:end_id]
    start_loc = ((start_id:end_id)-1)*frag_res+1
    end_loc = (start_id:end_id)*frag_res
    tmp_frag = data.frame(chr=chr,start=start_loc,end=end_loc)
    write.table(as.matrix(tmp_pd),
                paste0("./Domain_data/PD_domain_",i,'.txt'),
                col.names = F,
                row.names = F,
                sep="\t",
                quote=F
    )
    write.table(as.matrix(tmp_input_if),
                paste0("./Domain_data/IF_domain_",i,'.txt'),
                col.names = F,
                row.names = F,
                sep="\t",
                quote=F
    )
    write.table(tmp_frag,
                paste0("./Genomic_loc/Genomic_loc_domain_",i,'.txt'),
                col.names = F,
                row.names = F,
                sep="\t",
                quote=F
    )
    pb$tick()
  }
}

  # rebuilt.flamingo.reconstruct_structure()
rebuilt.flamingo.reconstruct_structure<- function(sw,
                                                  lambda,
                                                  max_dist,
                                                  nThread = 28){
  cl <- parallel::makeCluster(nThread)
  parallel::clusterCall(cl,
                        function() library(FLAMINGOr))
  parallel::clusterExport(cl, c('rebuilt.flamingo.reconstruct_structure_wraper',
                                'check_data_availability'
                                )
                          )
  file <- dir('./Domain_data/')
  n <- length(grep('PD',file))
  res = list()
  domain_id <- c()
  parallel::clusterExport(cl,
                          c("sw","lambda","max_dist"),
                          envir=environment()
                          )
  worker_res = parallel::parSapply(cl,
                                   1:n,
                                   function(x){
    rebuilt.flamingo.reconstruct_structure_wraper(x,
                                                  sw,
                                                  lambda,
                                                  max_dist,
                                                  nThread = 1
                                                  )
  })
  
  for(i in 1:length(worker_res)){
    if(!is.null(worker_res[[i]])){
      res[[i]] <- worker_res[[i]]
      domain_id <- c(domain_id,i)
    }
  }
  res <- res[lengths(res)>0]
  names(res) <- domain_id
  parallel::stopCluster(cl)
  return(res)
}

  # rebuilt.flamingo.reconstruct_structure_wraper()
rebuilt.flamingo.reconstruct_structure_wraper<- function(index,
                                                         sw,
                                                         lambda,
                                                         max_dist,
                                                         nThread = 28
                                                         ){
  
  pd <- read.table(paste0("./Domain_data/PD_domain_",index,'.txt'))
  input_if <- read.table(paste0("./Domain_data/IF_domain_",index,'.txt'))
  file <- dir('./Domain_data/')
  
  if(check_data_availability(pd)){
    res =  flamingo.reconstruct_structure_worker(input_if,
                                                 pd,
                                                 sw,
                                                 lambda,
                                                 max_dist,
                                                 nThread = 1
                                                 )
  }else{
    res = NULL
  }
  return(res)
}

  # check_data_availability()
check_data_availability <- function(x){
  x <- as.matrix(x)
  invalid_id <- which(apply(x,1,min)==Inf)
  if(length(invalid_id) == dim(x)[1]){
    return(F)
  }else{
    return(T)
  }
}
