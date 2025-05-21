# load packages ----------------------------------------------------------------
library(dplyr)
library(biomaRt)

# set working directory -------------------------------------------------------
#setwd()

# choose reference base --------------------------------------------------------------------
  # load data with biomaRt package

    ##hg38
ensembl_grch38 = useMart("ENSEMBL_MART_ENSEMBL",
                          dataset = "hsapiens_gene_ensembl")

gene_data_38 = getBM(attributes = c("ensembl_gene_id",
                                  "external_gene_name",
                                  "chromosome_name",
                                  "start_position"
                                  ),
                     mart = ensembl_grch38)

    ##hg19
ensembl_grch37 = useMart("ENSEMBL_MART_ENSEMBL",
                          dataset = "hsapiens_gene_ensembl",
                          host = "https://grch37.ensembl.org")

gene_data_37 = getBM(attributes = c("ensembl_gene_id",
                                    "external_gene_name",
                                    "chromosome_name",
                                    "start_position"
                                    ),
                  mart = ensembl_grch37) 

# upload files -------------------------------------------------
  # expression matrix  
expression_file = read.table("data/hic_bladder_results/geneExpr/hic_bladder_2_gene_expr.txt", header = TRUE, sep = "\t")

genes_to_collect = expression_file %>%
  transmute(gene_version = rownames(expression_file),
            gene_id = sub("\\..*" , "", rownames(expression_file)) #to remove the gene version
            )

  # raw coordinates
files = list.files(path = "data/hic_bladder_results/SCABER/", pattern = "res_coordinates_chrom_.*.txt") #to get all chromosomes
files = files[!grepl("format", files)] #to only get the unformatted files

# function ----------------------------------------------------------------
frag_to_gene = function(frag_file, #raw coordinate file
                        resolution, #resolution used for flamingoR extraction
                        gene_database, #genes from database hg38, hg19, etc.
                        genes_to_collect #genes in the expression matrix
                        ) {
  
  print(paste(frag_file, "started."));
  
  # data about genes of interest
  gene_data = inner_join(gene_database,
                         genes_to_collect,
                         join_by(ensembl_gene_id == gene_id)
                         );
  
  # chromosome name
  chrom = sub("res_coordinates_chrom_([0-9XYMT]+)\\.txt", "\\1",
              frag_file);
  
  # raw coordinates
  frag_df = read.table(paste0("data/hic_bladder_results/SCABER/",
                              frag_file),
                       header = TRUE,
                       sep = "\t") %>%
    mutate(frag_id = frag_id * resolution); #to estimate the position on the chromosome
  
  # initialization
  df = data.frame(gene_id = character(),
                  symbol = character(),
                  chromosome = character(),
                  x = numeric(), y = numeric(), z = numeric(),
                  gene_version = character(),
                  stringsAsFactors = FALSE);
  
  many = 0; #to count the number of genes with multiple matches
  none = 0; #to count the number of genes without matches
  
  # iterate through each fragment
  for (row in 1:nrow(frag_df)){
    
    frag = frag_df[row, ];
    
    # filtering the database to match gene start position (bp) to frag_id
    search = gene_data %>%
      filter(start_position > frag$frag_id, 
             start_position < frag$frag_id + resolution, 
             chromosome_name == chrom);
    
    # no match in the database
    if (nrow(search) == 0) {
      none = none +1;
      next;
    }
    
    # a unique match
    else if (nrow(search) == 1){
      df = df %>%
        add_row(gene_id = search$ensembl_gene_id,
                symbol = search$external_gene_name,
                chromosome = search$chromosome_name,
                x = frag$x, y = frag$y, z = frag$z,
                gene_version = search$gene_version
                );
    }
    
    # multiple matches
    else if (nrow(search) > 1){
      many = many + 1;
      #take the first one
      df = df %>% 
        add_row(gene_id = search[1,]$ensembl_gene_id,
                symbol = search[1,]$external_gene_name,
                chromosome = search[1,]$chromosome_name,
                x = frag$x, y = frag$y, z = frag$z,
                gene_version = search$gene_version
                );
    }
  }
  
  # gene identification indexed
  rownames(df) = df$gene_version;
  df = df %>%
    dplyr::select(-gene_id,
           -gene_version
           );
  
  # individual chromosome file
  write.table(df,
              paste0("data/new/", "res_coordinates_chrom_", chrom, "_format", ".txt"),
              sep = "\t"
              );
  
  print(paste(none, "genes had zero match."));
  print(paste(many, "genes had multiple matches."));
  print(paste(frag_file, "finished."));
  
  return(df)
}

# export -------------------------------------------------------------
  # 3D coordinates
df_files = lapply(files, #to apply the function to the list of files
                  frag_to_gene,
                  resolution = 1e4,
                  gene_data = gene_data_37,
                  genes_to_collect = genes_to_collect
                  ) %>%
  bind_rows() %>%
  write.table("data/new/res_coordinates_all_format.txt",
              sep = "\t")

  # 1D coordinates

position = inner_join(genes_to_collect,
                      gene_data_37,
                      join_by(gene_id == ensembl_gene_id)) %>%
  transmute(gene_id = gene_version,
            chromosome = chromosome_name,
            geneStart = start_position
            )

rownames(position) = position$gene_id #gene identification indexed
position = position %>%
  dplyr::select(-gene_id) %>% 
  write.table("data/new/bladder_1d_coords.txt",
              sep = "\t"
              )
