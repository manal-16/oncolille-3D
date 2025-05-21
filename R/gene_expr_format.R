# load packages ----------------------------------------------------------------
library(dplyr)
library(janitor)
library(tidyselect)
library(purrr)

# set working directory -------------------------------------------------------
#setwd()

# upload files -------------------------------------------------------------------
files = list.files(path = "data/hic_bladder_results/geneExpr/",
                   pattern = "genes.results") #or files = c("file.txt", "file2.txt")

# function ---------------------------------------------------------
file_to_df = function(file #expression file to format
                      ){
  
  # get the cell type name from the file name
  long_tumor_name = sub("^GSE148079_(.*)\\.genes\\.results.*$", #according to the file name structure
                        "\\1",
                        file
                        );
  tumor_name = gsub("_?RNASeq?", "",
                    long_tumor_name);
  
  # name standardization
  df_txt = read.table(paste0("data/hic_bladder_results/geneExpr/",
                             file),
                      header = TRUE,
                      sep = "\t"
                      ) %>%
    clean_names();

  if ("gene" %in% colnames(df_txt)){
    df_txt = df_txt %>%
      dplyr::rename(gene_id = gene)
  };
  
  # data frame with gene ensembl identification and TPM value
  df = df_txt %>%
    transmute(gene_id = gene_id,
              !!tumor_name := tpm #to fetch the tumor name
              );
  
  return(df)
  
}

# export ------------------------------------------------------------------
df_files = map(files, file_to_df) %>% #to apply the function to the list of files
  reduce(inner_join,
         by = "gene_id" #to only get genes common to all files
         )  

rownames(df_files) = df_files$gene_id #gene identification indexed
df_merged = df_files %>%
  dplyr::select(-gene_id) %>%
  write.table("data/new/gene_expr_mat.txt",
              sep = "\t"
              )

