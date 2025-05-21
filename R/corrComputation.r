# packages ----------------------------------------------------------------
library(dplyr)

# corrComputation function ------------------------------------------------
corrComputation <- function(genePos, geneExpr, n=10, ...) {
  
  # Initialise the two report vectors.
  geneNames <- vector();
  gene3Dcors <- vector();
  # Get data that is overlaping.
  genePos = genePos
  geneExpr = geneExpr
  overlapGenes <- intersect(rownames(geneExpr),
                            rownames(genePos)
                            ); 

  head(overlapGenes)
  genePos <- genePos[overlapGenes,]; 
  geneExpr <- geneExpr[overlapGenes,]; 
  
  # Calculate the distance and correlation matrices.
  geneDist <- as.matrix(dist(genePos[,c("x", "y", "z"), drop = FALSE],
                             diag=TRUE,
                             upper=TRUE)
                        ); 

  # Calculate the 3D gene expresion score.
  for (gene in rownames(geneDist)) {
    closeGenesNames <- names(sort(geneDist[gene,],
                                  method = "shell"
                                  )[2:(n+1)]
                             );
    gene3Dcorr = 0;
    geneExprTemp <- rbind(geneExpr[closeGenesNames,],
                          geneExpr[gene,]
                          )
    geneCorr <- cor(t(geneExprTemp),
                    method="spearman",
                    use = "complete.obs"
                    );
    
    for(name in closeGenesNames){
      corr = (abs(geneCorr[(gene),][(name)]));
      if(is.na(corr)) {corr = 0}
      gene3Dcorr = gene3Dcorr + corr;
    }
    
    geneNames <- append(geneNames, gene);
    gene3Dcors <- append(gene3Dcors, gene3Dcorr);
  }
  
  transMap3D <- data.frame(corr3D=gene3Dcors);
  rownames(transMap3D) <- geneNames;
  remove(gene3Dcors, geneNames, gene3Dcorr, geneCorr, geneExprTemp, geneExpr, geneDist, genePos, geneExpr, overlapGenes)
  return(transMap3D)
}


# transcriptomeMap1D function ------------------------------------------------------
transcriptomeMap1D <- function(geneStart, geneExpr, n=10, ...) {
  
  # Initilise report data frame.
  transMap1D <- data.frame(corr1D = NULL, chromosome = NULL);
  # Get data that is overlapping.
  overlapGenes <- intersect(rownames(geneExpr), rownames(geneStart)); 
  geneStart <- geneStart[overlapGenes,];
  geneExpr <- geneExpr[overlapGenes,];
  
  # Do the 1D analysis per chromosome.
  for (chr in unique(geneStart$chromosome)) {
    # Initialize the report vector
    gene1Dcors <- vector();
    gsChr <- subset(geneStart, chromosome == chr);
    genesChr <- rownames(gsChr);
    
    # Calculate the distance and correlation matrices.
    geneDist <- as.matrix(dist(gsChr[,"geneStart", drop=FALSE],
                               diag=TRUE, upper=TRUE)); 
    geneCorr <- cor(t(geneExpr[genesChr,]), method="spearman", use = "complete.obs");
    
    # Calculate the 3D gene expression score.
    if (length(genesChr)>n){
    for (gene in genesChr) {
      closeGenesNames <- names(sort(geneDist[gene,])[2:(n+1)]);

      for(close in closeGenesNames){
        if(is.na(geneCorr[gene,][closeGenesNames][close])){
          geneCorr[gene,][closeGenesNames][close] = 0
        }

      }

      gene1Dcorr <- sum(abs(geneCorr[gene,][closeGenesNames]));
      gene1Dcors <- append(gene1Dcors, gene1Dcorr);
    }
    
    transMapChrom <- data.frame(corr1D = gene1Dcors,
                                chromosome = chr
                                );
    rownames(transMapChrom) <- genesChr;
    transMap1D <- rbind(transMap1D,
                        transMapChrom
                        );
  }}
  
  return(transMap1D)
}


# correlatedGenes function ------------------------------------------------
correlatedGenes <- function(chr, matrix= geneStartExp) {
  
  corr1D <- matrix[which(matrix$chromosome.x == chr),c("corr1D","Row.names")]
  mean_value = mean(corr1D$corr1D, na.rm = TRUE)
  corr1D$corr1D[is.na(corr1D$corr1D)] = mean_value 
  MAD <- median(corr1D$corr1D) + 2 * (mad(corr1D$corr1D))
  namesGenes <- as.matrix(subset(corr1D$Row.names,
                                 corr1D$corr1D >= MAD
                                 )
                          )
  MyList <-list(MAD,namesGenes)
  names(x = MyList) = c("MAD","CorrelatedGenes")
  
  return(MyList)
}

# incorrect_format_table function -----------------------------------------

incorrect_table_format=function(table, type=c("coords", "expre", "pos")){
  if(type == "coords"){
    for (col in c("symbol","chromosome","x", "y","z")){if(col %in% colnames(table)){}else{return(TRUE)}};
    return(FALSE)
  }
  if(type == "expre"){
    for (col in c("symbol","chromosome","x", "y","z", "gene_id", "geneStart")){if(col %in% colnames(table)){return(TRUE)}};
    return(FALSE)   
  }
  if(type == "pos"){
    for (col in c("gene_id","chromosome","geneStart")){if(col %in% colnames(table)){}else{return(TRUE)}};
    return(FALSE)
  }
}
