
# header ------------------------------------------------------------------
#16/04/25
#prototypage nouvelle app : reagencement + reactivite supplementaire
#utilisation de bslib pour esthetisme
#17/04/25
#transposition avec techno ancienne app
#18/04/25
#ajout loading, essai pour telecharger les csv depuis l'app

# packages --------------------------------------------------------------
library(shiny)
library(shinyjs)
library(shinythemes)
library(shinyalert)
library(shinycssloaders)
library(bslib)
library(dplyr)
library(stringr)
library(plotly)
#library(plot3D)
library(Hmisc)
library(RColorBrewer)
library(fpc)
library(dbscan)
library(MatrixModels)
library(GWASTools) ; data(centromeres.hg38)

# working directory -------------------------------------------------------


# files -------------------------------------------------------------------
source("R/corrComputation.r") # correlation functions

# Rshiny: UI ----------------------------------------------------------------
ui <- shinyUI(page_fluid(shinyjs::useShinyjs(),
                         theme = bs_theme(version = 5,
                                          bootswatch = "cosmo"  # chosen palette theme (bslib)
                                          ),
                         
                         page_navbar(title = "3DCorr",
                                     navbar_options = navbar_options(theme = "dark"),
                                     
 # 1st page: Description
 
                                     nav_panel(title = "Introduction of the Tool",
                                               h3(strong("Purpose of the tool")),
                                               p("3Dcorr is made to study the relation between 3D 
                                                         space organisation of chromosomes and gene expression to give information 
                                                         about potential epigenetic deregulation regions involved in the progression of tumors.",
                                                 style = "font-size:20px;"
                                                 ),
                                               h3(strong("How to use")),
                                               p("3Dcorr is composed of various tools each grouped in 2 categories, 
                                                         one called Interchromosomal study in the Visualisation tab which studies 
                                                         the relation between chromosomes (exclusively 3D) and the other called Intrachromosomal study 
                                                         which studies the relation between genes inside of a chromosome (both 3D and 1D available).",
                                                 style = "font-size:20px;"
                                                 ),
                                               p("Each categories has various parameters to explore but they both need standardized files to work,
                                                         Interchromosomal study needs a 3D Coordinates txt file in the following format:",
                                                 style = "font-size:20px;"
                                                 ),
                                               tableOutput('TDCoordEx'), # sample table 1
                                               p("The rowname needs to give information about the name of the gene, 
                                                         both the symbol and chromosome column needs to be a string and the 
                                                         x, y, z coordinates needs to be an integer. If you don't have a 3D coordinates 
                                                         file with the coordinates but you have an hi-c file, check the last tab 'How to Compute 3D coordinates from Hi-C'.",
                                                 style = "font-size:20px;"
                                                 ),
                                               p("The Interchromosomal study tab also needs another file called gene Expression in a txt format file	
                                                         which represent the value of the expression depending on the cell line, the format looks like this:",
                                                 style = "font-size:20px;"
                                                 ),
                                               tableOutput('GeneExprEx'), # sample table 2
                                               p("The Intrachromosomal study needs both the 3D coordinates files and the Gene Expression 
                                                         file but also a 1D coordinates files in a txt file that has the following format",
                                                 style = "font-size:20px;"
                                                 ),
                                               tableOutput('ODCoordEx') # sample table 3
                                               ),
 
 # 2nd page: Visualisation (inter & intra)
 
                                        nav_menu(title = "Visualisation",
            # interchromosomal                                       
                                                 nav_panel(title = "Interchromosomal Study",
                                                           card(card_header("Uploads"),
                                                                card_body(height= "100px",
                                                                          layout_columns(fileInput("coord", "Upload Coordinates file as .txt",accept = ".txt"),
                                                                                         fileInput("expre", "Upload Expression Matrix file as .txt",accept = ".txt")
                                                                                         )
                                                                          )
                                                                ),
                                                           card(card_header("Computation Parameters"),
                                                                navset_tab(
                                                                  nav_panel(title = "Correlation",
                                                                            card_body(height = "100px",
                                                                                      layout_columns(numericInput("neighbours", "Number of Close Neighbours:", 20, min = 1, max = 1000),
                                                                                                     actionButton("compute3D", "Compute", icon = icon("redo")),
                                                                                                     col_widths = c(4, 2)
                                                                                                     )
                                                                                      ),
                                                                            card_footer("explanation (a rediger)")
                                                                            ),
                                                                  nav_panel(title = "Clustering",
                                                                            card_body(height = "195px",
                                                                                      layout_columns(numericInput("minpoints", "MinPoints for DBScan:", 15, min = 1, max = 100),
                                                                                                     numericInput("neighbours2", "Number of Close Neighbours:", 20, min = 1, max = 1000),
                                                                                                     actionButton("computeCluster", "Compute", icon = icon("redo")),
                                                                                                     numericInput("distance_group", "Distance of Circle:", 0.02, min = 0.001, max = 1, step = 0.001),
                                                                                                     sliderInput("neighbour_factor", "Neighbour Coefficient:", min = 0, max = 1, step  =0.05, value = 0.8),
                                                                                                     col_widths = c(4, 4, 2, 4, 4)
                                                                                                     )
                                                                                      ),
                                                                            card_footer("explanation (a rediger)")
                                                                            )
                                                                           ),
                                                                ),
                                                           layout_columns(withSpinner(uiOutput("conditional_card_corr")), # outputs if correlation button pressed
                                                                          withSpinner(uiOutput("conditional_card_clust")) # outputs if clustering button
                                                                          )
                                                           ),
            # intrachromosomal
                                                 nav_panel("Intrachromosomal Study",
                                                           card(card_header("Uploads"),
                                                                card_body(height= "100px",
                                                                          layout_columns(fileInput("coord1D", "Upload 1D Coordinates file as .txt",accept = ".txt"),
                                                                                         fileInput("expre1D", "Upload Expression file as .txt", accept = ".txt"),
                                                                                         fileInput("intrachrom3D", "Upload 3D Coordinates file as .txt",accept = ".txt")
                                                                                         )
                                                                          )
                                                                ),
                                                           card(card_header("Computation Parameters"),
                                                                card_body(height = "100px",
                                                                          layout_columns(numericInput("neighboursIntra", "Number of Close Neighbours:", 10, min = 1, max = 100),
                                                                                         actionButton("computeIntra","Compute", icon = icon("redo")),
                                                                                         col_widths = c(4, 2)
                                                                                         )
                                                                          ),
                                                                card_footer("explanation (a rediger)")
                                                                ),
                                                           withSpinner(uiOutput("conditional_card_1D"))
                                                           )
                                                 ),
 
 # 3rd page: FlamingoR
 
                                        nav_panel(title = "How to Compute 3D coordinates from Hi-C",
                                                  h3(strong("FlamingoR")),
                                                  p("3D coordinates of genes in chromomes are computed thanks to",
                                                    tags$a(href="https://github.com/wangjr03/FLAMINGO", "FlamingoR"),
                                                    "which is a tool available in R that uses Hi-C file in .hic and .mcool format to predict
                                                         the 3D space organisation of chromosomes thanks to Hi-C informations.",
                                                    style="font-size:20px;"
                                                    ),
                                                  p("The tool requires several information to be able to work properly, but samples are available below.")
                                                  ),
                                        nav_spacer(),
 # Useful links (optional)
                                        nav_menu(
                                          title = "Links",
                                          align = "right",
                                          nav_item(tags$a("Posit",
                                                          href = "https://posit.co"
                                                          )
                                                   )
                                                 )
                                     )
  
                         )
              )


# Rshiny: server ----------------------------------------------------------
server = function(input, output){
  
  # INITIALIZATION
  
  # parameters
  options(shiny.maxRequestSize = 30*1024^2)
  shinyjs::disable("checkBoxViewLevels")
  hideSpinner("conditional_card_corr") # hidden until "compute3D" button pressed
  hideSpinner("conditional_card_clust")
  hideSpinner("conditional_card_1D")
  hideSpinner("plot_corr_color") # hidden until "draw3D" button pressed
  hideSpinner("plot")
  hideSpinner("map1D")
  hideSpinner("map3D")
  
  # variables
  genePos = NULL
  genePosExp = NULL
  genePosExp_Intra = NULL
  genePosExp_Cluster = NULL
  geneExprSing = NULL
  clusters = NULL
  top_clusters = NULL
  
  #1st page  
  
  # sample table 1
  TDCoordSample <- data.frame("TSPAN6","X",0.0515775082452429,0.114731643806213,0.0363846874093609)
  colnames(TDCoordSample) <- c("symbol", "chromosome", "x", "y", "z")
  rownames(TDCoordSample) <- c("ENSG00000000003.10")
  
  # sample table 2
  GeneExprSample <- data.frame(21.28, 16.62, 57.1, 22.86)
  colnames(GeneExprSample) <- c("Tumor_T3", "HT1376_1", "RT4_1", "SCABER_1")
  rownames(GeneExprSample) <- c("ENSG00000000003.10")
  
  # sample table 3
  ODCoordSample <- data.frame("chrX", 99883667)
  colnames(ODCoordSample) <- c("chromosome", "geneStart")
  rownames(ODCoordSample) <- c("ENSG00000000003.10")
  
  output$TDCoordEx <- renderTable(TDCoordSample, rownames = TRUE)
  output$GeneExprEx <-renderTable(GeneExprSample, rownames = TRUE)
  output$ODCoordEx <- renderTable(ODCoordSample, rownames = TRUE)
  
  # 2nd page 
  
  # Inter chromosomal computation of correlation
  
  dataframe = observeEvent(input$compute3D,{
    if(is.null(input$coord) || is.null(input$expre)){
      return(NULL)
    }
    
    # if(is.na(input$neighbours)||
    #    input$neighbours<1||input$neighbours>1000|| # see min & max numericInput
    #    ceiling(input$neighbours)/floor(input$neighbours)!=1){ # only integer
    #   shinyalert( # a warning for the user, in order for the app to not shut down
    #     title = "Incorrect Neighbour Input",
    #     text = "Choose a number between 1 and 1000",
    #     size = "s", 
    #     closeOnEsc = TRUE,
    #     closeOnClickOutside = TRUE,
    #     html = FALSE,
    #     type = "error",
    #     showConfirmButton = TRUE,
    #     showCancelButton = FALSE,
    #     confirmButtonText = "OK",
    #     confirmButtonCol = "#AEDEF4",
    #     timer = 0,
    #     imageUrl = "",
    #     animation = TRUE
    #   );
    #   return(NULL)
    #   }
      
    showSpinner("conditional_card_corr") # shows spinner
    
    # card with the visualization options
    output$conditional_card_corr = renderUI({card(card_header("Plot Correlation Color"),
                                                  full_screen = TRUE,
                                                  layout_sidebar(sidebar = sidebar(bg = "lightgrey",
                                                                                   selectInput("choosenchrom_inter",
                                                                                               "Chromosome:",
                                                                                               c("All","Chr1","Chr2","Chr3","Chr4","Chr5",
                                                                                                 "Chr6","Chr7","Chr8","Chr9","Chr10","Chr11",
                                                                                                 "Chr12","Chr13","Chr14","Chr15","Chr16","Chr17",
                                                                                                 "Chr18","Chr19","Chr20","Chr21","Chr22","ChrX","ChrY"
                                                                                                 )
                                                                                               ),
                                                                                   sliderInput("cor_value", "Minimum Correlation Value:", min = 0, max = 100, value = 0),
                                                                                   sliderInput("x_coord_cut", "X coordinates Cut:", min = -1.0, max = 1.0, value = c(-1.0, 1.0), step = 0.1),
                                                                                   sliderInput("y_coord_cut", "Y coordinates Cut:", min = -1.0, max = 1.0, value = c(-1.0, 1.0),step = 0.1),
                                                                                   sliderInput("z_coord_cut", "Z coordinates Cut:", min = -1.0, max = 1.0, value = c(-1.0, 1.0),step = 0.1)
                                                                                   ),
                                                                 withSpinner(plotlyOutput("plot_corr_color")
                                                                             )
                                                                 )
                                                  )
                                             }
                                            )
    
    print("Computation Start")

    # creation df from uploaded files
    # geneExpr = read.table(input$expre$datapath[1], header = TRUE) # expression level
    # genePos <<- read.table(input$coord$datapath[1]) # 3D coordinates
    
  #   if (incorrect_table_format(geneExpr, type="expre") || incorrect_table_format(genePos, type = "coords")){
  #   shinyalert( # a warning for the user, in order for the app to not shut down
  #     title = "Incorrect Upload",
  #     text = "Wrong table format",
  #     size = "s",
  #     closeOnEsc = TRUE,
  #     closeOnClickOutside = TRUE,
  #     html = FALSE,
  #     type = "error",
  #     showConfirmButton = TRUE,
  #     showCancelButton = FALSE,
  #     confirmButtonText = "OK",
  #     confirmButtonCol = "#AEDEF4",
  #     timer = 0,
  #     imageUrl = "",
  #     animation = TRUE
  #   );
  #     return(NULL)
  #   }
  # 
    # # calculation correlation between expression level and position
    # genePosExpr = corrComputation(genePos, geneExpr, n = input$neighbours)
    # 
    # # merging by gene id, "Row.names" column created consequently
    # genePosExp <<- merge.data.frame(genePosExpr, genePos, by = "row.names")
    # 
    # print("Computation finished")
    withProgress(message = "Computation in progress", value = 0, {
      incProgress(0.1, detail = "Reading input files")
      geneExpr = read.table(input$expre$datapath[1], header = TRUE)
      genePos <<- read.table(input$coord$datapath[1])
      
      incProgress(0.4, detail = "Calculating correlation")
      genePosExpr = corrComputation(genePos, geneExpr, n = input$neighbours)
      
      incProgress(0.4, detail = "Merging data")
      genePosExp <<- merge.data.frame(genePosExpr, genePos, by = "row.names")
      
      incProgress(0.1, detail = "Finalizing")
    })

  })
  
  hideSpinner("conditional_card_corr") # hides spinner 

  # Inter chromosomal drawing of correlation 
  
  draw = observeEvent(input$compute3D,{ 
    output$plot_corr_color <- renderPlotly({
      if(is.null(genePosExp)){
        return(NULL)
      }
      
      showSpinner("plot_corr_color") # shows spinner
      
      # filter according to the user choices
      genePosExpDraw <- genePosExp %>%
        filter(corr3D >= input$cor_value,
               x >= input$x_coord_cut[1],
               x <= input$x_coord_cut[2],
               y >= input$y_coord_cut[1],
               y <= input$y_coord_cut[2],
               z >= input$z_coord_cut[1],
               z <= input$z_coord_cut[2]
        )
      print(head(genePosExpDraw))
      
      if(input$choosenchrom_inter != "All"){
        chrom = str_remove(input$choosenchrom_inter, "Chr") 
        genePosExpDraw <- genePosExpDraw %>%
          filter(chromosome == chrom) 
      }
      
      # plotly graph
      plot = plot_ly(genePosExpDraw, 
                     x = genePosExpDraw$x, 
                     y = genePosExpDraw$y, 
                     z = genePosExpDraw$z,
                     # color scale 
                     marker = list(color = genePosExpDraw$corr3D, 
                                   symbol = "circle",
                                   size = 3,
                                   colorscale = "Jet",
                                   showscale = TRUE
                     ),
                     # label
                     text = ~paste("Chromosome ",
                                   chromosome,
                                   "<br>Symbol",
                                   symbol,
                                   "<br>Correlation Score",
                                   corr3D
                     )
      )
      gc() # for the report on memory usage
      return(plot)
    }
    )
    
  })
  
  # Cluster computation of Inter chromosomal correlation for each genes
  
  dataframe_cluster = observeEvent(input$computeCluster,{
    if(is.null(input$coord) || is.null(input$expre)){
      return(NULL)
    }
    if(is.na(input$neighbours2)||
       input$neighbours2<1||input$neighbours2>1000|| # see min & max numericInput
       ceiling(input$neighbours2)/floor(input$neighbours2)!=1){ # only integer
      shinyalert( # a warning for the user, in order for the app to not shut down
        title = "Incorrect Neighbour Input",
        text = "Choose a number between 1 and 1000",
        size = "s", 
        closeOnEsc = TRUE,
        closeOnClickOutside = TRUE,
        html = FALSE,
        type = "error",
        showConfirmButton = TRUE,
        showCancelButton = FALSE,
        confirmButtonText = "OK",
        confirmButtonCol = "#AEDEF4",
        timer = 0,
        imageUrl = "",
        animation = TRUE
      );
      return(NULL)
    }
    if(is.na(input$minpoints)||
       input$minpoints<1||input$minpoints>100|| # see min & max numericInput
       ceiling(input$minpoints)/floor(input$minpoints)!=1){ # only integer
      shinyalert( # a warning for the user, in order for the app to not shut down
        title = "Incorrect MinPoint Input",
        text = "Choose a number between 1 and 100",
        size = "s", 
        closeOnEsc = TRUE,
        closeOnClickOutside = TRUE,
        html = FALSE,
        type = "error",
        showConfirmButton = TRUE,
        showCancelButton = FALSE,
        confirmButtonText = "OK",
        confirmButtonCol = "#AEDEF4",
        timer = 0,
        imageUrl = "",
        animation = TRUE
      );
      return(NULL)
    }
    if(is.na(input$distance_group)||
       input$distance_group<0.001||input$distance_group>1 # see min & max numericInput
       ){ 
      shinyalert( # a warning for the user, in order for the app to not shut down
        title = "Incorrect Distance of Circle Input",
        text = "Choose a number between 0.001 and 1",
        size = "s", 
        closeOnEsc = TRUE,
        closeOnClickOutside = TRUE,
        html = FALSE,
        type = "error",
        showConfirmButton = TRUE,
        showCancelButton = FALSE,
        confirmButtonText = "OK",
        confirmButtonCol = "#AEDEF4",
        timer = 0,
        imageUrl = "",
        animation = TRUE
      );
      return(NULL)
    }
    showSpinner("conditional_card_clust") # shows Spinner
    
    # card with the visualization options
    output$conditional_card_clust = renderUI({card(card_header("Plot Cluster Color"),
                                                     full_screen = TRUE,
                                                     layout_sidebar(sidebar = sidebar(bg = "lightgrey",
                                                                                      selectInput("choosenchrom_cluster",
                                                                                                  "Chromosome:",
                                                                                                  c("All","Chr1","Chr2","Chr3","Chr4","Chr5",
                                                                                                    "Chr6","Chr7","Chr8","Chr9","Chr10","Chr11",
                                                                                                    "Chr12","Chr13","Chr14","Chr15","Chr16","Chr17",
                                                                                                    "Chr18","Chr19","Chr20","Chr21","Chr22","ChrX","ChrY"
                                                                                                  )
                                                                                      ),
                                                                                      checkboxInput("checkBoxHideGenes", "Hide non-clustered Genes", value = FALSE),
                                                                                      checkboxInput("checkBoxViewLevels", "View by correlation level", value = FALSE),
                                                                                      selectInput("dataset", "Clustering Results (.csv)", choices = c("clusters", "top_clusters")),
                                                                                      downloadButton("downloadData", "Download")
                                                     ),
                                                     withSpinner(plotlyOutput("plot")
                                                     )
                                                     )
    )})
    
    print("Computation Start")
    
    # creating df from uploaded files
    geneExpr = read.table(input$expre$datapath[1], header = TRUE) 
    genePos = read.table(input$coord$datapath[1])
    if(incorrect_table_format(geneExpr, type = "expre")||incorrect_table_format(genePos, type = "coords"))
    {shinyalert( # a warning for the user, in order for the app to not shut down
      title = "Incorrect Upload",
      text = "Wrong table format",
      size = "s",
      closeOnEsc = TRUE,
      closeOnClickOutside = TRUE,
      html = FALSE,
      type = "error",
      showConfirmButton = TRUE,
      showCancelButton = FALSE,
      confirmButtonText = "OK",
      confirmButtonCol = "#AEDEF4",
      timer = 0,
      imageUrl = "",
      animation = TRUE
    );return(NULL)}
    
    # if(is.null(nrow(intersect(rownames(geneExpr),
    #                           rownames(genePos)
    # )))){
    #   shinyalert( # a warning for the user, in order for the app to not shut down
    #     title = "Incorrect File Content",
    #     text = "0 overlapping gene between the files",
    #     size = "s", 
    #     closeOnEsc = TRUE,
    #     closeOnClickOutside = TRUE,
    #     html = FALSE,
    #     type = "info",
    #     showConfirmButton = TRUE,
    #     showCancelButton = FALSE,
    #     confirmButtonText = "OK",
    #     confirmButtonCol = "#AEDEF4",
    #     timer = 0,
    #     imageUrl = "",
    #     animation = TRUE
    #   );
    #   return(NULL)
    # }
    # calculation correlation between expression level and position
    genePosExpr = corrComputation(genePos, geneExpr, n = input$neighbours2)
    
    # merging by gene id, "Row.names" column created consequently 
    genePosExp <<- merge.data.frame(genePosExpr, genePos, by = "row.names")
    
    # initialization of a cluster column in the merged df
    genePosExp_clu = genePosExp %>%
      mutate(cluster = 0)
    
    # definition of a significance level
    MAD_value = median(genePosExp_clu$corr3D) + 2 * (mad(genePosExp_clu$corr3D))
    genePosExp_rest = genePosExp_clu %>%
      filter(corr3D < (input$neighbour_factor * MAD_value)) # below the threshold
    genePosExp_clu = genePosExp_clu %>%
      filter(corr3D >= (input$neighbour_factor * MAD_value)) # above the threshold
    
    # deleting the "Row.names" column by indexing it 
    rownames(genePosExp_clu) = genePosExp_clu$Row.names 
    rownames(genePosExp_rest) = genePosExp_rest$Row.names
    genePosExp_clu = genePosExp_clu[,-1] 
    genePosExp_rest = genePosExp_rest[,-1]
    
    
    # keeping only overlapping genes
    overlapGenes <- Reduce(intersect, 
                           list(rownames(geneExpr),
                                rownames(genePos),
                                rownames(genePosExp_clu)
                           )
    ) 
    genePos <- genePos[overlapGenes,] 
    geneExpr <- geneExpr[overlapGenes,] 
    
    # density based clustering (dbscan)
    df <- as.matrix(genePosExp_clu[,4:6]) 
    print(head(df)) # "x", y", "z" columns
    db <- fpc::dbscan(df,
                      eps = input$distance_group,
                      MinPts = input$minpoints
    ) 
    genePosExp_clu$cluster = c(db$cluster) # assigning a gene to a cluster
    
    print("finished")
    
    # filtering the clusters
    new_d = genePosExp_clu %>%
      count(cluster, name = "gene_in_cluster") %>%
      filter(gene_in_cluster >= input$minpoints) 
    
    genePosExp_clu = filter(genePosExp_clu, 
                            cluster %in% new_d$cluster)
    
    genePosExp_clu <- genePosExp_clu %>%
      arrange(desc(cluster))
    
    # saving in a csv file the clustering
    cluster_to_file = genePosExp_clu %>%
      filter(cluster != 0) %>% 		  
      mutate(k_neighbours = input$neighbours,
             filtrage_facteur_mad = input$neighbour_factor,
             DBSCAN_minpoints = input$minpoints,
             DBSCAN_distance = input$distance_group) %>%
      relocate(k_neighbours, filtrage_facteur_mad, DBSCAN_minpoints, DBSCAN_distance, .before = corr3D)
    
    clusters <<- cluster_to_file
    
    data_vars = cluster_to_file %>%
      group_by(cluster) %>%
      summarise(mean_corr = mean(corr3D), 
                num_pop = n(), 
                n_chrom = n_distinct(chromosome)
      )
    data_fin = data_vars %>%
      group_by(cluster) %>%
      summarise(score = mean_corr * (num_pop/(max(data_vars$num_pop))) * (n_chrom/(max(data_vars$n_chrom)))
      ) %>% 
      arrange(desc(score))
    
    data_top_five = data_fin[1:5, ] 
    new_my_data = filter(cluster_to_file,
                         cluster %in% data_top_five$cluster
    )
    
    top_clusters <<- data_top_five
    
    
    # coloration of cluster
    genePosExp_clu = rbind(genePosExp_clu, genePosExp_rest) %>%
      mutate(opacity = 1.0,
             size = 50,
             color = "grey"
      )
    all_colors = unlist(mapply(brewer.pal,
                               brewer.pal.info$maxcolors,
                               rownames(brewer.pal.info)
    )
    )
    color_palette = sample(all_colors,
                           length(unique(genePosExp_clu$cluster))
    )
    print(head(genePosExp_clu))
    
    color_chooser = 0
    for(gene in rownames(genePosExp_clu)){
      if(genePosExp_clu[gene,]$cluster != 0){ # cluster 0 = grey
        temp = filter(genePosExp_clu,
                      cluster == genePosExp_clu[gene,]$cluster
        ) # focus on one cluster
        if(temp[1,]$color == "grey"){ # no color yet for the cluster
          color_chooser = color_chooser + 1 
          genePosExp_clu[gene,]$color = color_palette[color_chooser]
        }else{
          genePosExp_clu[gene,]$color = temp[1,]$color # already a color
        }
      }
    }
    
    print("pass")
    genePosExp_Cluster <<- genePosExp_clu 
    
  })
  
  # Download csv
  data <- reactive({
    get(input$dataset) # chosen file
  })
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste0(input$dataset, ".csv")
    },
    content = function(file) {
      write.csv2(data(), file)
    }
  )
  
  # Cluster drawing of Inter chromosomal correlation for each genes
  
  draw_cluster = observeEvent(input$computeCluster,{
    output$plot <- renderPlotly({
      if(is.null(genePosExp_Cluster)){
        return(NULL)
      }
      showSpinner("plot")
      print(input$checkBoxHideGenes)
      
      # filtering according to the user choices
      if(input$checkBoxHideGenes){
        genePosExp_Cluster = genePosExp_Cluster %>%
          filter(cluster != 0)
      }
      if(input$choosenchrom_cluster !="All"){
        chrom = str_remove(input$choosenchrom_cluster, "Chr")
        print(chrom)
        genePosExp_Cluster<- genePosExp_Cluster %>%
          filter(chromosome == chrom)
      }
      
      # plotly graph
      if(input$checkBoxViewLevels){ # correlation + cluster
        plot = plot_ly(genePosExp_Cluster,
                       x = genePosExp_Cluster$x,
                       y = genePosExp_Cluster$y,
                       z = genePosExp_Cluster$z,
                       size = genePosExp_Cluster$size,
                       opacity = genePosExp_Cluster$opacity,
                       marker = list(color = genePosExp_Cluster$corr3D,
                                     symbol = "circle",
                                     size = genePosExp_Cluster$size,
                                     colorscale = "Jet",
                                     showscale = TRUE,
                                     opacity = genePosExp_Cluster$opacity
                       ),
                       # label
                       text = ~paste("Chromosome ",
                                     chromosome, 
                                     "<br>Symbol", 
                                     symbol, 
                                     "<br>Correlation Score", 
                                     corr3D, 
                                     "<br> Cluster ", 
                                     cluster
                       )
        )
      }else{
        plot = plot_ly(genePosExp_Cluster, # cluster
                       x = genePosExp_Cluster$x, 
                       y = genePosExp_Cluster$y, 
                       z = genePosExp_Cluster$z, 
                       size = genePosExp_Cluster$size,
                       opacity = genePosExp_Cluster$opacity,
                       marker = list(color = genePosExp_Cluster$color,
                                     symbol = "circle",
                                     line = list(width = 0)
                       ),
                       text = ~paste("Cluster ",
                                     cluster,
                                     "<br> corr",
                                     corr3D,
                                     "<br> Chromosome ",
                                     chromosome,
                                     "<br> symbol",
                                     symbol
                       )
        )
      }
      
      
      return(plot)
      
    })
  })
  
  # SAFE CHECKS DRAWING OPTIONS
  
  observeEvent(input$checkBoxHideGenes,{ 
    if(input$checkBoxHideGenes){
      shinyjs::enable('checkBoxViewLevels')
    }else{
      shinyjs::reset('checkBoxViewLevels')
      shinyjs::disable('checkBoxViewLevels')
    }
  })
  
  # 3rd page
  
  # Intra chromosomal computation of correlation 
  
  dataframe_intra = observeEvent(input$computeIntra,{
    if(is.null(input$intrachrom3D) || is.null(input$expre1D)){
      return(NULL)
    }
    if(is.na(input$neighboursIntra)||
       input$neighboursIntra<1||input$neighboursIntra>1000|| # see min & max numericInput
       ceiling(input$neighboursIntra)/floor(input$neighboursIntra)!=1){ # only integer
      shinyalert( # a warning for the user, in order for the app to not shut down
        title = "Incorrect Neighbour Input",
        text = "Choose a number between 1 and 1000",
        size = "s", 
        closeOnEsc = TRUE,
        closeOnClickOutside = TRUE,
        html = FALSE,
        type = "error",
        showConfirmButton = TRUE,
        showCancelButton = FALSE,
        confirmButtonText = "OK",
        confirmButtonCol = "#AEDEF4",
        timer = 0,
        imageUrl = "",
        animation = TRUE
      );
      return(NULL)
    }
    showSpinner("conditional_card_1D") # shows spinner
    
    # card with the visualization parameters
    output$conditional_card_1D = renderUI({card(card_header("Plot Correlation Color"),
                                                full_screen = TRUE,
                                                layout_sidebar(sidebar = sidebar(bg = "lightgrey",
                                                                                 selectInput("choosenchrom",
                                                                                             "Chromosome:",
                                                                                             c("chr1","chr2","chr3","chr4","chr5",
                                                                                               "chr6","chr7","chr8","chr9","chr10","chr11",
                                                                                               "chr12","chr13","chr14","chr15","chr16","chr17",
                                                                                               "chr18","chr19","chr20","chr21","chr22","chrX","chrY"
                                                                                               )
                                                                                             ),
                                                                                 sliderInput("cor_value_intra", "Minimum Correlation Value:", min = 0, max = 20, value = 0),
                                                                                 sliderInput("x_coord_cut_intra", "X coordinates Cut:", min = -1.0, max = 1.0, value = c(-1.0, 1.0), step = 0.1),
                                                                                 sliderInput("y_coord_cut_intra", "Y coordinates Cut:", min = -1.0, max = 1.0, value = c(-1.0, 1.0),step = 0.1),
                                                                                 sliderInput("z_coord_cut_intra", "Z coordinates Cut:", min = -1.0, max = 1.0, value = c(-1.0, 1.0),step = 0.1),
                                                                                 ),
                                                               withSpinner(plotlyOutput("map1D")),
                                                               withSpinner(plotlyOutput("map3D"))
                                                               )
                                                )
                                           }
                                          )
    
    print("Computation Start")
    
    # creation df from uploaded files
    geneExpr = read.table(input$expre1D$datapath[1], header = TRUE)
    genePos = read.table(input$intrachrom3D$datapath[1])
    
    # correlation
    genePosExpr = corrComputation(genePos,geneExpr, n=input$neighboursIntra)
    genePosExp_Intra <<- merge.data.frame(genePosExpr, genePos, by = "row.names") 
    
    print("Computation finished")
    
  })
  
  hideSpinner("conditional_card_1D")
  
  # Intra chromosomal drawing of correlation 
  
  
  draw_intra = observeEvent(input$computeIntra,{
    
    showSpinner("map1D")
    showSpinner("map3D")
    
    # 3D graph
    output$map3D <- renderPlotly({
      if(is.null(genePosExp_Intra)){
        return(NULL)
      }
      
      genePosExpDraw <- genePosExp_Intra %>%
        filter(corr3D >= input$cor_value_intra,
               x >= input$x_coord_cut_intra[1],
               x <= input$x_coord_cut_intra[2],
               y >= input$y_coord_cut_intra[1],
               y <= input$y_coord_cut_intra[2],
               z >= input$z_coord_cut_intra[1],
               z <= input$z_coord_cut_intra[2]
        )
      
      plot = plot_ly(genePosExpDraw,
                     x = genePosExpDraw$x,
                     y = genePosExpDraw$y,
                     z = genePosExpDraw$z,
                     marker = list(color = genePosExpDraw$corr3D,symbol = "circle",
                                   size = 5,
                                   colorscale = "Jet",
                                   showscale = TRUE
                     ),
                     # label
                     text = ~paste("Chromosome ",
                                   chromosome,
                                   "<br>Symbol",
                                   symbol,
                                   "<br>Correlation Score",
                                   corr3D
                     )
      )
      gc()
      return(plot)
    }
    )
    
    # 1D graph
    output$map1D <- renderPlotly({
      if(is.null(input$coord1D) || is.null(input$expre1D)){
        return(NULL)
      }
      # GWASTools :: centromeres 
      centro = centromeres.hg38
      centro$chrom = paste0("chr",
                            centro$chrom
                            )
      colnames(centro) = c("chromosome",
                           "centromereStart",
                           "centromereEnd"
                           )
      
      
      # creation df from uploaded files
      geneExpr1D = read.table(input$expre1D$datapath[1], header = TRUE) # expression level
      genePos1D = read.table(input$coord1D$datapath[1], header = TRUE) # 1D coordinate
      
      # indexing the gene id
      row.names(genePos1D) = genePos1D[,1]
      genePos1D = genePos1D[,-1]

      # correlation 1D
      geneStartExpr = transcriptomeMap1D(genePos1D, geneExpr1D, n=input$neighboursIntra)
      geneStartExp = merge.data.frame(geneStartExpr, genePos1D, by = "row.names")
      

      
      MyList <- correlatedGenes(input$choosenchrom, geneStartExp)
      MAD <- MyList$MAD
      print(MAD)
      
      geneStartExpChr <-subset(geneStartExp,
                               geneStartExp$chromosome.x == input$choosenchrom)
      print(head(geneStartExpChr))
      
      # plotly graph		
      plot <- plot_ly(x=geneStartExpChr$geneStart,
                      y=geneStartExpChr$corr1D,
                      width = 900,
                      height = 350,
                      type = 'bar',
                      name = "Correlation Value"
      ) %>%
        add_segments(x = 0,
                     xend = max(geneStartExpChr$geneStart),
                     y = MAD,
                     yend = MAD, 
                     name = "MAD value", 
                     line=list(color="#6c25be")
        ) %>%
        # centromere's location
        add_segments(x = centro[which(centro$chromosome== input$choosenchrom),"centromereStart"],
                     xend= centro[which(centro$chromosome== input$choosenchrom),"centromereStart"],
                     y = 0,
                     yend = max(geneStartExpChr$corr1D),
                     line=list(color="Red"),
                     name = "Centromere Start"
        ) %>%
        add_segments(x = centro[which(centro$chromosome== input$choosenchrom),"centromereEnd"],
                     xend= centro[which(centro$chromosome== input$choosenchrom),"centromereEnd"],
                     y = 0, yend = max(geneStartExpChr$corr1D), 
                     line=list(color="Red"),
                     name = "Centromere End"
        ) %>%
        layout(title = "Intrachromosomal Study 1D Visualisation",
               yaxis= list(title="Correlation Value"),
               xaxis=list(title="Gene Range")
        )
      
      return(plot)
      
    })
    
    
  })
}


# Rshiny: execution -------------------------------------------------------
shinyApp(ui = ui, server = server)

