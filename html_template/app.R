# packages --------------------------------------------------------------
library(shiny)
library(shinyjs)
library(shinyalert)
library(shinycssloaders)
library(shinyWidgets)
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
library(htmltools)

#library(GWASTools) ; data(centromeres.hg38)

# working directory -------------------------------------------------------
#setwd("C:/Users/localadmin/Documents/STAGE/DATA_SCIENCE/MASTER_1")

# files -------------------------------------------------------------------
#source("essais_R/corrComputation.r") # correlation functions

# -------------------------------------------------------------------------
#ui seulement pour la partie visualisation
##Intra
#Upload
file_1D_coords = fileInput("coord1D",
                           "1D Coordinates (.txt)",
                           accept = ".txt"
                           )
file_expr_intra = fileInput("expre1D",
                            "Expression level (.txt)",
                            accept = ".txt"
                            )
redo_btn_intra = actionBttn("redo_intra",
                            icon = icon("redo") #personnalisation possible
                            )
#Computation
close_nbors_intra = numericInput("neighboursIntra",
                                 "Number of Close Neighbours:",
                                 10,
                                 min = 1,
                                 max = 100
                                 )
question_mrk_intra = circleButton(inputId = "quest_mk_intra",
                                  icon = icon("circle-question"),
                                  #style = "color:darkblue; margin-bottom: 0.9vh;"
                                  )#outline none
compute_intra = actionButton("computeIntra",
                             "Start"
                             )
chrom_barplot = selectInput("choosenchrom",
                            "Chromosome:",
                            c("chr1","chr2","chr3","chr4","chr5", #les ecrire en toute lettre
                              "chr6","chr7","chr8","chr9","chr10","chr11",
                              "chr12","chr13","chr14","chr15","chr16","chr17",
                              "chr18","chr19","chr20","chr21","chr22","chrX","chrY"
                            )
                            )
barplot_intra = withSpinner(plotlyOutput("map1D"))

##InterCorr
#Upload
file_3D_coords_corr = fileInput("coord",
                                "3D Coordinates (.txt)",
                                accept = ".txt"
                                )
file_expr_corr = fileInput("expre",
                           "Expression level (.txt)",
                           accept = ".txt"
                           )
redo_btn_corr = actionBttn("redo_corr",
                            icon = icon("redo") #personnalisation possible
)
close_nbors_corr = numericInput("neighbours",
                                "Number of Close Neighbours:",
                                20,
                                min = 1,
                                max = 1000
                                )
question_mrk_corr = circleButton(inputId = "quest_mk_corr",
                                  icon = icon("circle-question"),
                                  #style = "color:darkblue; margin-bottom: 0.9vh;"
                                 )#outline none
compute_intra = actionButton("compute3D",
                             "Start"
                             )
#plot

##Interclust
#Upload


# ui = shinyUI(page_navbar(shinyjs::useShinyjs(),
#                          file_1D_coords,
#                          file_expr_intra,
#                          redo_btn_intra,
#                          close_nbors_intra,
#                          question_mrk_intra,
#                          compute_intra
# )
# )
# shinyApp(ui = ui, server = function(input, output){})

ui = htmlTemplate(filename="index.html",
                  file_1D_coords = file_1D_coords,
                  file_expr_intra = file_expr_intra
                  )
                       
shinyApp(ui = ui, server = function(input, output, session){})