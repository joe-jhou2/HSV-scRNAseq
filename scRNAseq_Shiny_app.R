# Setup path
setwd("/Users/jhou2/Documents/GitHub/HSV434-HSViz")

# Load libraries
library(shiny)
library(Seurat)
library(dplyr)
library(tidydr)
library(tibble)
library(ggplot2)
library(patchwork)
library(purrr)
library(ComplexHeatmap)
library(circlize)
library(tidyverse)
library(ggh4x)

# Load functions
source("sources/computeUMAP.R")
source("sources/computeBarPlot_CellType_Status.R")
source("sources/computeVln.R")
source("sources/computeHeatmap.R")
source("sources/computeStackBarPlot_Subject.R")
source("sources/computeFeaturePert.R")
source("sources/seurat_analysis_utilities.R")

# Load your dataset
# load("data/HSV_T_Myeloid_combine_Figures.Rdata")
load("data/TClusters.Rdata")
load("data/MyeloidClusters.Rdata")
load("data/Cluster_color.Rdata")

# assign data
# seurat_object = HSV_combine
seurat_object_T = TClusters_integrated
seurat_object_Mye = MyeloidClusters_integrated

# Define UI
ui = fluidPage(
  tags$head(
    tags$style(HTML("
    .shiny-plot-output {
      margin-top: 20px;
      margin-bottom: 0px;
      margin-left: 0px;
    }
  "))
  ),
  titlePanel("Interactive HSV+ Skin Biopsy scRNAseq Data Visualization"),
  
  # Place the dataset selection above the tabs
  div(style = "display: inline-block; width: 100%; padding-bottom: 20px;",
      selectInput("tab_dataset_option", "Select DataSet:",
                  choices = c("", "T cell dataset", "Myeloid cell dataset"),
                  selected = "",
                  multiple = FALSE)
  ),
  
  sidebarLayout(
    sidebarPanel(width = 12,
                 tabsetPanel(
                   tabPanel("Cluster Discovery", value = 6,
                            br(),
                            div(style="display: inline-block;vertical-align:top; width: 19%;",
                                selectInput("celltype_tab_celltype_option", "Select Cluster:",
                                            choices = NULL, 
                                            selected = "All",
                                            multiple = TRUE)),
                            div(style="display: inline-block;vertical-align:top; width: 19%;",
                                selectInput("celltype_tab_subject_option", "Select Subject:",
                                            choices = NULL, 
                                            selected = "All",
                                            multiple = TRUE)),
                            actionButton("run_cluster", "Run Analysis", class = "btn-primary"),
                            tabsetPanel(
                              tabPanel("UMAP", value = 3,
                                       br(),mainPanel(
                                         fluidRow(
                                           column(8, plotOutput("UMAP_all", height = "300px")), # Main UMAP plot
                                         ),
                                         fluidRow(
                                           column(4, plotOutput("UMAP_Prior", height = "300px")), 
                                           column(4, plotOutput("UMAP_Lesion", height = "300px")),
                                           column(4, plotOutput("UMAP_Post", height = "300px"))
                                         )
                                       )
                              ),
                              tabPanel("Cell Type Stat", value = 3,
                                       br(),
                                       mainPanel(
                                         fluidRow(
                                         column(12, uiOutput("dynamicBarplot_CellType_Pert")),
                                         column(12, uiOutput("dynamicBarplot_CellType_Count"))
                                         )
                                       )
                              ),
                              tabPanel("Cell Type by Subject", value = 3,
                                       br(),
                                       mainPanel(
                                         uiOutput("dynamicBarplot")
                                       )
                              ),
                            ),
                   ),
                   tabPanel("Gene Discovery", value = 6,
                            br(),
                            div(style="display: inline-block;vertical-align:top; width: 19%;",
                                selectInput("feature_tab_celltype_option", "Select Cluster:",
                                            choices = NULL, 
                                            selected = "All",
                                            multiple = TRUE)),
                            div(style="display: inline-block;vertical-align:top; width: 19%;",
                                selectInput("feature_tab_subject_option", "Select Subject:",
                                            choices = NULL, 
                                            selected = "All",
                                            multiple = TRUE)),
                            div(style="display: inline-block;vertical-align:top; width: 19%;",
                                textInput("feature_tab_gene_input", "Select Gene (separated by commas):",
                                          value = "PTPRC, CD3D, CD4, CD8A, CD14, CD68, GZMB, IFNG")),
                            actionButton("run_gene", "Run Analysis", class = "btn-primary"),
                            tabsetPanel(
                              tabPanel("Feature Gene Highlight", value = 3,
                                       mainPanel(
                                         fluidRow(
                                           column(6, plotOutput("UMAP_feature", height = "400px")), # Main plot for Feature Highlight
                                           column(6, plotOutput("Dot_split", height = "400px"))
                                         ),
                                         uiOutput("dynamicFeaturePlotsUI") # Placeholders for dynamically generated plots
                                       )
                              ),
                              tabPanel("Heatmap",value = 3,
                                       br(),
                                       mainPanel(uiOutput("dynamicHeatmap"))
                                      ),
                              tabPanel("Violin",value = 3,
                                       br(),
                                       mainPanel(uiOutput("dynamicViolin"))
                                      ),
                              tabPanel("FeaturePert",value = 3,
                                       br(),
                                       mainPanel(uiOutput("dynamicFeaturePert"))
                              )
                            )
                   )
                 )
    ),
    mainPanel(width = 12)
  )
)

# Define server logic
server = function(input, output, session) {
  # Reactive value to store the current dataset
  current_dataset = reactiveVal("")
  # Reactive value to store the notification id
  notify_id = reactiveVal(NULL)
  
  # Reactively determine which Seurat object to use based on the current dataset
  seurat_object_to_use = reactive({
    req(current_dataset())  # Ensure the dataset is selected
    switch(current_dataset(),
           "T cell dataset" = seurat_object_T,
           "Myeloid cell dataset" = seurat_object_Mye,
           stop("No dataset selected")
    )
  })
  
  # Observe changes in the dataset selection and update the current dataset
  observeEvent(input$tab_dataset_option, {
    req(input$tab_dataset_option)  # Ensure there is a dataset selected
    
    # Set the current dataset
    current_dataset(input$tab_dataset_option)
    
    # Get the data based on the selection
    selected_data = getDataset(input$tab_dataset_option)
    
    # Update UI elements with new data
    updateUI(session, selected_data)
  })
  
  # Function to retrieve dataset based on selection
  getDataset = function(dataset_name) {
    if (dataset_name == "T cell dataset") {
      return(list(
        cell_types = levels(seurat_object_T$CellType_Fine),
        subjects = unique(seurat_object_T$Subject),
        data_object = seurat_object_T
      ))
    } else if (dataset_name == "Myeloid cell dataset") {
      return(list(
        cell_types = levels(seurat_object_Mye$CellType_Fine),
        subjects = unique(seurat_object_Mye$Subject),
        data_object = seurat_object_Mye
      ))
    } else {
      stop("Invalid dataset selected")
    }
  }
  
  # Function to update UI elements based on the selected dataset
  updateUI = function(session, data) {
    updateSelectInput(session, "celltype_tab_celltype_option",
                      choices = c("All", data$cell_types))
    updateSelectInput(session, "feature_tab_celltype_option",
                      choices = c("All", data$cell_types))
    updateSelectInput(session, "celltype_tab_subject_option",
                      choices = c("All", data$subjects))
    updateSelectInput(session, "feature_tab_subject_option",
                      choices = c("All", data$subjects))
  }
  
  # Validating Gene Names
  validateGenes = eventReactive(input$run_gene, {
    input_genes = strsplit(input$feature_tab_gene_input, ",\\s*")[[1]]
    input_genes = trimws(input_genes)
    valid_genes = input_genes[input_genes %in% rownames(getDataset(input$tab_dataset_option)$data_object@assays$RNA@data)]
    invalid_genes = setdiff(input_genes, rownames(getDataset(input$tab_dataset_option)$data_object@assays$RNA@data))
    
    if (length(invalid_genes) > 0) {
      showModal(modalDialog(
        title = "Invalid Gene Names Detected",
        paste("The following gene names are not found in the dataset and need correction:", paste(invalid_genes, collapse = ", ")),
        easyClose = TRUE,
        footer = modalButton("Close")
      ))
    }
    
    list(valid_genes = valid_genes, invalid_genes = invalid_genes)
  })
 
  #---------------------------------------------------------------------------#
  # Cluster Discovery Tab
  #---------------------------------------------------------------------------#
  # enable click and run function
  processedClusterData = eventReactive(input$run_cluster, {
    
    message("Run button clicked.")
    
    selected_celltype = if(any(input$celltype_tab_celltype_option == "All")) {
      levels(Idents(seurat_object_to_use()))
    } else {
      input$celltype_tab_celltype_option
    }
    num_selected_celltype = length(selected_celltype)
    
    selected_subject = if(any(input$celltype_tab_subject_option == "All")) {
      unique(seurat_object_to_use()$Subject)
    } else {
      input$celltype_tab_subject_option
    }

    num_selected_subject = length(selected_subject)

    # general theme for UMAP
    UMAPTab_customTheme = function(){
      theme(legend.position = "none",
            axis.title = element_text(size = 0), 
            plot.title = element_text(color="black", size = 16, hjust = 0.5, face = "bold")) 
    }

    list(
      # Plot whole UMAP
      UMAP_all = computeUMAP(seurat_object = seurat_object_to_use(), 
                             selectedCellType = input$celltype_tab_celltype_option, 
                             selectedStatus = "All",
                             selectedSubject = input$celltype_tab_subject_option,
                             CellType_color = Cluster_color),
      # Plot Prior UMAP
      UMAP_Prior = computeUMAP(seurat_object = seurat_object_to_use(), 
                               selectedCellType = input$celltype_tab_celltype_option, 
                               selectedStatus = "Prior",
                               selectedSubject = input$celltype_tab_subject_option,
                               CellType_color = Cluster_color) +
        UMAPTab_customTheme() + ggtitle("Prior"),
      
      # Plot Lesion UMAP
      UMAP_Lesion = computeUMAP(seurat_object = seurat_object_to_use(), 
                                selectedCellType = input$celltype_tab_celltype_option, 
                                selectedStatus = "Lesion",
                                selectedSubject = input$celltype_tab_subject_option,
                                CellType_color = Cluster_color) +
        UMAPTab_customTheme() + ggtitle("Lesion"),
      
      # Plot Post UMAP
      UMAP_Post = computeUMAP(seurat_object = seurat_object_to_use(), 
                              selectedCellType = input$celltype_tab_celltype_option, 
                              selectedStatus = "Post",
                              selectedSubject = input$celltype_tab_subject_option,
                              CellType_color = Cluster_color) +
        UMAPTab_customTheme() + ggtitle("Post"),
      
      # Plot barplot of celltype pert for subject
      Barplot_CellType_Pert = computeBarPlot_CellType_Status(seurat_object = seurat_object_to_use(), 
                                                             selectedSubject = input$celltype_tab_subject_option)$PertBar,
      
      # Plot barplot of celltype count for subject
      Barplot_CellType_Count = computeBarPlot_CellType_Status(seurat_object = seurat_object_to_use(), 
                                                              selectedSubject = input$celltype_tab_subject_option)$CountBar,
      
      # Plot Stacked barplot of celltype pert for subject
      Barplot_Subject_Cluster_Pert = computeStackBarPlot_Subject(seurat_object = seurat_object_to_use(),
                                                                 selectedSubject = input$celltype_tab_subject_option,
                                                                 CellType_color = Cluster_color),
      
      selected_celltype = selected_celltype,
      num_selected_celltype = num_selected_celltype,
      selected_subject = selected_subject,
      num_selected_subject = num_selected_subject
    )
  })
  
  # main UMAP
  output$UMAP_all = renderPlot({
    req(processedClusterData())  # Ensure analysisResults is computed
    plot = processedClusterData()$UMAP_all
    notify_id(showNotification("UMAP plot updated", type = "message", duration = 10))
    plot
    # removeNotification(notify_id())
  })
  
  # Prior Status UMAP
  output$UMAP_Prior = renderPlot({
    req(processedClusterData()) 
    processedClusterData()$UMAP_Prior
  })
  
  # Lesion Status UMAP
  output$UMAP_Lesion = renderPlot({
    req(processedClusterData()) 
    processedClusterData()$UMAP_Lesion
  })
  
  # Post Status UMAP
  output$UMAP_Post = renderPlot({
    req(processedClusterData()) 
    processedClusterData()$UMAP_Post
  })
  
  # Barplot, Celltype pert
  output$dynamicBarplot_CellType_Pert = renderUI({
    req(processedClusterData()) 
    numCellType = length(unique(getDataset(input$tab_dataset_option)$cell_types)) #processedClusterData()$num_selected_celltype
    numfacetRow = round(numCellType/6)
    dynamicHeight = numfacetRow * 150 #max(100, 100 + numfacetRow * 100)
    plotOutput("Barplot_CellType_Pert", height = paste0(dynamicHeight, "px"), width = "800px")
  })
  
  output$Barplot_CellType_Pert = renderPlot({
    req(processedClusterData())  
    plot = processedClusterData()$Barplot_CellType_Pert
    showNotification("Bar plot (Pert) updated", type = "message", duration = 10)
    plot
  })
  
  # Barplot, Celltype count
  output$dynamicBarplot_CellType_Count = renderUI({
    req(processedClusterData())  
    numCellType = length(unique(getDataset(input$tab_dataset_option)$cell_types)) #processedClusterData()$num_selected_celltype
    numfacetRow = round(numCellType/6)
    dynamicHeight = numfacetRow * 150 #max(100, 100 + numfacetRow * 100)
    plotOutput("Barplot_CellType_Count", height = paste0(dynamicHeight, "px"), width = "800px")
  })
  
  output$Barplot_CellType_Count = renderPlot({
    req(processedClusterData())  
    plot = processedClusterData()$Barplot_CellType_Count
    showNotification("Bar plot (Count) updated", type = "message", duration = 10)
    plot
  })
  
  # Barplot, by Subject
  output$dynamicBarplot = renderUI({
    req(processedClusterData())  
    numSubjects = processedClusterData()$num_selected_subject 
    dynamicHeight = max(150, 150 + numSubjects * 35)
    plotOutput("Barplot_Subject_Cluster_Pert", height = paste0(dynamicHeight, "px"))
  })
  
  output$Barplot_Subject_Cluster_Pert = renderPlot({
    req(processedClusterData())
    plot = processedClusterData()$Barplot_Subject_Cluster_Pert
    # showNotification("Bar plot updated", type = "message", duration = 60)
    plot
  })
  
  #----------------------------------------------------------------------------#
  # Gene Discovery Tab
  #----------------------------------------------------------------------------#
  # Use eventReactive to perform data processing when the button is clicked
  processedGeneData = eventReactive(input$run_gene, {
    
    message("Run button clicked.")
    
    # Trigger validateGenes and store the result
    validation_results = validateGenes()
    
    # Halt execution if there are invalid genes
    validate(
      need(length(validation_results$invalid_genes) == 0, "Please correct the invalid gene names before proceeding.")
    )
    
    # Get genes from valid gene list and get cluster from sharedinput
    valid_genes = validation_results$valid_genes
    
    # Calculate the number of gene plots to render
    num_valid_genes = length(valid_genes) #min(length(valid_genes), 8)

    selected_celltype = if(any(input$feature_tab_celltype_option == "All")) {
      levels(Idents(seurat_object_to_use()))
    } else {
      input$feature_tab_celltype_option
    }
    
    num_selected_celltype = length(selected_celltype)

    # Perform your data extraction or analysis function
    list(
      UMAP_feature = computeUMAP(seurat_object = seurat_object_to_use(), 
                                 selectedCellType = input$feature_tab_celltype_option, 
                                 selectedStatus = "All",
                                 selectedSubject = input$feature_tab_subject_option,
                                 CellType_color = Cluster_color),

      Dot_split = DotPlot(subset(seurat_object_to_use(), 
                                 subset = ("All" %in% input$feature_tab_subject_option | Subject %in% input$feature_tab_subject_option)), 
                          idents = selected_celltype, 
                          features = valid_genes, split.by = "Status", 
                          cols = c("green", "lightgrey", "red")) + RotatedAxis() + theme(axis.title = element_blank()),

      heatmap = computeHeatmap(seurat_object = seurat_object_to_use(), 
                               selectedCellType = input$feature_tab_celltype_option,
                               selectedStatus = "All",
                               selectedSubject = input$feature_tab_subject_option,
                               selectedFeature = valid_genes,
                               CellType_color = Cluster_color),

      vln = computeVln(seurat_object = seurat_object_to_use(), 
                       selectedCellType = input$feature_tab_celltype_option,
                       selectedStatus = "All",
                       selectedSubject = input$feature_tab_subject_option,
                       selectedFeature = valid_genes,
                       CellType_color = Cluster_color),
      
      feaPert = computeFeaturePert(seurat_object = seurat_object_to_use(), 
                                   selectedCellType = input$feature_tab_celltype_option,
                                   selectedSubject = input$feature_tab_subject_option,
                                   selectedFeature = valid_genes,
                                   CellType_color = Cluster_color),

      num_valid_genes = num_valid_genes,
      valid_genes = valid_genes,
      selected_celltype = selected_celltype,
      num_selected_celltype = num_selected_celltype
    )
  })
  
  # Tab: Feature
  output$UMAP_feature = renderPlot({
    req(processedGeneData())  # Ensure analysisResults is computed
    plot = processedGeneData()$UMAP_feature
    showNotification("UMAP plot updated", type = "message", duration = 10)
    plot
  })
  
  output$Dot_split = renderPlot({
    req(processedGeneData())
    processedGeneData()$Dot_split
  })
  
  # Dynamically create UI placeholders for feature plots
  output$dynamicFeaturePlotsUI = renderUI({
    req(processedGeneData()) # Ensure data is processed
    
    numGenes = processedGeneData()$num_valid_genes
    
    # Check if numGenes is NA or NULL
    if(is.na(numGenes) || is.null(numGenes)) {
      numGenes = 0
    }
    
    plotUIs = lapply(seq_len(numGenes), function(i) {
      column(3, plotOutput(outputId = paste0("Feature", i), height = "200px"))
    })
    
    do.call(fluidRow, plotUIs)
  })

  # Observe button click or any specific trigger
  observeEvent(input$run_gene, {
    req(processedGeneData())  # Ensure gene data is processed before proceeding
    
    # Access valid genes from processed data
    valid_genes = processedGeneData()$valid_genes

    
    # Loop through the valid genes to create plots
    for (i in seq_along(valid_genes)) {
      local({
        geneName = valid_genes[i]
        plotName = paste0("Feature", i)
        
        # Assign plot rendering to a dynamic output
        output[[plotName]] = renderPlot({
          
          # Isolate input dependencies to ensure plots only update on button click
          cellTypeOption = isolate(input$feature_tab_celltype_option)
          subjectOption = isolate(input$feature_tab_subject_option)
          
          # Subset the Seurat object based on isolated inputs
          Seurat_subset = subset(seurat_object_to_use(), subset = (
            ("All" %in% cellTypeOption | Idents(seurat_object_to_use()) %in% cellTypeOption) &
              ("All" %in% subjectOption | Subject %in% subjectOption)
          ))
          
          # Generate and return the plot
          FeaturePlot(object = Seurat_subset, features = geneName) + 
              theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
                    panel.background = element_rect(fill = NA), 
                    panel.grid = element_blank(),
                    axis.title = element_text(size = 0), 
                    axis.line = element_blank(),
                    axis.text = element_blank(),
                    axis.ticks = element_blank(),
                    legend.position = "right")

        })
      })
    }
    showNotification("Feature plot updated", type = "message", duration = 10)
  })
  
  # Tab: Heatmap
  output$dynamicHeatmap = renderUI({
    req(processedGeneData()) # Ensure data is processed
    numGene = processedGeneData()$num_valid_genes
    dynamicHeight = max(120, 120 + numGene * 30)  # Calculate height based on the number of subjects
    plotOutput("heatmap", height = paste0(dynamicHeight, "px"))
  })
  
  output$heatmap = renderPlot({
    req(processedGeneData())
    plot = processedGeneData()$heatmap
    showNotification("Heatmap plot updated", type = "message", duration = 10)
    plot
  })
  
  # Tab: Violin
  output$dynamicViolin = renderUI({
    req(processedGeneData()) 
    numGene = processedGeneData()$num_valid_genes
    numCellType = processedGeneData()$num_selected_celltype
    dynamicHeight = max(120, 120 + numGene * 70) 
    dynamicWidth = max(200, 200 + numCellType * 70)
    plotOutput("vln", height = paste0(dynamicHeight, "px"), width = paste0(dynamicWidth, "px"))
  })
  
  output$vln = renderPlot({
    req(processedGeneData())
    plot = processedGeneData()$vln
    showNotification("Violin plot updated", type = "message", duration = 10)
    plot
  })
  
  # Tab: FeaturePert
  output$dynamicFeaturePert = renderUI({
    req(processedGeneData()) 
    numGene = processedGeneData()$num_valid_genes
    numCellType = processedGeneData()$num_selected_celltype
    dynamicHeight = max(120, 120 + numGene * 70)  
    dynamicWidth = max(200, 200 + numCellType * 70)
    plotOutput("feaPert", height = paste0(dynamicHeight, "px"), width = paste0(dynamicWidth, "px"))
  })
  
  output$feaPert = renderPlot({
    req(processedGeneData())
    plot = processedGeneData()$feaPert
    showNotification("Feature Pert plot updated", type = "message", duration = 10)
    plot
  })
}

# Run the app
shinyApp(ui = ui, server = server)
