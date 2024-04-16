#' Create Violin Plots for Gene Expression in Seurat Object
#'
#' This function generates violin plots of selected gene expressions across specified 
#' clusters, statuses, and subjects within a Seurat object. Each plot facet represents a gene by 
#' cell type, with expression values on the y-axis and statuses as categories on the x-axis.
#'
#' @param seurat_object A Seurat object containing single-cell RNA-seq data.
#' @param selectedCellType Optional; a vector of clusters to include in the analysis. 
#'        If NULL, includes all clusters.
#' @param selectedFeature Optional; a vector of gene features to plot.
#'        If NULL, no features are plotted.
#' @param selectedStatus Optional; a vector of statuses to include in the analysis. 
#'        If NULL, includes all statuses.
#' @param selectedSubject Optional; a vector of subjects to include in the analysis.
#'        If NULL, includes all subjects.
#' @param CellType_color
#'  
#' @return A ggplot object representing the violin plot of log-normalized gene expression.
#'         The plot facets by gene and cell type, with expression values by status.
#' @import ggplot2
#' @importFrom reshape2 melt
#' 
computeVln = function(seurat_object, selectedCellType, selectedStatus, selectedSubject, selectedFeature, CellType_color = NULL){
  
  # Extract features
  df_Extraction = feature_extraction(seurat_object, selectedCellType = "All", selectedFeature)
  
  # Filtering based on parameters
  df_Extraction_sel = tryCatch({
    df_Extraction %>%
      filter((any(selectedCellType == "All") | CellType %in% selectedCellType) &
               (any(selectedStatus == "All") | Status %in% selectedStatus) &
               (any(selectedSubject == "All") | Subject %in% selectedSubject))
  },
  error = function(e) {
    # Handle the error by returning NULL or creating a specific output
    return(NULL)  # Returning NULL to indicate an error occurred during subsetting
  }
  )
  
  if (is.null(df_Extraction_sel)) {
    # Return a message or a blank plot
    return(ggplot() + 
             geom_blank() +
             annotate("text", x = 0.25, y = 0.25, 
                      label = "Data\nUnavailable", 
                      size = 6, hjust = 0.5, vjust = 0.5) +
             theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
                   panel.background = element_rect(fill = NA), 
                   panel.grid = element_blank(),
                   axis.text = element_blank(),
                   axis.ticks = element_blank(),
                   axis.title = element_blank(), 
                   legend.title = element_blank())
    )
  }
  
  # reshape data for Vln
  exp_meta_df_reshape = reshape2::melt(df_Extraction_sel %>% select(!c("orig.ident", "Subject")), 
                                       id.vars=c("CellType", "Status"), 
                                       variable.name = "Gene", value.name = "Exp")
  exp_meta_df_reshape$Exp = as.numeric(exp_meta_df_reshape$Exp)

  # change gene name order (optional), prior data already done, but confirm
  exp_meta_df_reshape$CellType = factor(exp_meta_df_reshape$CellType, levels = levels(seurat_object$CellType_Fine))
  exp_meta_df_reshape$Status = factor(exp_meta_df_reshape$Status, levels = c("Prior", "Lesion", "Post"))
  exp_meta_df_reshape$Gene = factor(exp_meta_df_reshape$Gene, levels = levels(exp_meta_df_reshape$Gene))
  
  exp_meta_df_reshape = exp_meta_df_reshape %>% droplevels()
  
  # plotting: each facet is gene by cluster, x is status, y is exp value, 
  g = ggplot(data = exp_meta_df_reshape,
             aes(x = Status, y = Exp, fill = CellType)) +
    geom_violin(scale = 'width',
                draw_quantiles = c(0.25, 0.5, 0.75),
                color = 'black',
                size = 0.3,
                alpha = 0.8) +
    stat_summary(fun = median, geom = "pointrange", color = "black") +
    stat_summary(fun = median, geom = "line", color = "black", aes(group = 1)) +
    scale_fill_manual(values = CellType_color) +
    facet_grid(Gene ~ CellType, scales = "free_x") +
    theme_bw() +
    theme(
      panel.spacing = unit(0.05, "lines"),
      panel.grid = element_blank(),
      strip.background = element_blank(),
      strip.text.x = element_text(size = 14, angle = 0, face = "plain", hjust = 0.5),
      strip.text.y = element_text(size = 14, angle = 0, face = "plain", hjust = 0.5),
      axis.title.x = element_text(size = 0, color = 'black'),
      axis.title.y = element_text(size = 14, color = 'black'),
      axis.text.x = element_text(size = 14, color = 'black', angle = 90, hjust = 1, vjust = 0.5),
      axis.text.y = element_text(size = 14, color = 'black'),
      legend.position = "none"
    ) +
    labs(y = 'Log Normalized Expression')
  
  return(g)
}





