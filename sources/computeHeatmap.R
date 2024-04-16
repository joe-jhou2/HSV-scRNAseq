
#' Compute and Plot a Heatmap of Gene Expression from a Seurat Object
#'
#' This function generates a heatmap visualizing the expression levels of selected features
#' across specified clusters and statuses in a Seurat object.
#'
#' @param seurat_object A Seurat object containing single-cell RNA-seq data.
#' @param selectedCellType Optional; a character vector of cluster identifiers to include
#'        in the analysis. If NULL, all clusters are included.
#' @param selectedStatus Optional; a character vector of statuses to include in the analysis.
#'        If NULL, all statuses are included.
#' @param selectedFeature Optional; a character vector of gene features to plot.
#'        If NULL, no features are plotted.
#' @param CellType_color 
#'
#' @return A ComplexHeatmap object representing the heatmap, which is also drawn to the current
#'         graphics device. The heatmap shows expression levels of selected features, grouped
#'         by cell type and status.
#'
#' @importFrom ComplexHeatmap Heatmap
#' @importFrom ComplexHeatmap HeatmapAnnotation
#' @importFrom circlize colorRamp2
#' @import reshape2
#' 

computeHeatmap = function(seurat_object, selectedCellType, selectedStatus, selectedSubject, selectedFeature, CellType_color = NULL) {
  
  # Ensure 'selectedCellType', 'selectedStatus', and 'selectedFeature' are NULL or character vectors
  # if(!is.null(selectedCellType) && !is.character(selectedCellType)) stop("'selectedCellType' must be NULL or a character vector.")
  # if(!is.null(selectedStatus) && !is.character(selectedStatus)) stop("'selectedStatus' must be NULL or a character vector.")
  # if(!is.null(selectedFeature) && !is.character(selectedFeature)) stop("'selectedFeature' must be NULL or a character vector.")
  # 
  # Use feature_extraction function to get data; assume this function is defined elsewhere
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
  
  # Prepare matrix for heatmap
  ht_data = as.matrix(t(df_Extraction_sel[, !(names(df_Extraction_sel) %in% c("orig.ident", "Subject", "CellType", "Status"))]))
  
  # Prepare annotation for heatmap
  ht_meta = df_Extraction_sel[, c("CellType", "Status"), drop = FALSE]
  ht_meta$CellType = factor(ht_meta$CellType, levels = levels(seurat_object$CellType_Fine))

  # re-order columns by celltype --> status
  column_order = order(df_Extraction_sel$CellType, df_Extraction_sel$Status)
  
  # Create column annotation for CellType
  anno_celltype_status = HeatmapAnnotation(CellType = ht_meta$CellType,
                                           Status = ht_meta$Status,
                                           annotation_name_gp = gpar(fontsize = 16), 
                                           border = TRUE,
                                           simple_anno_size = unit(0.75, "cm"), 
                                           annotation_legend_param = list(title_gp = gpar(fontsize = 16), 
                                                                          labels_gp = gpar(fontsize = 16),
                                                                          nrow = 3, by_row = TRUE),
                                           col = list(CellType = CellType_color,
                                                      Status = c("Prior" = "#0F9D58", 
                                                                 "Lesion" = "#DB4437",  
                                                                 "Post" = "#F4B400" )))
  
  # Create a color mapping for the heatmap
  c(min(ht_data), max(ht_data))
  col_fun = colorRamp2(c(0, 0.01, 7), c("gray","blue", "red"))
  
  # Generate the heatmap with column annotations
  htmp = Heatmap(ht_data, col = col_fun, row_names_side = "right", 
                 cluster_rows = FALSE, cluster_columns = FALSE, 
                 show_column_names = FALSE, column_order = column_order, 
                 top_annotation = anno_celltype_status,
                 row_split = seq(1, nrow(ht_data)),
                 column_title_gp = gpar(fontsize = 0), 
                 column_split = ht_meta$CellType, 
                 column_title_rot = 45, 
                 row_names_gp = gpar(fontsize = 17),
                 row_title = NULL,
                 row_gap = unit(1, "mm"),
                 border = TRUE, 
                 heatmap_legend_param = list(title = "Expression\nLevel", title_position = "topcenter", direction = "horizontal", 
                                             title_gp = gpar(col = "red", fontsize = 16), 
                                             labels_gp = gpar(col = "black", fontsize = 16)), use_raster = TRUE)
  
  return(
    draw(htmp, heatmap_legend_side="bottom", annotation_legend_side="bottom", merge_legend = TRUE)
  )
}