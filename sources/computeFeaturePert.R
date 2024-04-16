#' Visualize Cytokine Percentage in Cell Types and Subjects
#'
#' @param seurat_object A Seurat object containing single-cell RNA sequencing data.
#' @param selectedCellType A vector of cell types (clusters) to include in the analysis.
#' @param selectedFeature A vector of genes (features) to include in the analysis, typically cytokines.
#' @param CellType_color A named vector of colors corresponding to each cell type for plotting.
#'
#' @return A ggplot object visualizing the percentage of selected genes within specified cell types across different statuses.
#'
#' @description This function processes a Seurat object to extract expression data for specified genes and cell types,
#'              calculates the percentage of these genes in the cell types across different conditions (Prior, Lesion, Post),
#'              and visualizes the data using a jitter plot with error bars and lines indicating means. The visualization is faceted
#'              by gene and cell type, showing the variability and mean percentage of gene expression within each cell type across
#'              different statuses.

computeFeaturePert = function(seurat_object, selectedCellType, selectedSubject, selectedFeature, CellType_color = NULL){
  
  # Extract features
  exp_meta_df = feature_extraction(seurat_object, selectedCellType = "All", selectedFeature)

  # Calculate gene% in cell type and subject
  exp_pert_df = calculate_feature_pert(feature_df = exp_meta_df %>% filter((any(selectedSubject == "All") | Subject %in% selectedSubject)), 
                                       selectedFeature = selectedFeature, 
                                       selectedCellType = selectedCellType, 
                                       selectedStatus = c("Prior", "Lesion", "Post"))
  
  # reshape data
  exp_pert_df_long = exp_pert_df %>% pivot_longer(cols = !c("Subject", "CellType", "Status"), names_to = "Gene", values_to = "Percentage")
  
  # Ensure factors have the correct levels for plotting
  exp_pert_df_long$Status = factor(exp_pert_df_long$Status, levels = c("Prior", "Lesion", "Post"))
  exp_pert_df_long$CellType = factor(exp_pert_df_long$CellType, levels = levels(seurat_object$CellType_Fine))
  exp_pert_df_long$Gene = factor(exp_pert_df_long$Gene, levels = selectedFeature) 
  
  exp_pert_df_long = exp_pert_df_long %>% droplevels()
  
  # plotting
  g = ggplot(exp_pert_df_long, aes(x = Status, y = Percentage, group = Status)) + 
    geom_jitter(shape = 21, aes(fill = Status), color = "black", size = 2, alpha = 0.8,width = 0.2) +
    stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), geom = "errorbar", color = "#5C5C5C", width = 0.2) +
    stat_summary(fun = mean, geom = "line", aes(group = CellType, color = CellType), size = 1.5) +
    stat_summary(fun = mean, geom = "pointrange", color = "black") +
    stat_summary(fun = mean, geom = "line", color = "black", aes(group = 1)) +
    scale_y_continuous(limits = c(-5, 105),position = "left") + 
    scale_fill_manual(values = c("Prior" = "#0F9D58", 
                                 "Lesion" = "#DB4437",  
                                 "Post" = "#F4B400")) +
    scale_color_manual(values = CellType_color) +
    labs(x = NULL, y = NULL) +
    # stat_compare_means(aes(label = "p.signif"),
    #                    size = 6, hide.ns = TRUE, label.y = c(75, 85),
    #                    comparisons = list(c("Lesion", "Prior"))) +
    facet_grid2(Gene ~ CellType) +
    theme_bw() +
    theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1), 
          panel.spacing = unit(0.05, "lines"),
          panel.grid = element_blank(),
          strip.background = element_blank(),
          strip.text.x = element_text(size = 14, angle = 0, face = "plain", hjust = 0.5),
          strip.text.y = element_text(size = 14, angle = 0, face = "plain", hjust = 0.5),
          axis.title.x = element_text(size = 0, color = 'black'),
          axis.title.y = element_text(size = 14, color = 'black'),
          axis.text.x = element_text(size = 14, color = 'black', angle = 90,  hjust = 1, vjust = 0.5),
          axis.text.y = element_text(size = 14, color = 'black'),
          legend.position = "none") +
    labs(y = 'Cyto % in CellType and Subject')
  
  return(g)
}


