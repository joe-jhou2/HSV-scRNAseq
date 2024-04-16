#' Generate a Stacked Bar Plot for each subject from a Seurat Object
#'
#' This function processes a Seurat object to generate a stacked bar plot, illustrating the
#' proportion of different cell types within selected subjects across specified statuses.
#' It optionally allows for the selection of specific subjects and customization of cluster colors.
#'
#' @param seurat_object A Seurat object from which statistics are to be computed.
#' @param selectedSubject A character vector specifying the subjects to include in the plot.
#'        If `NULL`, all subjects are included. Default is `NULL`.
#' @param CellType_color A named vector of colors for the clusters to be used in the plot.
#'        Names of the vector should match the cell type names. If `NULL`, default ggplot2
#'        colors are used. Default is `NULL`.
#'
#' @return A ggplot object representing the stacked bar plot of cell type proportions
#'         across selected subjects and statuses.
#'
#' @importFrom ggplot2 ggplot geom_col scale_y_continuous xlab ylab facet_grid scale_fill_manual theme_bw theme
#' @importFrom ggplot2 panel_border panel_grid.major panel_grid.minor panel_spacing strip_background strip_text
#' @importFrom ggplot2 axis_text axis_title axis_ticks axis_line.x legend_position legend_title legend_text guides
#' @importFrom ggplot2 coord_flip scale_fill_manual guide_legend
#' @importFrom tidyr complete
#' @importFrom dplyr filter mutate
#' @importFrom forcats fct_rev
#' @import scales percent
#'
computeStackBarPlot_Subject = function(seurat_object, selectedSubject, CellType_color = NULL){
  
  # calculate first
  res_list = computeStat(seurat_object)
  
  # expand the data to include all subjects with entry and lesion for each cluster
  expanded_data = res_list$Summary_cluster_per_sample %>%
    tidyr::complete(CellType, Subject, Status, fill = list(Freq = 0)) %>%
    filter(if (any(selectedSubject == "All")) TRUE
           else Subject %in% selectedSubject) %>%
    mutate(Subject = fct_rev(factor(Subject, levels = c("Subject1", "Subject2", "Subject3", "Subject5", "Subject6", 
                                                        "Subject7", "Subject8", "Subject9", "Subject10", "Subject11","Subject12", 
                                                        "Subject13", "Subject14", "Subject15", "Subject16", "Subject17", "Subject18"))),
           Status = factor(Status,  levels = c("Prior", "Lesion", "Post")))
  
  number_of_subjects = if(any(selectedSubject == "All")) length(unique(seurat_object$Subject)) else length(selectedSubject)
  scaling_factor = 1 # Adjust this based on your needs and experimentation
  calculated_width = scaling_factor / number_of_subjects
  
  
  p = ggplot(expanded_data, aes(x = Subject, y = Percentage, fill = CellType)) +
    geom_col(position = position_fill(reverse = TRUE), color = "black") + # re-order on the stacked bar
    scale_y_continuous(labels = scales::percent)  +
    xlab("") +
    ylab("Proportion") +
    facet_grid(.~ Status) +
    scale_fill_manual(values = CellType_color) +
    theme_bw(base_size = 12) +
    theme(panel.border = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.spacing = unit(1.5, "lines"),  # Adjust spacing between panels
          strip.background = element_blank(),  # Remove strip background
          strip.text = element_text(size = 16, face = "bold"), 
          axis.text = element_text(size = 12, color = "black"),
          axis.title = element_text(size = 16, color = "black"),
          axis.ticks = element_blank(), 
          axis.line.x = element_line(), 
          legend.position = "bottom",
          legend.title = element_blank(),
          legend.text = element_text(size = 14)) + 
    guides(fill = guide_legend(ncol = 8, byrow = TRUE)) + 
    coord_flip()
  
  return(p)
}
