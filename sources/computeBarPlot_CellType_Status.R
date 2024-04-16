#' Compute Bar Plots of Cell Types Across Different Statuses or Subjects
#'
#' This function creates two bar plots based on the provided Seurat object.
#' One bar plot shows the count of cells within each cell type for selected subjects.
#' The other shows the percentage of cells within each cell type for selected statuses.
#'
#' @param seurat_object A Seurat object containing single-cell RNA-seq data.
#' @param selectedSubject A vector of selected subject IDs. If "All" or unspecified,
#'        the function uses all subjects.
#'
#' @return A list containing two ggplot objects: 'CountBar' for counts of cells,
#'         and 'PertBar' for the percentage of cells within each status.
#'

computeBarPlot_CellType_Status = function(seurat_object, selectedSubject){
  # calculate and pull results
  res_list = computeStat(seurat_object)
  
  # sort of all data or selected subjects
  if (any(selectedSubject == "All")){
    data_plot = res_list$Summary_cluster_per_status
  } else {
    data_plot = res_list$Summary_cluster_per_sample %>%
      tidyr::complete(CellType, Subject, Status, fill = list(Freq = 0)) %>%
      filter(if (any(selectedSubject == "All")) TRUE
             else Subject %in% selectedSubject) %>%
      mutate(Subject = fct_rev(factor(Subject, levels = c("Subject1", "Subject2", "Subject3", "Subject5", "Subject6", 
                                                          "Subject7", "Subject8", "Subject9", "Subject10", "Subject11","Subject12", 
                                                          "Subject13", "Subject14", "Subject15", "Subject16", "Subject17", "Subject18"))),
             Status = factor(Status,  levels = c("Prior", "Lesion", "Post"))) 
  }
  
  # define theme
  theme_BarCellType2Status = theme(panel.background = element_rect(color = "black", fill = NA),
                                   panel.grid.major = element_blank(),
                                   panel.grid.minor = element_blank(),
                                   strip.text = element_blank(),
                                   plot.title = element_text(size = 16, hjust = 0.5),
                                   axis.title.x = element_blank(),
                                   axis.title.y = element_text(size = 14, color = "black"),
                                   axis.text.x = element_text(size = 12, color = "black", angle = 0, hjust = 0.5, vjust = 1),
                                   axis.text.y = element_text(size = 12, color = "black"),
                                   legend.title = element_blank(), 
                                   legend.text = element_text(size = 14, color = "black"), 
                                   legend.position = "bottom") 
  
  # count plot
  p1 = ggplot(data_plot, aes(x = CellType, y = n, fill = Status)) +
    geom_bar(stat = "identity", position = position_dodge(width = 1)) +
    scale_fill_manual(values = c("Prior" = "#0F9D58", "Lesion" = "#DB4437", "Post" =  "#F4B400")) +  # Use a distinct color palette for clarity
    xlab("Cell Type") +
    ylab("Count") +
    ggtitle("Counts for each Status and CellType") +
    facet_wrap(~ CellType, scales = "free", ncol = 6) + 
    theme_minimal() +
    theme_BarCellType2Status
  
  # pert plot
  p2 = ggplot(data_plot, aes(x = CellType, y = Percentage, fill = Status)) +
    geom_bar(stat = "identity", position = position_dodge(width = 1)) +
    scale_fill_manual(values = c("Prior" = "#0F9D58", "Lesion" = "#DB4437", "Post" =  "#F4B400")) +  # Use a distinct color palette for clarity
    xlab("Cell Type") +
    ylab("Percentage(%)") +
    ggtitle("Percentage of each CellType in each Status") +
    facet_wrap(~ CellType, scales = "free", ncol = 6) + 
    theme_minimal() +
    theme_BarCellType2Status
  
  return(list(CountBar = p1,
              PertBar = p2))
}




