setwd("/Users/jhou2/Documents/GitHub/HSV-Viz-Shiny")

load("Cluster_color.Rdata")

StackedBar = function(seurat_object, selectedCellType = NULL, selectedStatus = NULL, selectedSubject = NULL){
  
  # Basic calculation
  Summary_cluster_per_status = seurat_object@meta.data %>%
    count(CellType_Fine, Status) %>%
    rename(Clusters = CellType_Fine) %>%
    # filter(Clusters %in% selectedCellType) %>%
    mutate(
      # Clusters = factor(Clusters, levels = cluster_order),
      Status = factor(Status, levels = c("Prior", "Lesion", "Post"))
    ) %>%
    # arrange(match(Clusters, cluster_order)) %>%
    dplyr::group_by(Status) %>%
    mutate(Total_Cells = sum(n)) %>%
    ungroup() %>%
    mutate(Percentage = n / Total_Cells * 100) %>%
    dplyr::select(-Total_Cells) %>% 
    as.data.frame()
  
  Status_over_cluster_filtered = Summary_cluster_per_status#[Summary_cluster_per_status$Clusters %in% c(selectedCellType), ]
  
  max_space = tapply(Status_over_cluster_filtered$Percentage, Status_over_cluster_filtered$Clusters, max)
  
  # Adjust the data to use the max space for each component
  Status_over_cluster_filtered$AdjustedPercentage = with(Status_over_cluster_filtered, max_space[Clusters])
  
  p = ggplot(Status_over_cluster_filtered, aes(x = Status, y = Percentage, fill = Clusters)) +
    geom_bar(stat = "identity") +
    # if want to show y axis label, lift up axis text y size = 0 limit
    scale_y_continuous(breaks = function(x) {
      c(0.8*round(max(x, na.rm = TRUE)))
    }) +
    geom_text(aes(label = round(Percentage, 1), y = AdjustedPercentage/2), position = position_stack(vjust = 0.5)) +
    theme_minimal() +
    scale_fill_manual(values = Cluster_color) +
    xlab("") +
    ylab("Proportion(%) in total T cells") +
    facet_wrap(~Clusters, ncol = 1, scales = "free_y") +
    theme(panel.border = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          strip.text.x = element_blank(),
          panel.spacing.x = unit(0, "lines"),
          axis.line.x = element_line(), 
          axis.line.y = element_line(), 
          axis.text.x = element_text(size = 14, color = "black", angle = 0, vjust = 0.5, hjust = 0.5),
          axis.text.y = element_text(size = 0, color = "black"), 
          axis.title = element_text(size = 14, color = "black"), 
          legend.position = "none") + 
    guides(fill=guide_legend(ncol = 1, byrow = TRUE, fill = FALSE))
  
  return(p)
  
}



