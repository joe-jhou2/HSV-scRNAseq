
#' Feature Extraction from Seurat Object
#'
#' This function extracts features, aka, genes, based on specified clusters and 
#' features from a Seurat object. It returns a data frame with the selected 
#' features for each cell, along with their cluster and status.
#'
#' @param seurat_object A Seurat object containing single-cell RNA-seq data.
#' @param selectedCellType A character of clusters for which to extract features.
#' @param selectedFeature A character of features/genes to extract.
#'
#' @return A data frame of extracted features, with columns for each of the 
#' selected features, CellType, and other group/category info.
#' The result will look similar to the table below:
#' \dontrun{
#'                         TNF     IFNG     PRF1      orig.ident  Subject CellType Status
#' AAACCTGAGATACACA-1_10   0 3.656766 0.000000 Subject5_Lesion Subject5  CD8 TRM Lesion
#' AAACCTGCAACTGCTA-1_5    0 0.000000 0.000000 Subject2_Lesion Subject2    ISG T Lesion
#' }

feature_extraction = function(seurat_object, selectedCellType, selectedFeature){
  
  # subset data
  selected_object_select = tryCatch({
    # If 'All' is one of the selected clusters, return the entire Seurat object
    if("All" %in% selectedCellType) {
      seurat_object
    } else {
      # Otherwise, subset the Seurat object based on the selected clusters
      subset(seurat_object, subset = Idents(seurat_object) %in% selectedCellType)
    }
  }, error = function(e) {
    # Optionally print the error message for debugging purposes
    message("Error during subsetting: ", e$message)
    # Return NULL to indicate an error occurred during subsetting
    NULL
  })
  
  # Fetch expression data and meta data
  exp_df = FetchData(object = selected_object_select, vars = selectedFeature)
  meta_df = selected_object_select@meta.data[, c('orig.ident', 'Subject', 'CellType_Fine', 'Status')]
  
  # Stitch expression and meta data
  feature_df = data.frame(merge(exp_df, meta_df, by = 'row.names', all = TRUE))
  
  # clean up row names
  rownames(feature_df) = feature_df$Row.names
  feature_df = feature_df[,-1]
  feature_df = feature_df %>% rename("CellType" = "CellType_Fine")
  
  # assign levels
  feature_df$CellType = factor(feature_df$CellType, levels = if ("All" %in% selectedCellType) unique(feature_df$CellType) else selectedCellType)
  feature_df$Status = factor(feature_df$Status, levels = c("Prior", "Lesion", "Post"))
  
  return(feature_df)
}


#' Calculate Percentage of Cells Expressing Selected Genes within Specific Clusters and Statuses
#'
#' This function returns a data frame with the percentage of cells expressing 
#' each of the selected features/genes for the given clusters and statuses. 
#' The data frame includes the percentage value of cell expression per feature,
#' grouped by Subject, CellType, and Status.
#'
#' @param feature_df A data frame extracted from a Seurat object or similar, 
#' containing columns for Subject, CellType, Status, and gene expression data.
#' @param selectedFeature A character vector of features/genes to calculate the 
#' expression percentage for.
#' @param selectedCellType A character vector of cell types/clusters of interest.
#' @param selectedStatus A character vector of statuses/time points of interest.
#'
#' @return A data frame where each row corresponds to a unique combination of 
#' Subject, CellType, and Status, and includes columns for each of the 
#' selected genes indicating the percentage of cells expressing the gene.
#' 
#' @examples
#' \dontrun{
#' Subject   CellType Status           TNF             IFNG            PRF1
#' <chr>     <fct>    <fct>           <dbl>           <dbl>           <dbl>
#' Subject1  ISG T    Prior            30.8            38.5            38.5
#' Subject1  ISG T    Lesion           50               0               0  
#' Subject1  CD8 TRM  Prior            42.2            53.1            77.1
#' }

calculate_feature_pert = function(feature_df, selectedFeature, selectedCellType, selectedStatus) {
  # Check for necessary columns in feature_df
  required_cols = c("Subject", "CellType", "Status")
  if (!all(required_cols %in% names(feature_df))) {
    stop("feature_df is missing one or more of the required columns: Subject, CellType, Status.")
  }
  
  # Ensure selectedFeature exist in feature_df
  if (!all(selectedFeature %in% names(feature_df))) {
    missing_genes = selectedFeature[!selectedFeature %in% names(feature_df)]
    stop("feature_df is missing the following genes: ", paste(missing_genes, collapse = ", "), ".")
  }
  
  # Process data
  result = feature_df %>%
    filter((any(selectedCellType == "All") | CellType %in% selectedCellType) &
             (any(selectedStatus == "All") | Status %in% selectedStatus)) %>%
    group_by(Subject, CellType, Status) %>%
    summarise(across(all_of(selectedFeature), ~mean(. > 0) * 100, .names = "{.col}")) %>%
    ungroup()
  
  return(result)
}

#' Compute Cluster Statistics for Seurat Object
#'
#' This function calculates statistics across different dimensions (Subject, Status, CellType)
#' of a Seurat object. It generates summaries per subject and status, per cluster and sample,
#' and per cluster and status. Optionally, it allows filtering of the clusters to be included
#' in the analysis based on a provided set of cluster identities.
#'
#' @param seurat_object A Seurat object containing single-cell RNA-seq data.
#'
#' @return A list containing three data frames:
#'         - `Summary_perSubject_perStatus`: Counts of cells per Subject and Status.
#'         - `Summary_cluster_per_sample`: Counts and percentages of cells by CellType and Sample.
#'         - `Summary_cluster_per_status`: Counts and percentages of cells by CellType and Status.
#' @import dplyr
#' @importFrom tidyr separate
#' @importFrom stringr str_extract
#' @export
#'
computeStat = function(seurat_object){
  
  # Conceptually, the proportion of each cluster is fixed, it should be reflected
  # in total population in each subject and time points, should NOT only calculated
  # in some selected clusters. But you can pull out specific cell type afterwards
  
  # seurat_object_select = if (is.null(selectedCellType) | "All" %in% selectedCellType) {
  #   # If 'selectedCellType' is NULL or contains "All", no subsetting by 'idents'
  #   seurat_object
  # } else {
  #   # Subset by the specified cluster identities
  #   subset(seurat_object, idents = selectedCellType)
  # }
  
  seurat_object_select = seurat_object
  
  # Summary1
  Summary_perSubject_perStatus = as.data.frame(seurat_object_select@meta.data %>% group_by(Subject,  Status) %>% tally())
  
  # Summary2 - by Clusters and Sample
  # Custom sort function
  custom_sort = function(samples, suffix_order) {
    ids = str_extract(samples, "^[^_]+_[^_]+")
    suffixes = str_extract(samples, "[^_]+$")
    suffix_order = match(suffixes, suffix_order)
    sorted = samples[order(ids, suffix_order)]
    return(sorted)
  }
  suffix_order = c("Prior", "Lesion", "Post")
  
  Summary_cluster_per_sample = seurat_object_select@meta.data %>%
    count(CellType_Fine, orig.ident) %>%
    rename(CellType = CellType_Fine, Sample = orig.ident) %>%
    # filter(if (length(selectedCellType) == 0) TRUE
    #        else CellType %in% selectedCellType) %>%
    tidyr::separate(Sample, into = c("Subject", "Status"), sep = "_", remove = FALSE) %>%
    mutate(
      Status = case_when(Status == "Entry" ~ "Prior",
                         Status == "8WPH" ~ "Post",
                         TRUE ~ Status)
    ) %>%
    mutate(
      CellType = factor(CellType, levels = levels(seurat_object_select$CellType_Fine)),
      Sample = factor(Sample, levels = custom_sort(unique(Sample), suffix_order)),
      Status = factor(Status, levels = c("Prior", "Lesion", "Post"))
    ) %>%
    # arrange(match(CellType, selectedCellType)) %>%
    dplyr::group_by(Sample) %>%
    mutate(Total_Cells = sum(n)) %>%
    ungroup() %>%
    mutate(Percentage = n / Total_Cells * 100) %>%
    dplyr::select(-Total_Cells) %>%
    select(Sample, Subject, Status, CellType, n, Percentage) %>% 
    as.data.frame()
  
  # Summary3 by Cluster and Status
  Summary_cluster_per_status = seurat_object_select@meta.data %>%
    count(CellType_Fine, Status) %>%
    rename(CellType = CellType_Fine) %>%
    # filter(if (length(selectedCellType) == 0) TRUE
    #        else CellType %in% selectedCellType) %>%
    mutate(
      CellType = factor(CellType, levels = levels(seurat_object_select$CellType_Fine)),
      Status = factor(Status, levels = c("Prior", "Lesion", "Post"))
    ) %>%
    # arrange(match(CellType, selectedCellType)) %>%
    dplyr::group_by(Status) %>%
    mutate(Total_Cells = sum(n)) %>%
    ungroup() %>%
    mutate(Percentage = n / Total_Cells * 100) %>%
    dplyr::select(-Total_Cells) %>% 
    as.data.frame()
  
  return(list(Summary_perSubject_perStatus = Summary_perSubject_perStatus,
              Summary_cluster_per_sample = Summary_cluster_per_sample,
              Summary_cluster_per_status = Summary_cluster_per_status))
}
