<h1 align="center">HSV+ Skin Biopsy scRNAseq Data Visualization</h1>Author: Joe Hou (Fred Hutch)

# <span style="font-size:15px;"> 
## Overview
This Shiny application provides an interactive interface for exploring 10X scRNAseq data derived from HSV+ skin biopsy samples. As the manuscript is currently under review for publication, the raw and processed data remain undisclosed. This repository currently includes only the code and demo images to demonstrate the capabilities of the app. The app link will be made publicly accessible once the manuscript is published.

## Features and Functionalities

### Interactive Exploration on Cell Type
- **Cell Type and Status Identification:** Users can view all available cell types and subjects, with the flexibility to select specific ones for detailed exploration.(others/App%20showcase/Slide2.png)
- **Dynamic UMAP Visualization:** Initial displays show holistic UMAP results, which can be refined by status. Selecting different subjects and cell types dynamically updates the UMAP to show specific clusters.
- **Cell Type Statistics:** Displays percentages and counts of cell types for selected or all subjects.
- **Cell Type by Subject:** Visualizes the percentage of each cell type within a given sample using stacked bar plots.

### Gene Expression Analysis
- **Feature Exploration:** Allows for the selection of cell types, subjects, and specific genes of interest.
- **Gene Highlighting:** Includes UMAP overlays, dot plots showing expression intensity across cell types, and a heatmap detailing expression levels at the single-cell level annotated by cell type and status.
- **Violin Plots:** Shows average expression levels across cell types and statuses, aiding in understanding temporal dynamics.
- **Feature Percentage Analysis:** Summarizes gene expression at the sample/subject level, indicating the proportion of cells expressing a specific gene within a cell type.

## Gallery
Here are some visual showcases of the application features:
- **Cell Type Exploration**
  ![UMAP Visualization](others/App%20showcase/Slide1.png)
  ![Subject and Cell Type Selector](others/App%20showcase/Slide4.png)

- **Statistical Views**
  ![Cell Type Statistics](others/App%20showcase/Slide5.png)
  ![Cell Type by Subject](others/App%20showcase/Slide6.png)

- **Gene Expression Analysis**
  ![Feature Gene Highlight](others/App%20showcase/Slide8.png)
  ![Expression Heatmap](others/App%20showcase/Slide9.png)
  ![Violin Plot](others/App%20showcase/Slide10.png)

- **Perturbation Analysis**
  ![Feature Perturbation](others/App%20showcase/Slide11.png)
</span># HSViz
