Figure 3
================
2024-07-03

- [R Markdown](#r-markdown)
- [Load libraries](#load-libraries)
- [Load the RNA-Seq data](#load-the-rna-seq-data)
- [Figure 3-A: Data Preparation for UpSet Plot for Upregulated
  genes](#figure-3-a-data-preparation-for-upset-plot-for-upregulated-genes)
- [Figure 3-A: UpSet plot for Upregulated
  genes](#figure-3-a-upset-plot-for-upregulated-genes)
- [Figure 3-B: Violin plots of SNCA in iNL2 and LRRK2 in
  oligodendroglial
  lineage](#figure-3-b-violin-plots-of-snca-in-inl2-and-lrrk2-in-oligodendroglial-lineage)
- [Figure 3-C: GO terms of upregulated
  DEGs](#figure-3-c-go-terms-of-upregulated-degs)
- [Figure 3-C: GO Heatmap for upregulated
  DEGs](#figure-3-c-go-heatmap-for-upregulated-degs)

------------------------------------------------------------------------

## R Markdown

This is an R Markdown document. Markdown is a simple formatting syntax
for authoring HTML, PDF, and MS Word documents. For more details on
using R Markdown see <http://rmarkdown.rstudio.com>.

When you click the **Knit** button a document will be generated that
includes both content as well as the output of any embedded R code
chunks within the document. You can embed an R code chunk like this:

## Load libraries

``` r
.libPaths( c( "/data/Common_Folder/R/Single_cell_packages/", .libPaths()) )
library("SCpubr")
library(Seurat)
library(patchwork)
library(gridExtra)
library(ggplot2)
library(future)
library(dplyr)
library(SeuratData)
library(glmGamPoi, lib.loc = "/data/nasser/R/packages_nasser/")
library(harmony)
library(enrichR)
library(gprofiler2)
library(data.table) 
library(speckle, lib.loc = "/data/Common_Folder/R/Single_cell_packages/")
library(ggpubr)
library(tibble)

library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(viridis)
library(UpSetR)

#devtools::install_github('immunogenomics/presto')
```

## Load the RNA-Seq data

``` r
pd = readRDS("/data/nasser/Manuscript/Strict_threshold/processedobject/ODC35_woClus8_subclust3_res0.15_NK")
```

<!-- ## Supplementary table S5 - -->
<!-- ```{r} -->
<!-- # Load necessary libraries -->
<!-- library(dplyr) -->
<!-- library(readr) -->
<!-- library(openxlsx) -->
<!-- pathto.outsupTable <- "/data/nasser/Manuscript/table/Supplementary Tables/" -->
<!-- # Define the path to the DEG files -->
<!-- #pathto.outTable <- "/data/nasser/Manuscript/table/Figure3/DEGs/" -->
<!-- # Define clusters to process -->
<!-- clusters <- c("iPPC_0", "iPPC_1", "iPPC_2", "iNL2", "iOPC", "iCEP", "iNL1", "iRGC", "iODC", "iINPC", "All") -->
<!-- # Initialize a new workbook -->
<!-- wb <- createWorkbook() -->
<!-- # Loop through each cluster to read the DEG files and add them to separate sheets -->
<!-- for (cluster in clusters) { -->
<!--   file_path <- paste0(pathto.outTable, cluster, "_DEG_LRRK2.csv") -->
<!--   if (file.exists(file_path)) { -->
<!--     deg_data <- read_csv(file_path) -->
<!--     addWorksheet(wb, cluster) -->
<!--     writeData(wb, cluster, deg_data) -->
<!--     header_style <- createStyle(fontSize = 11, fontName = "Arial", textDecoration = "bold") -->
<!--     addStyle(wb, cluster, header_style, rows = 1, cols = 1:ncol(deg_data), gridExpand = TRUE) -->
<!--   } -->
<!-- } -->
<!-- # Save the workbook -->
<!-- #output_file <- "/mnt/data/TableS5_combined.xlsx" -->
<!-- #saveWorkbook(wb, file = paste0(pathto.outsupTable,"TableS5_combined.xlsx"), overwrite = TRUE) -->
<!-- ``` -->
<!-- ## Supplementary table S6 - -->
<!-- GO terms downloaded from EnrichR are put together in a Table for Supplementary Table S6 -->
<!-- ```{r} -->
<!-- # Define the path to your directory containing the .txt files -->
<!-- path_to_files <- "/data/nasser/Manuscript/DEGs/GO terms/" -->
<!-- # List all the .txt files in the directory -->
<!-- txt_files <- list.files(path = path_to_files, pattern = "*.txt", full.names = TRUE) -->
<!-- # Create a new workbook -->
<!-- wb <- createWorkbook() -->
<!-- # Loop through each file and add it as a sheet in the workbook -->
<!-- for (file in txt_files) { -->
<!--   # Read the file -->
<!--   data <- read_tsv(file) -->
<!--   # Extract the base name of the file and modify it to create the sheet name -->
<!--   file_name <- basename(file) -->
<!--   sheet_name <- sub("GO_Biological_Process_2023_", "", file_name) -->
<!--   sheet_name <- sub(".txt", "", sheet_name) -->
<!--   # Add a sheet to the workbook with the modified sheet name -->
<!--   addWorksheet(wb, sheet_name) -->
<!--   # Write the data to the sheet -->
<!--   writeData(wb, sheet_name, data) -->
<!-- } -->
<!-- #pathto.outTable <- "/data/nasser/Manuscript/table/" -->
<!-- # Save the workbook to a file -->
<!-- #saveWorkbook(wb,  file = paste0(pathto.outTable,"S6.xlsx"), overwrite = TRUE) -->
<!-- ``` -->

## Figure 3-A: Data Preparation for UpSet Plot for Upregulated genes

``` r
# Define clusters and path to DEG files
clusters <- c("iPPC_0", "iPPC_1", "iPPC_2", "iNL2", "iOPC", "iCEP", "iNL1", "iRGC", "iODC", "iINPC")
pathto.inTable <- "/data/nasser/Manuscript/Strict_threshold/table/Figure3/DEGs_up/"

# Initialize an empty list to store DEGs for each cluster
deg_lists <- list()

# Loop through each cluster to read the DEG files
for (cluster in clusters) {
    file_path <- paste0(pathto.inTable, cluster, "_up.csv")
    if (file.exists(file_path)) {
        deg_df <- read.csv(file_path, row.names = 1)
        # Assuming the gene names are row names
        deg_lists[[cluster]] <- rownames(deg_df)
    }
}

#upset(fromList(deg_lists), order.by = "freq")

# Create a combination matrix using ComplexHeatmap's make_comb_mat function
comb_mat <- make_comb_mat(deg_lists)
```

## Figure 3-A: UpSet plot for Upregulated genes

``` r
#Upset plots data prep

# Define clusters and path to DEG files
clusters <- c("iPPC_0", "iPPC_1", "iPPC_2", "iNL2", "iOPC", "iCEP", "iNL1", "iRGC", "iODC", "iINPC")
pathto.inTable <- "/data/nasser/Manuscript/Strict_threshold/table/Figure3/DEGs_up/"

# Initialize an empty list to store DEGs for each cluster
deg_lists <- list()

# Loop through each cluster to read the DEG files
for (cluster in clusters) {
    file_path <- paste0(pathto.inTable, cluster, "_up.csv")
    if (file.exists(file_path)) {
        deg_df <- read.csv(file_path, header = FALSE)
        # Extract the gene names from the first column
        gene_names <- deg_df[, 1]
        deg_lists[[cluster]] <- gene_names
    }
}
 
# Create a combination matrix using ComplexHeatmap's make_comb_mat function
comb_mat <- make_comb_mat(deg_lists)

comb_mat <- comb_mat[comb_size(comb_mat) >= 5]
desired_celltype_order <- c('iODC', 'iOPC', 'iPPC_2', 'iPPC_1', 'iPPC_0', 'iCEP', 'iNL2', 'iNL1', 'iINPC', 'iRGC')

#pathto.outPlots = "/data/nasser/Manuscript/plots/figure3/"
#png(paste0(pathto.outPlots,"UpSet_up_new.png"), width=2500, height=1500, res = 300)
ht <- UpSet(comb_mat, top_annotation = upset_top_annotation(comb_mat, add_numbers = TRUE),
    left_annotation = upset_left_annotation(comb_mat, add_numbers = TRUE),
            set_order = desired_celltype_order,
            comb_order = order(comb_size(comb_mat), decreasing = TRUE))
ht  
```

![](Figure3_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

``` r
#dev.off()
```

## Figure 3-B: Violin plots of SNCA in iNL2 and LRRK2 in oligodendroglial lineage

``` r
#Vln plots
Idents(pd) = "CellType"
#pathto.outPlots= "/data/nasser/Manuscript/plots/figure3/"
celltype_subset <- subset(pd, idents = c("iOPC","iPPC_2", "iPPC_0", "iPPC_1", "iODC"))
Idents(celltype_subset) <- "Mutation"
#png(paste0(pathto.outPlots,"vlnplot_iOPC_LRRK2_PDGFRA_up_regulated.png"), width=2000, height=1100, res = 300)
VlnPlot(celltype_subset, features = "LRRK2")
```

![](Figure3_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

``` r
#dev.off()

celltype_subset <- subset(pd, idents = "iNL2")
Idents(celltype_subset) <- "Mutation"
#png(paste0(pathto.outPlots,"vlnplot_iNL2_SNCA_up_regulated.png"), width=2000, height=1100, res = 300)
VlnPlot(celltype_subset, features = c("SNCA") )
```

![](Figure3_files/figure-gfm/unnamed-chunk-8-2.png)<!-- -->

``` r
#dev.off()
```

## Figure 3-C: GO terms of upregulated DEGs

``` r
# Define the path to the directory containing the files
path <- "/data/nasser/Manuscript/Strict_threshold/table/Figure3/GOs/GO_ODC35/"

# List all files in the directory that match the pattern for downregulated genes
files_up <- list.files(path, pattern = "up.txt$", full.names = TRUE)

# Initialize an empty list to store data frames
go_data_list <- list()
top_terms_list <- list()

# Read each file and store in the list
for (file in files_up) {
  cell_type <- gsub(".*\\/|_up.txt", "", file) # Extract cell type from file name
  go_data <- read.delim(file, header = TRUE, sep = "\t") # Adjust if necessary
  
  # Add a column for cell type
  go_data$CellType <- cell_type
  
  # Append the full data to the list
  go_data_list[[cell_type]] <- go_data
  
  # Select top 5 terms by combined score
  top_go_data <- go_data %>%
    arrange(P.value) %>%
    slice_head(n = 5)
  
  # Append to the top terms list
  top_terms_list[[cell_type]] <- top_go_data
}

# Combine all top terms data frames into one
top_go_data_combined <- bind_rows(top_terms_list)

# Define the desired cell type order
desired_celltype_order <- c('GO_iODC', 'GO_iOPC', 'GO_iPPC_2', 'GO_iPPC_1', 'GO_iPPC_0', 'GO_iCEP', 'GO_iNL2', 'GO_iNL1', 'GO_iINPC', 'GO_iRGC')


# Arrange the dataframe with the specified order of CellType column
top_go_data_combined <- top_go_data_combined %>%
  mutate(CellType = factor(CellType, levels = desired_celltype_order)) %>%
  arrange(CellType)

# Combine all full data frames into one
full_go_data_combined <- bind_rows(go_data_list)

# Create a data frame ensuring all top terms are represented across cell types
heatmap_data <- expand.grid(Term = unique(top_go_data_combined$Term), CellType = unique(full_go_data_combined$CellType)) %>%
  left_join(full_go_data_combined %>% select(Term, Combined.Score, CellType), by = c("Term", "CellType")) %>%
  spread(key = CellType, value = Combined.Score)

# Ensure the order of Terms is preserved
heatmap_data <- heatmap_data %>%
  mutate(Term = factor(Term, levels = unique(top_go_data_combined$Term))) %>%
  arrange(Term)

# Re-order columns in heatmap_data
heatmap_data <- heatmap_data %>%
  select(Term, all_of(desired_celltype_order))

# Convert to a matrix for the heatmap
heatmap_matrix <- as.matrix(heatmap_data %>% select(-Term))
rownames(heatmap_matrix) <- heatmap_data$Term

# Replace NA values with a sensible default, such as a value slightly above the threshold for significance
heatmap_matrix[is.na(heatmap_matrix)] <- 0
```

## Figure 3-C: GO Heatmap for upregulated DEGs

``` r
# Cap the maximum value to a specified threshold
scaled_mat = t(scale(t(heatmap_matrix)))
#scaled_mat[is.na(scaled_mat)] <- 0



#pathto.outPlots = "/data/nasser/Manuscript/plots/figure3/"
#png(paste0(pathto.outPlots,"GO_up_heatmap_scaled_top5_pvalue_all.png"), width=6000, height=4000, res = 300)
# Generate the heatmap with capped values
Heatmap(scaled_mat,
        row_names_max_width = max_text_width(rownames(scaled_mat)),cluster_columns = FALSE, cluster_rows = FALSE,
        border_gp = gpar(col = "black"),
        rect_gp = gpar(col = "white") ,
        name = "Capped Combined Score",
        row_names_gp = gpar(fontsize = 8),
        column_names_gp = gpar(fontsize = 10),
        heatmap_legend_param = list(title = "Combined Score"),
        column_title = "CellType", 
        row_title = "Gene Ontology",
        col = viridis(5) ,
        row_names_rot = 35
)
```

![](Figure3_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

``` r
#dev.off()
```
