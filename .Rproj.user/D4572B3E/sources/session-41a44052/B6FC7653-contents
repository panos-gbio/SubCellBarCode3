# installed in R version 4.4.2

# CRAN packages 
install.packages("renv")
install.packages(c("tidyverse","data.table","mclust","dplyr","igraph",
                   "Hmisc","ggplot2","ggrepel","patchwork","mgcv","reshape2","circlize",
                   "RColorBrewer","httr","Seurat","readxl"))



# Bioconductor installation at 4.4.2 
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.20")



# Packages requiring Bioconductor 
BiocManager::install("SubCellBarCode")
BiocManager::install("DEqMS")
BiocManager::install("ComplexHeatmap")
BiocManager::install("rtracklayer")
BiocManager::install("Gviz")
BiocManager::install("biomaRt")
BiocManager::install("PCAtools")
BiocManager::install("Biostrings")
BiocManager::install("GenomicRanges")
BiocManager::install("clusterProfiler")

R.version
packageStatus()


