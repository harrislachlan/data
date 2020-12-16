#  Load packages and define data ----
library(reticulate)
library(umap)
library(Matrix)
library(lubridate)
library(glue)
library(RColorBrewer)
library(gridExtra)
library(gdata)
library(slingshot)
library(clusterProfiler)
library(org.Mm.eg.db)
library(DOSE)
library(gridExtra)
library(cowplot)
library(SummarizedExperiment)
library(Seurat)
library(tidyverse)

#load in dataset of GSE95315 - remove metadata to make readable by seurat
hochgerner <- read_delim("hoch_raw_data/GSE95315_10X_expression_data_v2.tab", 
  delim = '\t')
hochgerner <- hochgerner[3:14547, 1:5455]
hochgerner<- column_to_rownames(hochgerner, var = "cellid")

#read into seurat and analyse
seurat <- CreateSeuratObject(hochgerner)
seurat <- NormalizeData(seurat) %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors %>%
  FindClusters()
seurat <- RunUMAP(seurat, reduction = "pca", dims = 1:30)

#plot data
DimPlot(seurat, reduction = "umap", label = TRUE)

#Subset to only include neurogenic lineage (NSCs, IPCs, Neurons)
hoch_nsc_lineage <- subset(seurat, idents = c("8", "15", "3", "6",
                                              "2", "1", "0", "4"))
#re-analyse in Seurat
hoch_nsc_lineage <- NormalizeData(hoch_nsc_lineage) %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors %>%
  FindClusters()
hoch_nsc_lineage <- RunUMAP(hoch_nsc_lineage, reduction = "pca", dims = 1:30)

#Plot data
DimPlot(hoch_nsc_lineage, reduction = "umap", label = TRUE)

hoch_nsc_lineage <- RenameIdents(hoch_nsc_lineage, 
                               '0' = "granule neurons",
                               "1" = "granule neurons",
                               "2" = "granule neurons",
                               "3" = "granule neurons",
                               "4" = "granule neurons",
                               "5" = "granule neurons",
                               "6" = "granule neurons",
                               "7" = "granule neurons",
                               "8" = "nsc",
                               "9"= "ipcs")

#Plot data
DimPlot(hoch_nsc_lineage, reduction = "umap", label = TRUE)

gdata::keep(hoch_nsc_lineage, sure = TRUE)
save.image(file = "hoch_nsc_lineage.RData")
rm(hoch_nsc_lineage)

