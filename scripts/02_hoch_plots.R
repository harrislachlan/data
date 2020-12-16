load("hoch_nsc_lineage.RData")
cbPalette <- c( "#009E73","#0072B2","#D55E00")

A <- DimPlot(hoch_nsc_lineage, reduction = "umap") + scale_fill_manual(values=cbPalette)
A <- A + theme_minimal() 
A


#Featureplots of different genes
hoch_nsc_lineage <- NormalizeData(hoch_nsc_lineage, assay = "RNA")
DefaultAssay(hoch_nsc_lineage) <- "RNA"

B <- FeaturePlot(hoch_nsc_lineage, features = "Hopx", max.cutoff = 2)
B <- B + theme_minimal()
B
C <- FeaturePlot(hoch_nsc_lineage, features = "Eomes", max.cutoff = 2)
C <- C + theme_minimal()
C
D <- FeaturePlot(hoch_nsc_lineage, features = "Tubb3", max.cutoff = 2)
D <- D + theme_minimal() + scale_fill_manual(values=cbPalette)
D

ggsave(file="plots/classification.png", A, width = 5, height = 3, scale = 1)


E <- FeaturePlot(hoch_nsc_lineage , features = "Bace1", max.cutoff = 2, order = TRUE)
E <- E + theme_minimal() 
G <- FeaturePlot(hoch_nsc_lineage , features = "Fabp7", max.cutoff = 2, order = TRUE)
G <- G + theme_minimal() 
H <- FeaturePlot(hoch_nsc_lineage , features = "Psen1", max.cutoff = 1, order = TRUE)
H <- H + theme_minimal() 
I <- FeaturePlot(hoch_nsc_lineage , features = "Psen2", max.cutoff = 1, order = TRUE)
I <- I + theme_minimal()
J <- FeaturePlot(hoch_nsc_lineage , features = "Psenen", max.cutoff = 4, order = TRUE)
J <- J + theme_minimal() 
K <- FeaturePlot(hoch_nsc_lineage , features = "App", max.cutoff = 10, order = TRUE)
K <- K + theme_minimal() 
L  <- FeaturePlot(hoch_nsc_lineage , features = "Adam10", max.cutoff = 10, order = TRUE)
L <- L + theme_minimal() 
M <- FeaturePlot(hoch_nsc_lineage , features = "Bace1", max.cutoff = 10, order = TRUE)
M <- M + theme_minimal() 

g <- arrangeGrob(E, G, H, I, J, K, L, M, nrow = 2)
ggsave(file="plots/AD_genes.png", g, width = 10, height = 5, scale = 1) 



E <- VlnPlot(hoch_nsc_lineage , features = "Bace1", pt.size = 0.1,sort = "increasing")
E <- E + theme_minimal() + scale_fill_manual(values=cbPalette)
G <- VlnPlot(hoch_nsc_lineage , features = "Fabp7", pt.size = 0.1,sort = "increasing")
G <- G + theme_minimal() + scale_fill_manual(values=cbPalette)
H <- VlnPlot(hoch_nsc_lineage , features = "Psen1", pt.size = 0.1,sort = "decreasing")
H <- H + theme_minimal()+ scale_fill_manual(values=cbPalette)
I <- VlnPlot(hoch_nsc_lineage , features = "Psen2", pt.size = 0.1, sort = "increasing")
I <- I + theme_minimal()+ scale_fill_manual(values=cbPalette)
J <- VlnPlot(hoch_nsc_lineage , features = "Psenen", pt.size = 0.1, sort = "increasing")
J <- J + theme_minimal()+ scale_fill_manual(values=cbPalette)
K <- VlnPlot(hoch_nsc_lineage , features = "App", pt.size = 0.001, sort = "decreasing")
K <- K + theme_minimal() + scale_fill_manual(values=cbPalette)
L <- VlnPlot(hoch_nsc_lineage , features = "Adam10", pt.size = 0.1, sort = "decreasing")
L <- L + theme_minimal()+ scale_fill_manual(values=cbPalette)
M <- VlnPlot(hoch_nsc_lineage , features = "Bace1", pt.size = 0.1, sort = "decreasing")
M <- M + theme_minimal() + scale_fill_manual(values=cbPalette)

g <- arrangeGrob(E, G, H, I, J, K, L, M, nrow = 2)
ggsave(file="plots/AD_genes_violin.png", g, width = 13, height = 4, scale = 1) 

a <- ls()
rm(list = a)
rm(a)
