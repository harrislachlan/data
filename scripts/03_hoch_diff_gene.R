load("hoch_nsc_lineage.RData")

DefaultAssay(hoch_nsc_lineage) <- "RNA"
hoch_nsc_lineage <- NormalizeData(hoch_nsc_lineage, assay = "RNA")
A <- FindMarkers(hoch_nsc_lineage, ident.1 = "nsc", ident.2 = "granule neurons", 
                 min.pct = 0.0, assay = "RNA", logfc.threshold = 0.0)

sig_genes_pvalue <- A$p_val
FDR <- p.adjust(sig_genes_pvalue, "fdr")
A$p_val_adj <- FDR
A <- A %>% filter(pct.1 >= 0.1|pct.2 >= 0.1)

#granule neurons here are a mix of neuroblasts and more mature neurons
write.csv(A, "nscs_neurob.csv")


