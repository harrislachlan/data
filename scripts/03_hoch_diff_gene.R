load("hoch_nsc_lineage.RData")

hoch_nsc_lineage <- NormalizeData(hoch_nsc_lineage)
DefaultAssay(hoch_nsc_lineage) <- "RNA"
A <- FindMarkers(hoch_nsc_lineage, ident.1 = "nsc", min.diff.pct = 0.0, ident.2 = "granule neurons", 
                 min.pct = 0.1, assay = "RNA", logfc.threshold = 0.25)

#determine no genes to correct p-value
length(rownames(hoch_nsc_lineage@assays$RNA@data))

sig_genes_pvalue <- A$p_val
FDR <- p.adjust(sig_genes_pvalue, "fdr", n = 14545)
A$p_val_adj <- FDR

#granule neurons here are a mix of neuroblasts and more mature neurons
write.csv(A, "nscs_neurob.csv")

