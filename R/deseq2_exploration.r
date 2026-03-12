library(DESeq2)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(enrichplot)
library(clusterProfiler)

# ── Parameters ────────────────────────────────────────────────────────────────
pca_intgroup      <- "Timepoint"          # Column in coldata for PCA colouring
pca_title         <- "PCA Plot"           # PCA plot title
pca_colors        <- c("group1" = "#f8766d", "group2" = "#00bfc4")  # Named by group levels
pca_labels        <- c("group1" = "Group 1", "group2" = "Group 2")  # Legend labels

volcano_title     <- "Volcano Plot"       # Volcano plot title
volcano_lfc_thr   <- 1.0                  # log2FC threshold
volcano_padj_thr  <- 0.05                 # Adjusted p-value threshold
genes_of_interest <- c("GENE1", "GENE2") # Genes to highlight in volcano plot

heatmap_top_n     <- 50                   # Number of top genes for heatmap
heatmap_cutree    <- 4                    # Number of column clusters in heatmap

gsea_lfc_thr      <- 1.0                  # log2FC filter for GSEA
gsea_padj_thr     <- 0.1                  # Adjusted p-value filter for GSEA
gsea_orgdb        <- "org.Hs.eg.db"       # OrgDb (e.g. org.Mm.eg.db for mouse)
# ──────────────────────────────────────────────────────────────────────────────


# ── VST Normalisation ─────────────────────────────────────────────────────────
vsd <- vst(dds, blind = FALSE)


# ── PCA ───────────────────────────────────────────────────────────────────────
pcaData <- plotPCA(vsd, intgroup = pca_intgroup, returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"), 1)

pcaData$SampleNames <- rownames(pcaData)

ggplot(pcaData, aes(x = PC1, y = PC2, color = .data[[pca_intgroup]])) +
  geom_point(size = 3) +
  geom_text(aes(label = SampleNames), vjust = -1, hjust = 0.5, size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle(pca_title) +
  scale_color_manual(values = pca_colors, labels = pca_labels) +
  labs(color = pca_intgroup) +
  coord_fixed()


# ── Sample Distance Heatmap ───────────────────────────────────────────────────
colors_dist      <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
sample_dists     <- dist(t(assay(vsd)))
sample_dist_mat  <- as.matrix(sample_dists)
rownames(sample_dist_mat) <- names(vsd$sizeFactor)
colnames(sample_dist_mat) <- NULL

pheatmap(sample_dist_mat,
         clustering_distance_rows = sample_dists,
         clustering_distance_cols = sample_dists,
         col = colors_dist,
         main = "Sample Distances - Euclidean")


# ── Sample Clustering Dendrogram ──────────────────────────────────────────────
plot(hclust(sample_dists))


# ── MA Plot ───────────────────────────────────────────────────────────────────
resApeT <- lfcShrink(dds, coef = 2, type = "normal", lfcThreshold = 1)
plotMA(resApeT, ylim = c(-3, 3), cex = 0.8)
abline(h = c(-1, 1), col = "dodgerblue", lwd = 2)


# ── Gene Count Plot ───────────────────────────────────────────────────────────
# Replace with gene of interest (Ensembl ID with version)
plotCounts(dds, gene = which(rownames(dds) == "ENSG00000000000.00"), intgroup = "Condition")


# ── Volcano Plot ──────────────────────────────────────────────────────────────
volcano_data <- data.frame(
  log2FoldChange = res$log2FoldChange,
  log10pvalue    = -log10(res$padj),
  Gene           = rownames(res)
)

volcano_data$Significance <- ifelse(
  abs(volcano_data$log2FoldChange) > volcano_lfc_thr & volcano_data$log10pvalue > -log10(volcano_padj_thr),
  "Significant", "Not significant"
)

volcano_data$Highlight <- ifelse(
  volcano_data$Gene %in% genes_of_interest, "Highlighted", volcano_data$Significance
)

vol_colors <- c("Significant" = "#f8766d", "Not significant" = "#181C14", "Highlighted" = "#00FF00")

plot(volcano_data$log2FoldChange, volcano_data$log10pvalue,
     pch = 20, cex = 0.7,
     col = vol_colors[volcano_data$Highlight],
     main = volcano_title,
     xlab = "log2(Fold Change)",
     ylab = "-log10(adjusted p-value)")

abline(v = c(-volcano_lfc_thr, volcano_lfc_thr), h = -log10(volcano_padj_thr), col = "blue", lty = 2)

# Labels for genes of interest (cycling positions to reduce overlap)
pos_cycle <- c(4, 2, 3, 1)
for (i in seq_along(genes_of_interest)) {
  gene <- genes_of_interest[i]
  idx  <- volcano_data$Gene == gene
  if (any(idx)) {
    text(volcano_data$log2FoldChange[idx], volcano_data$log10pvalue[idx],
         labels = gene, col = "black", pos = pos_cycle[(i - 1) %% 4 + 1], font = 1)
  }
}

legend("topright",
       legend = c("Significant", "Not significant", "Highlighted"),
       col = vol_colors, pch = 20, cex = 0.8)


# ── Expression Heatmap (Top N genes by padj) ──────────────────────────────────
top_genes <- head(order(resOrdered$padj), heatmap_top_n)
mat       <- counts(dds)[top_genes, ]
mat       <- t(scale(t(mat)))

anno <- data.frame(Condition = colData(dds)[, "Condition"])
rownames(anno) <- colnames(mat)

heat_plot <- pheatmap(mat,
                      annotation_col = anno,
                      cutree_cols    = heatmap_cutree,
                      main           = paste0("Clustering Heatmap (Top ", heatmap_top_n, " Genes by padj)"),
                      show_rownames  = FALSE)

ggsave("Clustering_Heatmap.pdf", plot = heat_plot, device = "pdf", width = 12, height = 9)


# ── GSEA (Gene Ontology - Biological Process) ─────────────────────────────────
for (k in seq_along(res_list)) {

  res_filt <- as.data.frame(res_list[[k]]) %>%
    dplyr::filter(abs(log2FoldChange) >= gsea_lfc_thr & padj <= gsea_padj_thr) %>%
    tibble::rownames_to_column(var = "Geneid") %>%
    left_join(annot, by = "Geneid")

  if (nrow(res_filt) > 10) {

    print(names(res_list[k]))

    gene_list <- res_filt$log2FoldChange
    names(gene_list) <- res_filt$Gene.name
    gene_list <- sort(gene_list[names(gene_list) != ""], decreasing = TRUE)

    gsea <- gseGO(geneList     = gene_list,
                  ont          = "BP",
                  keyType      = "SYMBOL",
                  minGSSize    = 5,
                  maxGSSize    = 800,
                  pvalueCutoff = gsea_padj_thr,
                  verbose      = TRUE,
                  OrgDb        = gsea_orgdb,
                  pAdjustMethod = "fdr")

    if (nrow(gsea@result) > 0) {

      dplot <- enrichplot::dotplot(gsea, showCategory = 10, split = ".sign") +
        facet_grid(. ~ .sign)

      write.table(gsea@result,
                  file      = paste0("gsea_", names(res_list[k]), ".txt"),
                  sep       = "\t", quote = FALSE, row.names = FALSE)

      ggsave(filename = paste0("dotplot_", names(res_list[k]), ".pdf"),
             plot = dplot, device = "pdf")
    }
  }
}
