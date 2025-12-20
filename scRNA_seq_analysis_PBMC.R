################################################################################
# SINGLE-CELL RNA-SEQ ANALYSIS: PBMC 3K DATASET
################################################################################
# 
# Author: Akshayata Naidu, PhD
# Date: December 2024
# GitHub: https://github.com/AkshNaidu
# 
# Description:
# Complete end-to-end analysis of peripheral blood mononuclear cells (PBMCs)
# single-cell RNA-seq data from 10X Genomics (3k cells). Workflow includes
# quality control, normalization, dimensionality reduction, clustering, cell
# type annotation, and differential expression analysis.
#
# Dataset: PBMC 3k from 10X Genomics
# Technology: Chromium Single Cell 3' v1
# Cells: ~2,700 high-quality PBMCs
# URL: https://cf.10xgenomics.com/samples/cell/pbmc3k/
#
################################################################################

# Load required libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

# Set working directory
setwd("~/scrna-analysis")

################################################################################
# 1. DATA LOADING
################################################################################

cat("\n=== LOADING DATA ===\n")

# Load 10X data
pbmc.data <- Read10X(data.dir = "data/pbmc3k/filtered_gene_bc_matrices/hg19/")

# Create Seurat object
# min.cells = 3: Include features detected in at least 3 cells
# min.features = 200: Include cells with at least 200 detected features
pbmc <- CreateSeuratObject(counts = pbmc.data, 
                           project = "pbmc3k", 
                           min.cells = 3, 
                           min.features = 200)

cat("Initial dataset dimensions:\n")
cat("Genes:", nrow(pbmc), "\n")
cat("Cells:", ncol(pbmc), "\n")

################################################################################
# 2. QUALITY CONTROL
################################################################################

cat("\n=== QUALITY CONTROL ===\n")

# Calculate QC metrics
# Mitochondrial percentage - indicator of cell stress/dying cells
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Ribosomal percentage
pbmc[["percent.ribo"]] <- PercentageFeatureSet(pbmc, pattern = "^RP[SL]")

# Visualize QC metrics BEFORE filtering
p1 <- VlnPlot(pbmc, 
              features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
              ncol = 3, 
              pt.size = 0.1) +
  plot_annotation(title = "QC Metrics Before Filtering")

ggsave("figures/01_qc_before_filtering.png", p1, width = 14, height = 5, dpi = 300)

# Feature-feature relationships
p2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt") +
  ggtitle("UMI Count vs Mitochondrial %")

p3 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  ggtitle("UMI Count vs Gene Count")

p_scatter <- p2 + p3
ggsave("figures/02_qc_scatter_plots.png", p_scatter, width = 12, height = 5, dpi = 300)

# Apply quality filters
# nFeature_RNA > 200: Remove cells with very few genes (likely empty droplets)
# nFeature_RNA < 2500: Remove cells with too many genes (likely doublets)
# percent.mt < 5: Remove cells with high mitochondrial content (stressed/dying)
pbmc <- subset(pbmc, 
               subset = nFeature_RNA > 200 & 
                        nFeature_RNA < 2500 & 
                        percent.mt < 5)

cat("\nAfter filtering:\n")
cat("Cells retained:", ncol(pbmc), "\n")
cat("Mean features per cell:", round(mean(pbmc$nFeature_RNA), 1), "\n")
cat("Median features per cell:", median(pbmc$nFeature_RNA), "\n")

# Visualize QC metrics AFTER filtering
p4 <- VlnPlot(pbmc, 
              features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
              ncol = 3, 
              pt.size = 0.1) +
  plot_annotation(title = "QC Metrics After Filtering")

ggsave("figures/03_qc_after_filtering.png", p4, width = 14, height = 5, dpi = 300)

################################################################################
# 3. NORMALIZATION
################################################################################

cat("\n=== NORMALIZATION ===\n")

# Normalize data using log-normalization
# Global-scaling normalization: Normalizes feature expression for each cell by
# total expression, multiplies by scale factor (10,000), and log-transforms
pbmc <- NormalizeData(pbmc, 
                      normalization.method = "LogNormalize", 
                      scale.factor = 10000)

# Identify highly variable features (genes)
# These genes show high cell-to-cell variation (highly expressed in some cells, 
# lowly expressed in others). Focusing on these genes highlights biological signal.
pbmc <- FindVariableFeatures(pbmc, 
                              selection.method = "vst", 
                              nfeatures = 2000)

# Identify top 10 most variable genes
top10 <- head(VariableFeatures(pbmc), 10)
cat("\nTop 10 highly variable genes:\n")
print(top10)

# Plot variable features
p5 <- VariableFeaturePlot(pbmc)
p6 <- LabelPoints(plot = p5, points = top10, repel = TRUE)
ggsave("figures/04_variable_features.png", p6, width = 12, height = 7, dpi = 300)

# Scale data
# Linear transformation prior to PCA
# Shifts expression so mean = 0 and variance = 1
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

cat("Normalization complete. Variable features:", length(VariableFeatures(pbmc)), "\n")

################################################################################
# 4. DIMENSIONALITY REDUCTION (PCA)
################################################################################

cat("\n=== DIMENSIONALITY REDUCTION ===\n")

# Run Principal Component Analysis (PCA)
# Uses highly variable genes identified above
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc), verbose = FALSE)

# Visualize PCA results
p7 <- DimPlot(pbmc, reduction = "pca") + 
  ggtitle("PCA of PBMC Dataset")
ggsave("figures/05_pca.png", p7, width = 8, height = 7, dpi = 300)

# Elbow plot to determine dimensionality
# Helps determine how many PCs to use for downstream analysis
p8 <- ElbowPlot(pbmc, ndims = 50) +
  ggtitle("Elbow Plot - PC Selection")
ggsave("figures/06_elbow_plot.png", p8, width = 10, height = 6, dpi = 300)

cat("PCA complete. Using first 10 PCs for downstream analysis.\n")

################################################################################
# 5. CLUSTERING
################################################################################

cat("\n=== CLUSTERING ===\n")

# Construct k-nearest neighbor (KNN) graph
# Uses first 10 principal components
pbmc <- FindNeighbors(pbmc, dims = 1:10)

# Cluster cells using Louvain algorithm
# Resolution parameter sets granularity (higher = more clusters)
# 0.4-1.2 typically returns good results for ~3K cells
pbmc <- FindClusters(pbmc, resolution = 0.5, verbose = FALSE)

# Check cluster sizes
cat("\nCluster distribution:\n")
print(table(Idents(pbmc)))

# Run UMAP for visualization
# Non-linear dimensionality reduction
pbmc <- RunUMAP(pbmc, dims = 1:10, verbose = FALSE)

# Visualize clusters on UMAP
p9 <- DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) +
  ggtitle("UMAP - Cell Clusters") +
  theme(legend.position = "right")
ggsave("figures/07_umap_clusters.png", p9, width = 10, height = 8, dpi = 300)

cat("Clustering complete. Identified", length(unique(Idents(pbmc))), "clusters.\n")

################################################################################
# 6. FIND CLUSTER MARKERS
################################################################################

cat("\n=== IDENTIFYING CLUSTER MARKERS ===\n")

# Find markers for all clusters
# only.pos = TRUE: Return only positive markers
# min.pct = 0.25: Only test genes detected in at least 25% of cells
# logfc.threshold = 0.25: Only test genes with log2FC > 0.25
pbmc.markers <- FindAllMarkers(pbmc, 
                                only.pos = TRUE, 
                                min.pct = 0.25, 
                                logfc.threshold = 0.25,
                                verbose = FALSE)

# Get top 10 markers per cluster
top10_markers <- pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

# Save marker results
write.csv(pbmc.markers, "results/cluster_markers_all.csv", row.names = FALSE)
write.csv(top10_markers, "results/cluster_markers_top10.csv", row.names = FALSE)

# Print top 3 markers per cluster
cat("\nTop 3 marker genes per cluster:\n")
for(i in sort(unique(pbmc.markers$cluster))) {
  cat("\n=== Cluster", i, "===\n")
  top_genes <- pbmc.markers %>% 
    filter(cluster == i) %>% 
    arrange(desc(avg_log2FC)) %>%
    head(3)
  print(top_genes[, c("gene", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")])
}

# Heatmap of top markers
p10 <- DoHeatmap(pbmc, features = top10_markers$gene) + 
  NoLegend() +
  ggtitle("Top 10 Marker Genes per Cluster")
ggsave("figures/08_marker_heatmap.png", p10, width = 14, height = 18, dpi = 300)

################################################################################
# 7. CELL TYPE ANNOTATION
################################################################################

cat("\n=== CELL TYPE ANNOTATION ===\n")

# Annotate clusters based on canonical markers
# Based on expression patterns and known immune cell markers
new.cluster.ids <- c("Naive CD4 T",    # Cluster 0: IL7R+, CCR7+
                     "CD14+ Mono",      # Cluster 1: CD14+, LYZ+
                     "Memory CD4 T",    # Cluster 2: IL7R+, S100A4+
                     "B",               # Cluster 3: MS4A1+
                     "CD8 T",           # Cluster 4: CD8A+
                     "FCGR3A+ Mono",    # Cluster 5: FCGR3A+, MS4A7+
                     "NK",              # Cluster 6: GNLY+, NKG7+
                     "DC",              # Cluster 7: FCER1A+, CST3+
                     "Platelet")        # Cluster 8: PPBP+

names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)

# UMAP with cell type labels
p11 <- DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + 
  NoLegend() +
  ggtitle("PBMC Cell Types")
ggsave("figures/09_umap_celltypes.png", p11, width = 10, height = 8, dpi = 300)

# Canonical marker expression
canonical_markers <- c("IL7R", "CCR7",      # Naive CD4+ T
                       "CD14", "LYZ",        # CD14+ Mono
                       "MS4A1",              # B cells
                       "CD8A",               # CD8+ T
                       "FCGR3A", "MS4A7",    # FCGR3A+ Mono
                       "GNLY", "NKG7",       # NK
                       "FCER1A", "CST3",     # DC
                       "PPBP")               # Platelets

p12 <- FeaturePlot(pbmc, features = canonical_markers[1:9], 
                   ncol = 3, pt.size = 0.1)
ggsave("figures/10_canonical_markers.png", p12, width = 15, height = 15, dpi = 300)

# Dot plot of marker expression
p13 <- DotPlot(pbmc, features = canonical_markers) + 
  RotatedAxis() +
  ggtitle("Canonical Marker Expression Across Cell Types")
ggsave("figures/11_dotplot_markers.png", p13, width = 14, height = 6, dpi = 300)

cat("\nCell type annotation complete!\n")
cat("\nCell type counts:\n")
print(table(Idents(pbmc)))

################################################################################
# 8. DIFFERENTIAL EXPRESSION ANALYSIS
################################################################################

cat("\n=== DIFFERENTIAL EXPRESSION ANALYSIS ===\n")

# Compare monocyte populations: CD14+ vs FCGR3A+
mono_de <- FindMarkers(pbmc, 
                       ident.1 = "CD14+ Mono", 
                       ident.2 = "FCGR3A+ Mono",
                       min.pct = 0.25,
                       verbose = FALSE)

mono_de$gene <- rownames(mono_de)

# Flag significant genes
mono_de$significant <- ifelse(mono_de$p_val_adj < 0.05 & 
                               abs(mono_de$avg_log2FC) > 0.5, 
                               "Significant", "Not Significant")

# Save results
write.csv(mono_de, "results/monocyte_differential_expression.csv", row.names = FALSE)

# Get top differentially expressed genes
top_de <- mono_de %>% 
  filter(significant == "Significant") %>% 
  arrange(desc(abs(avg_log2FC))) %>% 
  head(10)

cat("\nTop 10 differentially expressed genes (CD14+ vs FCGR3A+ Monocytes):\n")
print(top_de[, c("gene", "avg_log2FC", "p_val_adj")])

# Volcano plot
library(ggrepel)

p14 <- ggplot(mono_de, aes(x = avg_log2FC, y = -log10(p_val_adj), 
                            color = significant)) +
  geom_point(alpha = 0.5, size = 1.5) +
  geom_text_repel(data = top_de, aes(label = gene), 
                  size = 3, max.overlaps = 20, box.padding = 0.5) +
  scale_color_manual(values = c("grey60", "red3")) +
  theme_minimal() +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey30") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "grey30") +
  labs(title = "Differential Expression: CD14+ vs FCGR3A+ Monocytes",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-value",
       color = "") +
  theme(legend.position = "top",
        plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("figures/12_volcano_plot.png", p14, width = 10, height = 8, dpi = 300)

# Violin plots for top DE genes
top5_genes <- head(top_de$gene, 5)
p15 <- VlnPlot(pbmc, features = top5_genes, 
               idents = c("CD14+ Mono", "FCGR3A+ Mono"),
               ncol = 5, pt.size = 0) +
  plot_annotation(title = "Top 5 Differentially Expressed Genes")

ggsave("figures/13_de_violin_plots.png", p15, width = 18, height = 4, dpi = 300)

cat("\nDifferential expression complete!\n")
cat("Significant genes (|log2FC| > 0.5, padj < 0.05):", 
    sum(mono_de$significant == "Significant"), "\n")

################################################################################
# 9. CELL TYPE PROPORTIONS
################################################################################

cat("\n=== CELL TYPE PROPORTIONS ===\n")

# Calculate cell type frequencies
cell_counts <- table(Idents(pbmc))
cell_props <- prop.table(cell_counts) * 100

prop_df <- data.frame(
  CellType = names(cell_props),
  Proportion = as.numeric(cell_props),
  Count = as.numeric(cell_counts)
)

# Save proportions
write.csv(prop_df, "results/cell_type_proportions.csv", row.names = FALSE)

# Bar plot
p16 <- ggplot(prop_df, aes(x = reorder(CellType, -Proportion), 
                            y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(aes(label = sprintf("%.1f%%", Proportion)), 
            vjust = -0.5, size = 3.5) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(title = "Cell Type Composition - PBMC Dataset",
       x = "Cell Type",
       y = "Percentage (%)")

ggsave("figures/14_cell_proportions.png", p16, width = 12, height = 7, dpi = 300)

cat("\nCell type proportions:\n")
print(prop_df)

################################################################################
# 10. SAVE FINAL OBJECT
################################################################################

cat("\n=== SAVING RESULTS ===\n")

# Save final Seurat object
saveRDS(pbmc, "results/pbmc_final_annotated.rds")

cat("\n✓ Analysis complete!\n")
cat("Results saved in: results/\n")
cat("Figures saved in: figures/\n")
cat("\nTotal cells analyzed:", ncol(pbmc), "\n")
cat("Cell types identified:", length(unique(Idents(pbmc))), "\n")

################################################################################
# SESSION INFO
################################################################################

cat("\n=== SESSION INFORMATION ===\n")
sessionInfo()

################################################################################
# END OF SCRIPT
################################################################################
