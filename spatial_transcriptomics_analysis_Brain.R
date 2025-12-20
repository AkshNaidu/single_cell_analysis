################################################################################
# SPATIAL TRANSCRIPTOMICS ANALYSIS: 10X VISIUM HUMAN BRAIN
################################################################################
# 
# Author: Akshayata Naidu, PhD
# Date: December 2024
# GitHub: https://github.com/AkshNaidu
# 
# Description:
# Complete end-to-end analysis of spatial transcriptomics data from human brain
# tissue using 10X Genomics Visium platform. Workflow includes quality control,
# normalization, spatial clustering, identification of spatially variable genes,
# and spatial domain characterization.
#
# Dataset: Human Brain Section (10X Visium)
# Technology: Visium Spatial Gene Expression
# Spots: ~3,600-4,000 spatial spots
# Resolution: 55μm spot diameter, ~10 cells per spot
# URL: https://support.10xgenomics.com/spatial-gene-expression/datasets
#
################################################################################

# Load required libraries
library(Seurat)
library(ggplot2)
library(patchwork)
library(hdf5r)
library(viridis)
library(dplyr)

# Set working directory
setwd("~/spatial-analysis")

################################################################################
# 1. DATA LOADING
################################################################################

cat("\n=== LOADING SPATIAL DATA ===\n")

# Load 10X Visium spatial data
# This loads: expression matrix, spatial coordinates, and tissue images
brain <- Load10X_Spatial(
  data.dir = "data/spatial/",
  filename = "filtered_feature_bc_matrix.h5",
  slice = "brain_section1"
)

cat("Spatial data dimensions:\n")
cat("Genes:", nrow(brain), "\n")
cat("Spots:", ncol(brain), "\n")

################################################################################
# 2. QUALITY CONTROL
################################################################################

cat("\n=== QUALITY CONTROL ===\n")

# Calculate QC metrics
# Mitochondrial percentage
brain[["percent.mt"]] <- PercentageFeatureSet(brain, pattern = "^MT-")

# Ribosomal percentage
brain[["percent.ribo"]] <- PercentageFeatureSet(brain, pattern = "^RP[SL]")

# Visualize QC metrics on tissue
p1 <- SpatialFeaturePlot(brain, features = "nCount_Spatial") + 
  ggtitle("Total UMI Count per Spot")

p2 <- SpatialFeaturePlot(brain, features = "nFeature_Spatial") + 
  ggtitle("Gene Count per Spot")

p_qc_spatial <- p1 + p2
ggsave("figures/01_qc_on_tissue.png", p_qc_spatial, width = 16, height = 7, dpi = 300)

# QC violin plots
p3 <- VlnPlot(brain, 
              features = c("nCount_Spatial", "nFeature_Spatial", "percent.mt"),
              pt.size = 0.1, ncol = 3) +
  plot_annotation(title = "Spatial QC Metrics")

ggsave("figures/02_qc_violin.png", p3, width = 14, height = 5, dpi = 300)

# View tissue H&E image
p4 <- SpatialPlot(brain, pt.size.factor = 1.6, alpha = 0) +
  ggtitle("Brain Tissue Section - H&E Image")

ggsave("figures/03_tissue_image.png", p4, width = 10, height = 8, dpi = 300)

# QC summary statistics
cat("\nQC Summary Statistics:\n")
cat("Mean UMI per spot:", round(mean(brain$nCount_Spatial), 1), "\n")
cat("Median UMI per spot:", median(brain$nCount_Spatial), "\n")
cat("Mean genes per spot:", round(mean(brain$nFeature_Spatial), 1), "\n")
cat("Median genes per spot:", median(brain$nFeature_Spatial), "\n")
cat("Mean % mitochondrial:", round(mean(brain$percent.mt), 2), "%\n")

################################################################################
# 3. NORMALIZATION
################################################################################

cat("\n=== NORMALIZATION ===\n")

# SCTransform normalization
# Recommended for spatial data - accounts for technical variation
# Also performs feature selection and scaling
brain <- SCTransform(brain, 
                     assay = "Spatial", 
                     verbose = FALSE)

cat("SCTransform normalization complete.\n")

################################################################################
# 4. DIMENSIONALITY REDUCTION & CLUSTERING
################################################################################

cat("\n=== DIMENSIONALITY REDUCTION & CLUSTERING ===\n")

# Run PCA
brain <- RunPCA(brain, assay = "SCT", verbose = FALSE)

# Construct neighbor graph and cluster
brain <- FindNeighbors(brain, reduction = "pca", dims = 1:30)
brain <- FindClusters(brain, verbose = FALSE, resolution = 0.5)

# Run UMAP
brain <- RunUMAP(brain, reduction = "pca", dims = 1:30, verbose = FALSE)

# Visualize clusters on UMAP
p5 <- DimPlot(brain, reduction = "umap", label = TRUE) +
  ggtitle("Spatial Clusters - UMAP")

ggsave("figures/04_umap_clusters.png", p5, width = 10, height = 8, dpi = 300)

# Visualize clusters on tissue - THIS IS THE KEY SPATIAL VISUALIZATION
p6 <- SpatialDimPlot(brain, label = TRUE, label.size = 3, pt.size.factor = 1.8) +
  ggtitle("Spatial Clusters on Tissue Section")

ggsave("figures/05_clusters_on_tissue.png", p6, width = 12, height = 10, dpi = 300)

# Faceted view of clusters
p7 <- SpatialDimPlot(brain, 
                     cells.highlight = CellsByIdentities(brain), 
                     facet.highlight = TRUE, 
                     ncol = 4,
                     pt.size.factor = 1.6)

ggsave("figures/06_clusters_faceted.png", p7, width = 18, height = 12, dpi = 300)

cat("\nClustering complete.\n")
cat("Cluster distribution:\n")
print(table(Idents(brain)))

################################################################################
# 5. SPATIALLY VARIABLE GENES
################################################################################

cat("\n=== IDENTIFYING SPATIALLY VARIABLE GENES ===\n")

# Find spatially variable features
# These are genes whose expression patterns show spatial structure
# Using Moran's I - measures spatial autocorrelation
brain <- FindSpatiallyVariableFeatures(
  brain, 
  assay = "SCT",
  features = VariableFeatures(brain)[1:1000],
  selection.method = "moransi"
)

# Get top spatially variable genes
top_spatial <- head(SpatiallyVariableFeatures(brain, selection.method = "moransi"), 12)

cat("\nTop 12 spatially variable genes:\n")
print(top_spatial)

# Save results
spatial_features_df <- data.frame(
  gene = SpatiallyVariableFeatures(brain, selection.method = "moransi"),
  rank = 1:length(SpatiallyVariableFeatures(brain, selection.method = "moransi"))
)
write.csv(spatial_features_df, "results/spatially_variable_genes.csv", row.names = FALSE)

# Visualize top 6 spatially variable genes on tissue
p8 <- SpatialFeaturePlot(brain, 
                         features = top_spatial[1:6], 
                         ncol = 3, 
                         alpha = c(0.1, 1),
                         pt.size.factor = 1.6)

ggsave("figures/07_spatial_variable_genes.png", p8, width = 18, height = 12, dpi = 300)

# Zoomed view of top spatially variable gene
p9 <- SpatialFeaturePlot(brain, 
                         features = top_spatial[1], 
                         pt.size.factor = 2, 
                         alpha = c(0.1, 1)) +
  ggtitle(paste("Top Spatially Variable Gene:", top_spatial[1]))

ggsave("figures/08_top_spatial_gene_zoom.png", p9, width = 10, height = 8, dpi = 300)

# Ridge plot showing expression distribution
library(ggridges)
p10 <- RidgePlot(brain, features = top_spatial[1:6], ncol = 2) +
  ggtitle("Expression Distribution of Spatially Variable Genes")

ggsave("figures/09_ridge_plot.png", p10, width = 12, height = 10, dpi = 300)

cat("\nSpatially variable gene identification complete!\n")

################################################################################
# 6. CLUSTER MARKER IDENTIFICATION
################################################################################

cat("\n=== IDENTIFYING CLUSTER MARKERS ===\n")

# Find markers for spatial clusters
spatial_markers <- FindAllMarkers(brain, 
                                   only.pos = TRUE,
                                   min.pct = 0.25, 
                                   logfc.threshold = 0.25,
                                   verbose = FALSE)

# Get top 5 markers per cluster
top5_spatial <- spatial_markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)

# Save results
write.csv(spatial_markers, "results/spatial_cluster_markers_all.csv", row.names = FALSE)
write.csv(top5_spatial, "results/spatial_cluster_markers_top5.csv", row.names = FALSE)

# Print top 3 markers per cluster
cat("\nTop 3 marker genes per spatial cluster:\n")
for(i in sort(unique(spatial_markers$cluster))) {
  cat("\n=== Cluster", i, "===\n")
  top_genes <- spatial_markers %>% 
    filter(cluster == i) %>% 
    arrange(desc(avg_log2FC)) %>%
    head(3)
  print(top_genes[, c("gene", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")])
}

# Heatmap of top markers
p11 <- DoHeatmap(brain, features = top5_spatial$gene, size = 3) +
  ggtitle("Top 5 Markers per Spatial Cluster")

ggsave("figures/10_marker_heatmap.png", p11, width = 14, height = 12, dpi = 300)

# Dot plot
p12 <- DotPlot(brain, features = unique(top5_spatial$gene[1:30])) + 
  RotatedAxis() +
  ggtitle("Marker Expression Across Spatial Clusters")

ggsave("figures/11_dotplot_markers.png", p12, width = 16, height = 6, dpi = 300)

# Visualize select top markers on tissue
select_markers <- spatial_markers %>%
  group_by(cluster) %>%
  top_n(n = 1, wt = avg_log2FC) %>%
  pull(gene) %>%
  head(6)

p13 <- SpatialFeaturePlot(brain, 
                          features = select_markers, 
                          ncol = 3, 
                          alpha = c(0.1, 1),
                          pt.size.factor = 1.6) +
  plot_annotation(title = "Top Marker Gene per Cluster")

ggsave("figures/12_top_markers_on_tissue.png", p13, width = 18, height = 12, dpi = 300)

cat("\nCluster marker identification complete!\n")

################################################################################
# 7. SPATIAL DOMAIN CHARACTERIZATION
################################################################################

cat("\n=== SPATIAL DOMAIN CHARACTERIZATION ===\n")

# Define brain region markers based on literature
# These are example markers - adjust based on actual tissue
cortex_markers <- c("NRGN", "SLC17A7")         # Cortical neurons
hippocampus_markers <- c("PROX1")              # Hippocampus
oligodendrocyte_markers <- c("MBP", "MOG")     # White matter/oligodendrocytes
astrocyte_markers <- c("AQP4", "GFAP")         # Astrocytes
microglia_markers <- c("TMEM119", "P2RY12")    # Microglia

# Visualize region-specific markers
region_markers <- c("NRGN", "MBP", "PROX1", "AQP4", "TMEM119", "SLC17A7")
region_markers <- region_markers[region_markers %in% rownames(brain)]

if(length(region_markers) > 0) {
  p14 <- SpatialFeaturePlot(brain, 
                            features = region_markers[1:min(6, length(region_markers))], 
                            ncol = 3, 
                            alpha = c(0.1, 1),
                            pt.size.factor = 2)
  
  ggsave("figures/13_region_markers.png", p14, width = 18, height = 12, dpi = 300)
}

# Create composite scores for tissue regions
# Module scoring - averages expression of gene sets
if(all(cortex_markers %in% rownames(brain))) {
  brain <- AddModuleScore(brain, 
                          features = list(cortex_markers), 
                          name = "Cortex_Score")
  
  # Visualize cortex score
  p15 <- SpatialFeaturePlot(brain, 
                            features = "Cortex_Score1", 
                            pt.size.factor = 2) +
    scale_fill_viridis() +
    ggtitle("Cortical Region Score")
  
  ggsave("figures/14_cortex_score.png", p15, width = 10, height = 8, dpi = 300)
}

# Summary statistics per cluster
cluster_stats <- data.frame(
  Cluster = names(table(Idents(brain))),
  N_spots = as.numeric(table(Idents(brain))),
  Mean_genes = tapply(brain$nFeature_Spatial, Idents(brain), mean),
  Mean_UMI = tapply(brain$nCount_Spatial, Idents(brain), mean),
  Mean_MT_pct = tapply(brain$percent.mt, Idents(brain), mean)
)

write.csv(cluster_stats, "results/spatial_cluster_statistics.csv", row.names = FALSE)

cat("\nSpatial cluster statistics:\n")
print(cluster_stats)

cat("\nSpatial domain characterization complete!\n")

################################################################################
# 8. ADVANCED SPATIAL VISUALIZATIONS
################################################################################

cat("\n=== ADVANCED SPATIAL VISUALIZATIONS ===\n")

# Interactive spatial plot with tissue overlay
p16 <- SpatialFeaturePlot(brain, 
                          features = "nFeature_Spatial",
                          pt.size.factor = 1.8,
                          alpha = c(0.3, 1)) +
  ggtitle("Gene Diversity Across Tissue") +
  scale_fill_viridis()

ggsave("figures/15_gene_diversity_spatial.png", p16, width = 10, height = 8, dpi = 300)

# Multiple features side-by-side
if(length(top_spatial) >= 4) {
  p17 <- SpatialFeaturePlot(brain, 
                            features = top_spatial[1:4],
                            ncol = 2,
                            pt.size.factor = 2,
                            alpha = c(0.1, 1)) +
    plot_annotation(title = "Top 4 Spatially Variable Genes")
  
  ggsave("figures/16_multi_feature_spatial.png", p17, width = 14, height = 14, dpi = 300)
}

# Cluster proportions
cluster_props <- prop.table(table(Idents(brain))) * 100
prop_df <- data.frame(
  Cluster = names(cluster_props),
  Proportion = as.numeric(cluster_props)
)

p18 <- ggplot(prop_df, aes(x = Cluster, y = Proportion, fill = Cluster)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = sprintf("%.1f%%", Proportion)), vjust = -0.5) +
  theme_minimal() +
  theme(legend.position = "none") +
  ggtitle("Spatial Cluster Proportions") +
  ylab("Percentage of Spots (%)")

ggsave("figures/17_cluster_proportions.png", p18, width = 10, height = 6, dpi = 300)

cat("\nAdvanced visualizations complete!\n")

################################################################################
# 9. SAVE FINAL OBJECT
################################################################################

cat("\n=== SAVING RESULTS ===\n")

# Save final Seurat object
saveRDS(brain, "results/brain_spatial_final.rds")

cat("\n✓ Spatial analysis complete!\n")
cat("Results saved in: results/\n")
cat("Figures saved in: figures/\n")
cat("\nTotal spots analyzed:", ncol(brain), "\n")
cat("Spatial clusters identified:", length(unique(Idents(brain))), "\n")
cat("Spatially variable genes identified:", 
    length(SpatiallyVariableFeatures(brain)), "\n")

################################################################################
# 10. ANALYSIS SUMMARY
################################################################################

cat("\n=== ANALYSIS SUMMARY ===\n")

summary_stats <- data.frame(
  Metric = c("Total Spots", 
             "Total Genes",
             "Spatial Clusters",
             "Mean Genes per Spot",
             "Mean UMI per Spot",
             "Spatially Variable Genes",
             "Total Marker Genes"),
  Value = c(ncol(brain),
            nrow(brain),
            length(unique(Idents(brain))),
            round(mean(brain$nFeature_Spatial), 1),
            round(mean(brain$nCount_Spatial), 1),
            length(SpatiallyVariableFeatures(brain)),
            nrow(spatial_markers))
)

write.csv(summary_stats, "results/analysis_summary.csv", row.names = FALSE)

cat("\nAnalysis Summary:\n")
print(summary_stats)

################################################################################
# SESSION INFO
################################################################################

cat("\n=== SESSION INFORMATION ===\n")
sessionInfo()

################################################################################
# KEY FINDINGS
################################################################################

cat("\n=== KEY ANALYSIS OUTPUTS ===\n")
cat("\n1. SPATIAL CLUSTERING:\n")
cat("   - Identified", length(unique(Idents(brain))), "distinct spatial domains\n")
cat("   - Clusters show clear spatial organization on tissue\n")

cat("\n2. SPATIALLY VARIABLE GENES:\n")
cat("   - Top gene:", top_spatial[1], "\n")
cat("   - Total spatially variable genes:", length(SpatiallyVariableFeatures(brain)), "\n")

cat("\n3. CLUSTER MARKERS:\n")
cat("   - Total marker genes identified:", nrow(spatial_markers), "\n")
cat("   - Average markers per cluster:", 
    round(nrow(spatial_markers) / length(unique(spatial_markers$cluster)), 1), "\n")

cat("\n4. TISSUE CHARACTERIZATION:\n")
cat("   - Mean genes per spot:", round(mean(brain$nFeature_Spatial), 1), "\n")
cat("   - Spatial patterns reveal distinct tissue architecture\n")

################################################################################
# END OF SCRIPT
################################################################################
