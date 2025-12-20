# Single-Cell and Spatial Transcriptomics Analysis Portfolio

**Author:** Akshayata Naidu, PhD  
**Contact:** akshyata.naidu@gmail.com  
**LinkedIn:** [linkedin.com/in/naidu-systemsvaccinology](https://www.linkedin.com/in/naidu-systemsvaccinology/)

---

## Overview

This repository contains comprehensive analysis pipelines for both **single-cell RNA-seq** and **spatial transcriptomics** data, demonstrating proficiency in modern transcriptomic analysis methods.

### Key Analyses:
1. **Single-Cell RNA-seq**: PBMC 3k dataset - cell type identification and characterization
2. **Spatial Transcriptomics**: Human brain section - spatial domain mapping and spatially variable genes

---

## Repository Structure

```
.
├── scRNA_seq_analysis_PBMC.R              # Complete single-cell analysis script
├── spatial_transcriptomics_analysis_Brain.R  # Complete spatial analysis script
├── README.md                               # This file
├── data/                                   # Data directory (download separately)
├── figures/                                # Generated figures
└── results/                                # Analysis results (CSV, RDS files)
```

---

## Analysis 1: Single-Cell RNA-seq (PBMC 3k)

### Dataset Information
- **Source:** 10X Genomics
- **Technology:** Chromium Single Cell 3' v1
- **Sample:** Peripheral Blood Mononuclear Cells (PBMCs)
- **Cells:** ~2,700 high-quality cells after QC
- **Genes:** ~13,700 detected

### Analysis Workflow

**Script:** `scRNA_seq_analysis_PBMC.R`

#### Steps:
1. **Data Loading & QC**
   - Quality filtering (200-2,500 genes/cell, <5% mitochondrial)
   - Retained 2,638 high-quality cells

2. **Normalization & Feature Selection**
   - LogNormalize with scale factor 10,000
   - Identified 2,000 highly variable genes

3. **Dimensionality Reduction**
   - PCA on variable features
   - UMAP for visualization

4. **Clustering**
   - Louvain algorithm (resolution 0.5)
   - Identified 9 distinct cell clusters

5. **Cell Type Annotation**
   - Used canonical immune cell markers
   - Identified cell types:
     - Naive CD4+ T cells (IL7R+, CCR7+)
     - Memory CD4+ T cells (IL7R+, S100A4+)
     - CD8+ T cells (CD8A+)
     - B cells (MS4A1+)
     - CD14+ Monocytes (CD14+, LYZ+)
     - FCGR3A+ Monocytes (FCGR3A+, MS4A7+)
     - NK cells (GNLY+, NKG7+)
     - Dendritic cells (FCER1A+, CST3+)
     - Platelets (PPBP+)

6. **Differential Expression**
   - Compared CD14+ vs FCGR3A+ monocyte populations
   - Identified significant marker genes

7. **Cell Type Proportions**
   - Quantified composition of immune cell populations

### Key Outputs
- **14 figures** showing QC, clustering, markers, and differential expression
- **CSV files** with marker genes and cell proportions
- **Annotated Seurat object** for further analysis

### Representative Results

| Cell Type | Percentage | Top Markers |
|-----------|------------|-------------|
| Naive CD4 T | ~30% | IL7R, CCR7 |
| CD14+ Mono | ~25% | CD14, LYZ |
| Memory CD4 T | ~15% | IL7R, S100A4 |
| B cells | ~12% | MS4A1 |
| CD8 T | ~10% | CD8A |

---

## Analysis 2: Spatial Transcriptomics (Human Brain)

### Dataset Information
- **Source:** 10X Genomics
- **Technology:** Visium Spatial Gene Expression
- **Sample:** Human brain tissue section
- **Spots:** ~3,600-4,000 spatial spots
- **Resolution:** 55μm spots (~10 cells/spot)

### Analysis Workflow

**Script:** `spatial_transcriptomics_analysis_Brain.R`

#### Steps:
1. **Data Loading**
   - Loaded expression matrix, spatial coordinates, and H&E tissue image

2. **Quality Control**
   - Visualized QC metrics spatially on tissue
   - Assessed gene/UMI counts per spot

3. **Normalization**
   - SCTransform - accounts for technical variation
   - Recommended for spatial data

4. **Spatial Clustering**
   - PCA (30 components)
   - Graph-based clustering
   - Identified spatially distinct domains

5. **Spatially Variable Genes**
   - **Moran's I statistic** for spatial autocorrelation
   - Identified genes with spatial expression patterns
   - Visualized top genes on tissue sections

6. **Cluster Marker Identification**
   - Found marker genes for each spatial domain
   - Characterized tissue regions

7. **Spatial Domain Characterization**
   - Identified cortical, white matter, and other regions
   - Used brain region-specific markers (NRGN, MBP, PROX1, etc.)

8. **Advanced Visualizations**
   - Expression patterns overlaid on tissue
   - Multi-feature comparisons
   - Spatial organization of cell types

### Key Outputs
- **17 figures** showing spatial patterns, clusters, and gene expression
- **CSV files** with spatially variable genes and cluster markers
- **Spatial Seurat object** with annotations

### Spatially Variable Genes
Top genes showing distinct spatial patterns were identified using Moran's I, revealing:
- Regional gene expression differences
- Tissue architecture organization
- Spatial domains corresponding to anatomical structures

---

## Technical Skills Demonstrated

### Core Competencies
✅ Single-cell RNA-seq data processing and QC  
✅ Spatial transcriptomics analysis  
✅ Dimensionality reduction (PCA, UMAP)  
✅ Unsupervised clustering algorithms  
✅ Cell type annotation using canonical markers  
✅ Differential expression analysis  
✅ Spatially variable gene identification (Moran's I)  
✅ Spatial domain mapping and characterization  
✅ Data visualization and interpretation  
✅ Reproducible workflow development  

### Software & Tools
- **R**: Seurat, dplyr, ggplot2, patchwork
- **Methods**: LogNormalize, SCTransform, Louvain clustering, Moran's I
- **Visualization**: UMAP, spatial plots, heatmaps, volcano plots

---

## Data Download Instructions

### Single-Cell PBMC 3k
```bash
wget https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz
tar -xzf pbmc3k_filtered_gene_bc_matrices.tar.gz
mv filtered_gene_bc_matrices data/pbmc3k/
```

### Spatial Brain Data
```bash
wget https://cf.10xgenomics.com/samples/spatial-exp/1.1.0/V1_Human_Brain_Section_1/V1_Human_Brain_Section_1_filtered_feature_bc_matrix.h5
wget https://cf.10xgenomics.com/samples/spatial-exp/1.1.0/V1_Human_Brain_Section_1/V1_Human_Brain_Section_1_spatial.tar.gz

tar -xzf V1_Human_Brain_Section_1_spatial.tar.gz
mkdir -p data/spatial
mv V1_Human_Brain_Section_1_filtered_feature_bc_matrix.h5 data/spatial/filtered_feature_bc_matrix.h5
mv spatial data/spatial/
```

---

## Running the Analyses

### Prerequisites
```r
# Install required packages
install.packages(c('Seurat', 'dplyr', 'ggplot2', 'patchwork', 'hdf5r', 'viridis', 'ggridges', 'ggrepel'))
```

### Execute Single-Cell Analysis
```r
source("scRNA_seq_analysis_PBMC.R")
```

### Execute Spatial Analysis
```r
source("spatial_transcriptomics_analysis_Brain.R")
```

Both scripts will:
- Create `figures/` directory with all plots
- Create `results/` directory with CSV files and RDS objects
- Print progress and summary statistics

---

## Key Findings

### Single-Cell Analysis
- Successfully identified 9 immune cell populations in PBMCs
- Characterized monocyte heterogeneity (CD14+ vs FCGR3A+)
- Cell composition analysis revealed expected PBMC frequencies

### Spatial Analysis  
- Mapped distinct spatial domains in brain tissue
- Identified spatially variable genes showing region-specific patterns
- Characterized tissue architecture with cortical and subcortical regions

---

## Applications

These analysis pipelines can be adapted for:
- **Immunology**: Immune cell profiling, vaccine response analysis
- **Oncology**: Tumor microenvironment mapping, immune infiltration
- **Neuroscience**: Brain region characterization, disease modeling
- **Developmental Biology**: Tissue organization, cell fate mapping
- **Biomarker Discovery**: Identification of diagnostic/prognostic markers

---

## Methods Summary

### Statistical Approaches
- **Normalization**: LogNormalize (scRNA-seq), SCTransform (spatial)
- **Differential Expression**: Wilcoxon rank-sum test with Bonferroni correction
- **Spatial Analysis**: Moran's I for spatial autocorrelation
- **Clustering**: Graph-based Louvain algorithm

### Quality Control
- Minimum gene cutoffs (200-2,500 genes/cell)
- Mitochondrial percentage filtering (<5% for scRNA, <10% for spatial)
- Feature-count relationship assessment

---

## Computational Environment

**Software:**
- R version 4.x
- Seurat v4/v5
- Analysis runtime: ~30-60 minutes total

**Hardware:**
- These analyses run on standard laptops
- No HPC required for these dataset sizes

---

## Future Directions

Potential extensions of these analyses:
- Integration of scRNA-seq and spatial data
- Trajectory inference and RNA velocity
- Cell-cell interaction analysis
- Deconvolution of spatial spots
- Integration with GWAS/disease data

---

## References

**Datasets:**
- 10X Genomics: https://www.10xgenomics.com/resources/datasets

**Methods:**
- Seurat: Hao et al., Cell 2021
- Spatial Transcriptomics: Ståhl et al., Science 2016
- Moran's I: Svensson et al., Nature Methods 2018

---

## License

This is an educational/portfolio project demonstrating computational biology skills.  
Data are publicly available from 10X Genomics.

---

## Contact

**Akshayata Naidu, PhD**  
Computational Immunologist  
📧 akshyata.naidu@gmail.com  
🔗 [LinkedIn](https://www.linkedin.com/in/naidu-systemsvaccinology/)  
💻 [GitHub](https://github.com/AkshNaidu)

---

**Last Updated:** December 2024
