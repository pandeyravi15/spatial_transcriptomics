# 🧠 Spatial Transcriptomics Analysis: 10x Genomics Visium Human Brain

[![Python](https://img.shields.io/badge/Python-3.10%2B-blue?logo=python)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Analysis](https://img.shields.io/badge/Analysis-Spatial%20Transcriptomics-purple)]()
[![Platform](https://img.shields.io/badge/Platform-10x%20Visium-orange)]()

> **End-to-end spatial transcriptomics analysis pipeline for 10x Genomics Visium data, demonstrating tissue architecture deconvolution, spatially-resolved marker gene expression, and Moran's I spatial autocorrelation in human brain FFPE tissue.**

---

## 📋 Table of Contents

- [Overview](#overview)
- [Dataset](#dataset)
- [Pipeline Overview](#pipeline-overview)
- [Key Results](#key-results)
- [Repository Structure](#repository-structure)
- [Installation & Usage](#installation--usage)
- [Methods](#methods)
- [Figures](#figures)
- [Future Directions](#future-directions)
- [Author](#author)

---

## Overview

This repository demonstrates a complete, reproducible spatial transcriptomics analysis workflow applied to 10x Genomics Visium human brain FFPE data. The pipeline covers all critical steps from raw count matrix ingestion through quality control, normalization, dimensionality reduction, unsupervised clustering, spatially-resolved marker gene visualization, differential expression analysis, and spatial autocorrelation (Moran's I).

The analysis identifies **8 distinct spatial domains** within human brain tissue — including tumor core, neuronal layers, oligodendrocytes, astrocytes, microglia, and vascular endothelium — and characterizes their molecular signatures using established marker genes.

**Key technical skills demonstrated:**
- Production-grade NGS analysis pipeline design
- Single-cell/spatial RNA-seq preprocessing & QC
- Dimensionality reduction (PCA, UMAP/t-SNE)
- Graph-based clustering (KNN + Leiden-like)
- Differential expression with multiple testing correction (Benjamini-Hochberg)
- Spatial autocorrelation analysis (Moran's I)
- Biologically interpretable visualization

---

## Dataset

| Property | Value |
|----------|-------|
| **Source** | 10x Genomics Visium HD — Human Brain FFPE |
| **Technology** | Visium Spatial Gene Expression |
| **Tissue** | Human Brain (FFPE) |
| **Download** | [10x Genomics Dataset Portal](https://www.10xgenomics.com/datasets) |
| **Spots** | ~2,400 tissue-covered spots |
| **Genes** | ~5,000 (post-filtering) |

> **Note:** This analysis uses a biologically realistic simulation of Visium data for reproducibility and portability. To run on actual SpaceRanger output, replace the data loading step with:
> ```python
> import scanpy as sc
> adata = sc.read_visium("path/to/spaceranger/output/")
> ```

---

## Pipeline Overview

```
Raw Count Matrix
      │
      ▼
┌─────────────────────┐
│  1. Quality Control  │  Filter low-quality spots (UMI, genes, % mito)
└─────────────────────┘
      │
      ▼
┌─────────────────────────────────┐
│  2. Normalization & Log1p        │  Library-size normalize → 10,000 CPM → log1p
└─────────────────────────────────┘
      │
      ▼
┌────────────────────────────────────┐
│  3. HVG Selection (n=2,000)         │  Dispersion-based; retain all marker genes
└────────────────────────────────────┘
      │
      ▼
┌───────────────────────────────────────┐
│  4. Dimensionality Reduction           │  PCA (50 PCs) → UMAP/t-SNE (top 20 PCs)
└───────────────────────────────────────┘
      │
      ▼
┌──────────────────────────────────────────┐
│  5. Graph-Based Clustering                │  KNN (k=15) → Leiden/KMeans (8 clusters)
└──────────────────────────────────────────┘
      │
      ▼
┌──────────────────────────────────────────────┐
│  6. Spatial Visualization                      │  Cluster maps + marker gene expression
└──────────────────────────────────────────────┘
      │
      ▼
┌──────────────────────────────────────────────────┐
│  7. Differential Expression                        │  Mann-Whitney U + BH correction
└──────────────────────────────────────────────────┘
      │
      ▼
┌────────────────────────────────────────────────────────┐
│  8. Spatial Autocorrelation (Moran's I)                 │  Per-gene spatial structure
└────────────────────────────────────────────────────────┘
      │
      ▼
   Results & Figures
```

---

## Key Results

### Spatial Domains Identified

| Cluster | Cell Type / Region | Key Markers | N Spots |
|---------|-------------------|-------------|---------|
| 0 | Tumor Core | MKI67, TOP2A, PCNA, CDK1 | ~216 |
| 1 | Proliferating Glia | GFAP, VIM, NES, SOX2 | ~375 |
| 2 | Neurons (L2/3) | RBFOX3, CUX1, CUX2, LAMP5 | ~282 |
| 3 | Neurons (L5/6) | RORB, TLE4, FEZF2, BCL11B | ~617 |
| 4 | Oligodendrocytes | MBP, MOG, MOBP, PLP1 | ~188 |
| 5 | Astrocytes | GFAP, AQP4, ALDH1L1, SLC1A2 | ~179 |
| 6 | Microglia / Immune | AIF1, CSF1R, TREM2, P2RY12 | ~127 |
| 7 | Vascular / Endothelial | CLDN5, PECAM1, VWF, FLT1 | ~279 |

### Spatial Autocorrelation (Moran's I)

Genes with high Moran's I are spatially restricted to specific tissue regions:

| Gene | Moran's I | Interpretation |
|------|-----------|---------------|
| CLDN5 | 0.71 | Highly spatially restricted (vascular) |
| AIF1 | 0.60 | Microglia spatially clustered |
| MBP | 0.55 | White matter tracts |
| GFAP | 0.54 | Reactive glia zones |
| AQP4 | 0.45 | Astrocyte endfeet |
| RBFOX3 | 0.30 | Neuronal layers |
| MKI67 | 0.23 | Proliferating zone |

---

## Repository Structure

```
spatial_transcriptomics_portfolio/
├── notebooks/
│   └── 01_visium_spatial_analysis.py    # Main analysis pipeline
├── src/
│   ├── qc_utils.py                       # QC helper functions
│   ├── spatial_utils.py                  # Spatial analysis utilities
│   └── plotting.py                       # Custom visualization functions
├── results/
│   ├── figures/
│   │   ├── 01_QC_metrics.png
│   │   ├── 02_HVG_selection.png
│   │   ├── 03_PCA_elbow.png
│   │   ├── 04_UMAP_clusters.png
│   │   ├── 05_spatial_clusters_and_markers.png
│   │   ├── 06_marker_spatial_grid.png
│   │   ├── 07_dotplot_markers.png
│   │   └── 08_morans_I.png
│   ├── tables/
│   │   ├── DE_cluster_markers.csv
│   │   └── morans_I_spatial_autocorrelation.csv
│   └── analysis_summary.json
├── data/
│   └── README.md                         # Data download instructions
├── environment.yml                       # Conda environment specification
├── requirements.txt                      # pip dependencies
└── README.md
```

---

## Installation & Usage

### 1. Clone the repository
```bash
git clone https://github.com/pandeyravi15/spatial-transcriptomics-visium.git
cd spatial-transcriptomics-visium
```

### 2. Set up environment

**Option A — Conda (recommended):**
```bash
conda env create -f environment.yml
conda activate spatial-tx
```

**Option B — pip:**
```bash
pip install -r requirements.txt
```

### 3. Download the data
```bash
# Download Visium HD Human Brain dataset from 10x Genomics
# https://www.10xgenomics.com/datasets/human-brain-cancer-11-mm-capture-area-ffpe-2-standard
# Place SpaceRanger output in data/visium_human_brain/
```

### 4. Run the analysis
```bash
# Full pipeline
python notebooks/01_visium_spatial_analysis.py

# Or interactively in Jupyter
jupyter lab notebooks/01_visium_spatial_analysis.py
```

### 5. View results
```
results/figures/   → All publication-quality figures (PNG, 150 DPI)
results/tables/    → DE markers and Moran's I statistics (CSV)
```

---

## Methods

### Quality Control
Spots were filtered based on three criteria:
- **Total UMI counts:** 500 – 50,000 per spot
- **Genes detected:** ≥ 200 per spot
- **Mitochondrial reads:** ≤ 20% (marker of spot quality/cell damage)

### Normalization
Library-size normalization to 10,000 counts per spot followed by log1p transformation (`log(X + 1)`), consistent with standard scRNA-seq/Visium best practices (Luecken & Theis, 2019).

### HVG Selection
Top 2,000 highly variable genes selected by dispersion (variance/mean), with all known marker genes force-retained to ensure biological interpretability.

### Dimensionality Reduction
PCA on scaled HVG expression matrix (50 components). UMAP/t-SNE computed on top 20 PCs based on elbow plot inspection. In production, UMAP computed via `umap-learn` for superior global structure preservation.

### Clustering
KNN graph (k=15) constructed in PCA space. Community detection via Leiden algorithm (Traag et al., 2019) or KMeans proxy. Resolution parameter tuned to recover biologically meaningful cell-type boundaries.

### Differential Expression
Per-cluster marker genes identified using Mann-Whitney U test (one vs. all), with Benjamini-Hochberg correction at FDR 5%. Log fold-change threshold ≥ 0.25.

### Spatial Autocorrelation
Moran's I computed per gene using KNN-weighted spatial adjacency matrix (k=6 neighbors), quantifying the degree to which gene expression is spatially structured rather than randomly distributed across tissue.

---

## Figures

| Figure | Description |
|--------|-------------|
| `01_QC_metrics.png` | Violin plots of UMI counts, genes per spot, % mitochondrial reads |
| `02_HVG_selection.png` | Mean-variance relationship with HVGs highlighted |
| `03_PCA_elbow.png` | Variance explained per PC; cumulative variance |
| `04_UMAP_clusters.png` | 2D embedding colored by cluster and total UMI |
| `05_spatial_clusters_and_markers.png` | Spatial cluster map + GFAP expression overlay |
| `06_marker_spatial_grid.png` | 8-panel spatial expression grid for key markers |
| `07_dotplot_markers.png` | Dot plot: expression × detection rate per cluster |
| `08_morans_I.png` | Spatial autocorrelation statistics per marker gene |

---

## Future Directions

- [ ] Cell-type deconvolution with RCTD or Cell2location
- [ ] Spatially variable gene analysis (SpatialDE, SPARK-X)
- [ ] Ligand-receptor interaction analysis (CellChat, NicheNet)
- [ ] Integration with paired scRNA-seq reference atlas
- [ ] Trajectory/pseudotime analysis of spatially adjacent clusters
- [ ] Xenium sub-cellular resolution analysis

---

## References

1. Stahl, P.L. et al. (2016). Visualization and analysis of gene expression in tissue sections by spatial transcriptomics. *Science*, 353(6294), 78–82.
2. Luecken, M.D. & Theis, F.J. (2019). Current best practices in single-cell RNA-seq analysis. *Mol Syst Biol*, 15(6), e8746.
3. Traag, V.A. et al. (2019). From Louvain to Leiden: guaranteeing well-connected communities. *Sci Rep*, 9, 5233.
4. Moran, P.A.P. (1950). Notes on continuous stochastic phenomena. *Biometrika*, 37(1/2), 17–23.
5. 10x Genomics Visium Spatial Gene Expression. https://www.10xgenomics.com/products/spatial-gene-expression

---

## Author

**Ravi Shanker Pandey, Ph.D.**
Computational Biologist | The Jackson Laboratory for Genomic Medicine
- 🔗 [LinkedIn](https://www.linkedin.com/in/ravi-pandey-50943980/)
- 🐙 [GitHub](https://github.com/pandeyravi15)
- 📧 ravisp1587@gmail.com
- 🎓 [Google Scholar](https://scholar.google.com/citations?user=FiCtSqQAAAAJ)

---

*This analysis was developed as part of a computational biology portfolio demonstrating spatial transcriptomics expertise. All code is original and reproducible.*
