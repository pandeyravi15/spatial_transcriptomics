# 🧠 Spatial Transcriptomics Analysis: 10x Genomics Visium Human Brain

[![Python](https://img.shields.io/badge/Python-3.10%2B-blue?logo=python)](https://www.python.org/)
[![Scanpy](https://img.shields.io/badge/Scanpy-1.9%2B-green?logo=python)](https://scanpy.readthedocs.io/)
[![Squidpy](https://img.shields.io/badge/Squidpy-1.3%2B-teal)](https://squidpy.readthedocs.io/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Analysis](https://img.shields.io/badge/Analysis-Spatial%20Transcriptomics-purple)]()
[![Platform](https://img.shields.io/badge/Platform-10x%20Visium-orange)]()

> **End-to-end spatial transcriptomics analysis pipeline for 10x Genomics Visium data, demonstrating tissue architecture deconvolution, spatially-resolved marker gene expression, neighborhood enrichment analysis, and Moran's I spatial autocorrelation in human brain glioblastoma FFPE tissue.**

---

## 📋 Table of Contents

- [Overview](#overview)
- [Notebooks](#notebooks)
- [Dataset](#dataset)
- [Pipeline Overview](#pipeline-overview)
- [Key Results](#key-results)
- [Repository Structure](#repository-structure)
- [Installation & Usage](#installation--usage)
- [Methods](#methods)
- [Figures](#figures)
- [Future Directions](#future-directions)
- [References](#references)
- [Author](#author)

---

## Overview

This repository demonstrates a complete, reproducible spatial transcriptomics analysis workflow applied to 10x Genomics Visium human brain glioblastoma FFPE data. The repository contains **two complementary analysis notebooks**:

1. A **fully portable simulation pipeline** (`01_visium_spatial_analysis.py`) that runs end-to-end without any data download — ideal for testing and demonstrating pipeline logic
2. A **real data analysis notebook** (`02_visium_real_data_analysis.ipynb`) applied to the actual 10x Genomics Human Brain Cancer Visium dataset, with rendered outputs and full biological interpretation

Together, these demonstrate the full spectrum of spatial transcriptomics expertise: from pipeline engineering and methodology to hands-on biological discovery on real tissue data.

**Key technical skills demonstrated:**
- Production-grade NGS analysis pipeline design
- Single-cell/spatial RNA-seq preprocessing, QC, and normalization
- Highly variable gene selection and dimensionality reduction (PCA + UMAP)
- Graph-based clustering with Leiden algorithm
- Spatially-resolved marker gene visualization on H&E tissue
- Differential expression with multiple testing correction (Wilcoxon + Benjamini-Hochberg)
- Spatial neighborhood enrichment analysis (Squidpy)
- Spatial autocorrelation (Moran's I) — identifying spatially structured gene programs
- Reproducible research practices (seed setting, layer preservation, versioned outputs)

---

## Notebooks

| Notebook | Description | Data Required | Outputs |
|----------|-------------|---------------|---------|
| `notebooks/01_visium_spatial_analysis.py` | Full 9-step pipeline on biologically realistic simulated Visium data. Runs immediately, no download needed. | None | 8 figures + 2 tables in `results/simulated/` |
| `notebooks/02_visium_real_data_analysis.ipynb` | End-to-end analysis on real 10x Genomics Human Brain Glioblastoma FFPE dataset. Rendered outputs included — view directly on GitHub without running. | [10x Visium dataset](https://www.10xgenomics.com/datasets/human-brain-cancer-11-mm-capture-area-ffpe-2-standard) | 9 figures + 2 tables in `results/real_data/` |

> **Tip:** The notebook (`02`) includes rendered figures and cell outputs — you can view the full analysis results directly on GitHub without running any code.

---

## Dataset

| Property | Value |
|----------|-------|
| **Source** | 10x Genomics Visium CytAssist — Human Brain Cancer (FFPE) |
| **Technology** | Visium Spatial Gene Expression |
| **Tissue** | Human Brain — Glioblastoma multiforme |
| **Download** | [10x Genomics Dataset Portal](https://www.10xgenomics.com/datasets/human-brain-cancer-11-mm-capture-area-ffpe-2-standard) |
| **Spots under tissue** | 10,878 |
| **Median UMI per spot** | 8,339 |
| **Median genes per spot** | 4,600 |
| **Sequencing depth** | ~805M reads (Illumina NovaSeq 6000) |
| **License** | CC BY 4.0 |

---

## Pipeline Overview

```
Raw SpaceRanger Output (.h5 + spatial/)
              │
              ▼
┌─────────────────────────────────────────┐
│  1. Data Loading                         │  sc.read_visium() — counts + coordinates + H&E
└─────────────────────────────────────────┘
              │
              ▼
┌─────────────────────────────────────────┐
│  2. Quality Control                      │  Filter: UMI ≥500, genes ≥200, MT% ≤20
└─────────────────────────────────────────┘
              │
              ▼
┌─────────────────────────────────────────┐
│  3. Normalization & Log1p               │  Library-size normalize → 10,000 CPM → log(X+1)
└─────────────────────────────────────────┘
              │
              ▼
┌─────────────────────────────────────────┐
│  4. HVG Selection (n=3,000)             │  Seurat dispersion method
└─────────────────────────────────────────┘
              │
              ▼
┌─────────────────────────────────────────┐
│  5. Dimensionality Reduction            │  PCA (50 PCs) → KNN graph → UMAP (top 20 PCs)
└─────────────────────────────────────────┘
              │
              ▼
┌─────────────────────────────────────────┐
│  6. Leiden Clustering                   │  Community detection on KNN graph (res=0.5)
└─────────────────────────────────────────┘
              │
              ▼
┌─────────────────────────────────────────┐
│  7. Spatial Visualization               │  Cluster maps + marker expression on H&E tissue
└─────────────────────────────────────────┘
              │
              ▼
┌─────────────────────────────────────────┐
│  8. Differential Expression             │  Wilcoxon rank-sum + BH correction per cluster
└─────────────────────────────────────────┘
              │
              ▼
┌─────────────────────────────────────────┐
│  9. Spatial Neighborhood Analysis       │  Squidpy: neighborhood enrichment + Moran's I
└─────────────────────────────────────────┘
              │
              ▼
       Results & Figures + Saved AnnData (.h5ad)
```

---

## Key Results

### Spatial Domains Identified (Simulated pipeline)

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

### Spatial Autocorrelation — Moran's I (Simulated data)

| Gene | Moran's I | Interpretation |
|------|-----------|---------------|
| CLDN5 | 0.71 | Highly spatially restricted (vascular) |
| AIF1 | 0.60 | Microglia spatially clustered |
| MBP | 0.55 | White matter tracts |
| GFAP | 0.54 | Reactive glia zones |
| AQP4 | 0.45 | Astrocyte endfeet |
| RBFOX3 | 0.30 | Neuronal layers |
| MKI67 | 0.23 | Proliferating zone |

> **Real data results:** See rendered outputs in `02_visium_real_data_analysis.ipynb`

---

## Repository Structure

```
spatial-transcriptomics-visium/
├── notebooks/
│   ├── 01_visium_spatial_analysis.py               # Portable simulation pipeline (no data needed)
│   └── 02_visium_real_data_analysis.ipynb          # Real data — rendered outputs included
├── src/
│   ├── qc_utils.py                                  # QC helper functions
│   ├── spatial_utils.py                             # Spatial analysis utilities
│   └── plotting.py                                  # Custom visualization functions
├── results/
│   ├── simulated/                                   # Outputs from simulation pipeline
│   │   ├── figures/
│   │   │   ├── 01_QC_metrics.png
│   │   │   ├── 02_HVG_selection.png
│   │   │   ├── 03_PCA_elbow.png
│   │   │   ├── 04_UMAP_clusters.png
│   │   │   ├── 05_spatial_clusters_and_markers.png
│   │   │   ├── 06_marker_spatial_grid.png
│   │   │   ├── 07_dotplot_markers.png
│   │   │   └── 08_morans_I.png
│   │   ├── tables/
│   │   │   ├── DE_cluster_markers.csv
│   │   │   └── morans_I_spatial_autocorrelation.csv
│   │   └── analysis_summary.json
│   └── real_data/                                   # Outputs from real Visium dataset
│       ├── figures/
│       │   ├── 01_QC_violin.png
│       │   ├── 02_HVG_selection.png
│       │   ├── 03_PCA_elbow.png
│       │   ├── 04_UMAP_clusters.png
│       │   ├── 05_spatial_clusters.png              ← clusters on H&E tissue
│       │   ├── 06_spatial_markers.png               ← marker expression on tissue
│       │   ├── 07_marker_dotplot.png
│       │   ├── 08_marker_heatmap.png
│       │   └── 09_neighborhood_enrichment.png       ← Squidpy spatial analysis
│       ├── tables/
│       │   ├── DE_marker_genes_per_cluster.csv
│       │   └── morans_I_spatial_autocorrelation.csv
│       └── visium_brain_analyzed.h5ad               ← saved AnnData (gitignored)
├── data/
│   └── README.md                                    # Data download instructions
├── environment.yml                                  # Conda environment specification
├── requirements.txt                                 # pip dependencies
├── .gitignore
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

**Option C — Install from within Jupyter Lab:**
```python
import sys
!{sys.executable} -m pip install scanpy[leiden] squidpy umap-learn anndata
```
> Restart the kernel after installation before importing packages.

### 3. Run Notebook 1 — Simulated pipeline (no data download needed)
```bash
conda activate spatial-tx
cd notebooks
python 01_visium_spatial_analysis.py
# Outputs saved to results/simulated/
```

### 4. Run Notebook 2 — Real Visium data

**Download the dataset first:**
```
1. Go to: https://www.10xgenomics.com/datasets/human-brain-cancer-11-mm-capture-area-ffpe-2-standard
2. Download: filtered_feature_bc_matrix.h5
3. Download: spatial/ folder (tissue_positions.csv, scalefactors_json.json, tissue images)
4. Place all files in: data/visium_human_brain/
```

**Run in Jupyter Lab:**
```bash
jupyter lab
# Open: notebooks/02_visium_real_data_analysis.ipynb
# Run all cells top to bottom
```

### 5. View results
```
results/simulated/figures/   → Simulated pipeline outputs (8 figures)
results/real_data/figures/   → Real data outputs (9 figures including H&E overlays)
results/real_data/tables/    → DE markers + Moran's I statistics (CSV)
```

---

## Methods

### Quality Control
Spots filtered on three standard criteria: total UMI ≥ 500, genes detected ≥ 200, and mitochondrial reads ≤ 20%. Raw counts preserved in `adata.layers['counts']` before any transformation.

### Normalization
Library-size normalization to 10,000 counts per spot followed by log1p transformation `log(X + 1)`, consistent with Visium best practices (Luecken & Theis, 2019).

### HVG Selection
Top 3,000 highly variable genes using the Seurat dispersion method. HVG subset used exclusively for PCA to reduce noise from lowly expressed genes.

### Dimensionality Reduction
PCA on scaled HVG matrix (50 components). UMAP computed on top 20 PCs selected by elbow plot, via `umap-learn` for superior global structure preservation.

### Clustering
KNN graph (k=15) in PCA space. **Leiden algorithm** (Traag et al., 2019) at resolution=0.5 for community detection. Resolution tunable for coarser or finer boundaries.

### Differential Expression
**Wilcoxon rank-sum test** (one cluster vs. all others) with Benjamini-Hochberg FDR correction — recommended for non-normally distributed count data.

### Spatial Neighborhood Analysis (Squidpy)
- **Neighborhood enrichment:** Permutation-based test identifying cluster pairs significantly co-localized in tissue space
- **Moran's I:** Per-gene spatial autocorrelation using KNN-weighted spatial adjacency (k=6 hexagonal neighbors), quantifying spatially structured expression patterns

---

## Figures

### Simulated Pipeline (`results/simulated/figures/`)

| Figure | Description |
|--------|-------------|
| `01_QC_metrics.png` | Violin plots — UMI counts, genes per spot, % mitochondrial reads |
| `02_HVG_selection.png` | Mean-variance relationship with HVGs highlighted |
| `03_PCA_elbow.png` | Variance explained per PC + cumulative variance |
| `04_UMAP_clusters.png` | 2D embedding colored by cluster and total UMI |
| `05_spatial_clusters_and_markers.png` | Spatial cluster map + GFAP expression overlay |
| `06_marker_spatial_grid.png` | 8-panel spatial expression grid for key cell-type markers |
| `07_dotplot_markers.png` | Dot plot: mean expression × % spots expressing per cluster |
| `08_morans_I.png` | Spatial autocorrelation (Moran's I) bar chart per marker gene |

### Real Data (`results/real_data/figures/`)

| Figure | Description |
|--------|-------------|
| `01_QC_violin.png` | Per-spot QC distributions on real tissue data |
| `02_HVG_selection.png` | HVG selection on real gene expression data |
| `03_PCA_elbow.png` | PCA variance explained — real Visium data |
| `04_UMAP_clusters.png` | UMAP embedding — Leiden clusters + QC metrics |
| `05_spatial_clusters.png` | **Leiden clusters overlaid on H&E tissue image** |
| `06_spatial_markers.png` | **MKI67, GFAP, MBP, AIF1 expression on tissue** |
| `07_marker_dotplot.png` | Top 5 marker genes per cluster — dot plot |
| `08_marker_heatmap.png` | Marker gene heatmap across clusters |
| `09_neighborhood_enrichment.png` | **Squidpy spatial neighborhood enrichment matrix** |

---

## Future Directions

- [ ] Cell-type deconvolution with RCTD or Cell2location
- [ ] Spatially variable gene analysis (SpatialDE, SPARK-X)
- [ ] Ligand-receptor interaction analysis (CellChat, NicheNet)
- [ ] Integration with paired scRNA-seq reference atlas (Human Cell Atlas)
- [ ] Trajectory / pseudotime analysis of spatially adjacent clusters
- [ ] Xenium sub-cellular resolution analysis
- [ ] Multi-sample integration and batch correction

---

## References

1. Stahl, P.L. et al. (2016). Visualization and analysis of gene expression in tissue sections by spatial transcriptomics. *Science*, 353(6294), 78–82.
2. Luecken, M.D. & Theis, F.J. (2019). Current best practices in single-cell RNA-seq analysis. *Mol Syst Biol*, 15(6), e8746.
3. Traag, V.A. et al. (2019). From Louvain to Leiden: guaranteeing well-connected communities. *Sci Rep*, 9, 5233.
4. Palla, G. et al. (2022). Squidpy: a scalable framework for spatial omics analysis. *Nature Methods*, 19, 171–178.
5. Moran, P.A.P. (1950). Notes on continuous stochastic phenomena. *Biometrika*, 37(1/2), 17–23.
6. 10x Genomics Visium Spatial Gene Expression. https://www.10xgenomics.com/products/spatial-gene-expression

---

## Author

**Ravi Shanker Pandey, Ph.D.**
Computational Biologist | The Jackson Laboratory for Genomic Medicine
- 🔗 [LinkedIn](https://www.linkedin.com/in/ravi-pandey-50943980/)
- 🐙 [GitHub](https://github.com/pandeyravi15)
- 📧 ravisp1587@gmail.com
- 🎓 [Google Scholar](https://scholar.google.com/citations?user=FiCtSqQAAAAJ)

---

*This repository was developed as part of a computational biology portfolio demonstrating spatial transcriptomics expertise. All code is original and reproducible under the MIT License.*
