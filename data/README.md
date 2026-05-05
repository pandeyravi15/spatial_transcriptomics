# Data Download Instructions

## 10x Genomics Visium Human Brain FFPE

### Option 1 — Download via 10x Genomics Portal (Recommended)

1. Go to: https://www.10xgenomics.com/datasets/human-brain-cancer-11-mm-capture-area-ffpe-2-standard
2. Download the **SpaceRanger output** folder (filtered_feature_bc_matrix + spatial/)
3. Place in: `data/visium_human_brain/`

### Option 2 — Commandline download (curl)

```bash
# Create data directory
mkdir -p data/visium_human_brain

# Download filtered feature matrix
curl -O https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/Visium_Human_Breast_Cancer/Visium_Human_Breast_Cancer_filtered_feature_bc_matrix.h5

# Download spatial coordinates
curl -O https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/Visium_Human_Breast_Cancer/Visium_Human_Breast_Cancer_spatial.tar.gz
tar -xzf Visium_Human_Breast_Cancer_spatial.tar.gz
```

### Option 3 — Load directly via scanpy

```python
import scanpy as sc

# Load from SpaceRanger output directory
adata = sc.read_visium("data/visium_human_brain/")
print(adata)
# AnnData object with n_obs × n_vars = 3493 × 33538
```

### Expected directory structure after download

```
data/
└── visium_human_brain/
    ├── filtered_feature_bc_matrix/
    │   ├── barcodes.tsv.gz
    │   ├── features.tsv.gz
    │   └── matrix.mtx.gz
    └── spatial/
        ├── tissue_positions.csv
        ├── scalefactors_json.json
        ├── tissue_hires_image.png
        └── tissue_lowres_image.png
```

## Other Public Visium Datasets

| Dataset | Tissue | Link |
|---------|--------|------|
| Human Brain (FFPE) | Cortex, White Matter | 10x Genomics Portal |
| Mouse Brain (Coronal) | Full brain section | 10x Genomics Portal |
| Human Heart | Cardiac tissue | 10x Genomics Portal |
| DLPFC (Human) | Dorsolateral prefrontal cortex | spatialLIBD package |

## DLPFC Dataset (Highly recommended for benchmarking)

The DLPFC dataset (Maynard et al. 2021, *Nature Neuroscience*) has manually annotated cortical layers — ideal for validating clustering:

```r
# In R
BiocManager::install("spatialLIBD")
library(spatialLIBD)
spe <- fetch_data(type = "spe")
```
