# Single-Nucleus RNA-seq Analysis of Hepatoblastoma      
   
This project analyzes single-nucleus RNA-seq (snRNA-seq) data from hepatoblastoma tumors, PDX models, and normal samples using the Seurat and Harmony frameworks. The pipeline includes data preprocessing, batch correction, clustering, manual annotation using GeneCards, differential expression analysis, and cell composition visualization.

All steps are reproducible via Docker.  
 
---   

## Project Structure

```
snrna_hb_project/
├── data/                        # Raw 10X filtered_feature_bc_matrix folders (user-provided)
├── results/                     # Output figures, heatmaps, DEGs
├── scripts/
│   ├── 01_initial_qc_normalization_harmony.R
│   ├── 02_annotation_and_differential_expression.R
│   ├── 03_celltype_composition_analysis.R
│   
├── envs/
│   └── Dockerfile
├── .dockerignore
└── README.md
```

---

## Setup

### Requirements (if not using Docker)

- R version 4.3.1
- R packages:
  - Seurat (v5.3.0)
  - harmony (v1.2.3)
  - ggplot2 (v3.5.2)
  - dplyr (v1.1.4)
  - patchwork (v1.3.0)
  - cowplot (v1.1.3)
  - tidyverse (v2.0.0)
  - EnhancedVolcano

---

## Usage

### 1. Prepare Input
  
**Important**: Update the following line in `01_initial_qc_normalization_harmony.R` to point to your data location:

```r
basedir <- "/your/absolute/path/to/data"
```

### 2. Run Scripts Manually

```r
source("scripts/01_initial_qc_normalization_harmony.R")
source("scripts/02_annotation_and_differential_expression.R")
source("scripts/03_celltype_composition_analysis.R")
```

---

## Reproducibility with Docker

This project includes a `Dockerfile` that builds a complete R environment with all required packages pinned to exact versions.

### Build the Docker Image

```bash
docker build -f envs/Dockerfile -t snrna-seurat .
```

### Run the Pipeline

```bash
docker run --rm \
  -v $(pwd)/data:/app/data \
  -v $(pwd)/results:/app/results \
  snrna-seurat
```

This will run `run_all.R` and save output plots and DE tables to `results/`.

---

## Results Summary

### Clustering and Annotation

- 10 clusters were identified using unsupervised clustering (Seurat + Harmony) on UMAP embeddings.
- Top markers were extracted using `FindAllMarkers()` and visualized with `DoHeatmap()`.
- Marker genes were cross-referenced with **GeneCards** to determine tissue-specific expression.
- Clusters were manually annotated to represent diverse cell states in hepatoblastoma and its microenvironment.

### Focused Differential Expression (Cluster 1)

- Cluster 1 was annotated as **fetal-like hepatoblastoma**, based on high expression of **AFP**, **SLC22A9**, and **CYP3A7**.
- Within this cluster, **tumor** vs **PDX** cells were compared using `FindMarkers()`.

#### Reference:
- Biological interpretation is guided by [PMID: 34497364](https://pubmed.ncbi.nlm.nih.gov/34497364), which describes transcriptomic shifts between fetal hepatocytes, tumor tissue, and xenograft models.

---

## Interpretation

- **Tumor cells** are exposed to immune cells, stress signals, and treatment pressures. As a result, they upregulate genes related to:
  - **Immune signaling** (e.g., **CFH** – complement system)
  - **Stress/drug metabolism** (e.g., **CYP3A5**)

- **PDX cells**, grown in immune-deficient mice, retain a purer fetal-like liver identity:
  - Higher expression of **AFP**, a classical fetal liver marker
  - Lower activation of immune/stress pathways

### Key Findings

- **CFH**: Complement-related gene, more expressed in tumor cells under immune pressure.
- **CYP3A5**: Drug-metabolizing enzyme, may be induced by tumor environment.
- **AFP**: Elevated in PDX cells, supporting fetal identity preservation.

### Biological Inference

- **Tumor samples** reflect a more complex and reactive microenvironment.
- **PDX samples** offer a simplified, fetal-like baseline, valuable for controlled experimental models but lacking immune/stress dynamics.
- Therefore, **PDX models are useful**, but **caution is needed** when interpreting immune-related pathways or therapy responses.

---
## Citation & Acknowledgments

- Seurat: Hao et al., Cell 2021.
- Harmony: Korsunsky et al., Nat Methods 2019.
- GeneCards: Weizmann Institute of Science.
- Data interpretation guided by [PMID: 34497364](https://pubmed.ncbi.nlm.nih.gov/34497364).

---

