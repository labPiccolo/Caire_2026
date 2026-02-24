# Caire et al., 2026

This directory contains metadata and scripts used in `Caire et al. 2026`, required to replicate the analysis described in the paper.

## 1. Contents

### 1.1 Template script for the analysis of scRNA-seq samples

[scRNAseq_example_workflow.R](scRNAseq_example_workflow.R) script provides an end-to-end single-sample scRNA-seq preprocessing and QC workflow. It reads CellRanger output, performs QC and filtering (including Scrublet doublet detection), normalizes data, detects highly variable features, scales and runs PCA, computes clustering and dimensionality reductions (tSNE/UMAP), generates marker and QC plots, runs `SingleR` annotation. For human data only, the script proceeds with `infercnv` for detection of malignant cells and tumor-subset reclustering. After the analysis, the tumor clusters for every single human sample, cleaned of residual contaminations, were marked as **tumor subset**; tumor subset cells are reported in [tumor subsets metadata](Tumor_subsets_metadata.csv) for reproducibility.  
**Important!** This template script is intended for human scRNA-seq samples. For mouse samples the gene symbols used through the analysis must be changed accordingly. Every other modification required to work on mouse data is reported in the comments. 

- Inputs: raw CellRanger output directory (containing raw data as downloaded by GEO) and sample identifier (`sample_name`), Scrublet Python executable and script paths (`RETICULATE_PYTHON`, `scrublet_path`), and gene position file for `infercnv`.
- Outputs: preprocessed Seurat objects and diagnostic PDFs placed under `prepro/`, `clustering/`, `inferCNV/`, and `tumor_subset/` (QC plots, variable-gene plots, PCA/UMAP/tSNE, violin/dotplots, PC heatmaps, marker plots); inferCNV results are saved in `inferCNV/`.
- Implementation notes: uses `popsicleR` (v 0.2.1) for initial QC, `Seurat` (v 3.1.5) for normalization, feature selection, scaling, PCA, clustering and plotting, `SingleR` (v 1.0.6) for annotation, `Scrublet` (via `reticulate`) for doublet detection, and `infercnv` (v 1.2.1) for CNV analysis; the workflow is parallelized with `future` and written for `R` (v3.6.3).  

### 1.2 Integration of scRNA-seq samples

#### 1.2.1 Integration of human normal mammary glands

[hMG_reference_integration.R](hMG_reference_integration.R) script integrates the three normal mammary gland scRNA-seq samples ([hMG1](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8621910), [hMG2](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8621911), and [hMG3](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8621912)) to construct an integrated human reference. To ensure reproducibility, epithelial cells used in this integration, for each of the samples, are specified in [hMG metadata](hMG_epi_metadata.csv).

- Inputs: raw Seurat objects for the three samples (`hMG1.rds`, `hMG2.rds`, `hMG3.rds`), and epithelial-cell barcode lists used to subset each sample (`hMG1_epi_cells.rds`, `hMG2_epi_cells.rds`, `hMG3_epi_cells.rds`).
- Outputs: integrated Seurat object `hMG_epithelialIntegration.rds`; preprocessed integrated object `hMG_PD_epithelialIntegration_prepro.rds`; diagnostic PDF plots written to `preprocessing/` (dimensional reductions, feature plots, PC heatmaps) and `clustering/` (cluster/marker plots and sample distribution PDFs) directories.
- Implementation notes: uses `Seurat` CCA-based integration (`FindIntegrationAnchors` + `IntegrateData`, *dims=1:30*, *anchor.features=2000*); the script is written for `R` (v3.6.3) and requires `Seurat`, `ggplot2`, `dplyr`, `cowplot`, and `future`.

#### 1.2.2 Integration of mouse developmental mammary glands

[mouse_MG_development_integration.R](mouse_MG_development_integration.R)  script integrates five scRNA-seq samples of mouse MG (all samples from [GSE111113](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE111113), and puberty sample replicate 1; [GSM2759554](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2759554) from [GSE103272](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103272)) to construct a reference for mouse MG development.

- Input: raw Seurat objects for the five samples (`E16.rds`, `E18.rds`, `P4.rds`, `WK5.rds`, `WK12.rds`), and epithelial-cell barcode lists used to subset each sample (`E16_epi_cells.rds`, `E18_epi_cells.rds`, `P4_epi_cells.rds`, `WK5_epi_cells.rds`, `WK12_epi_cells.rds`).
- Outputs: integrated Seurat object `mouseMG_epithelial_Integration.rds`; preprocessed integrated object `mouseMG_epithelial_Integration_prepro.rds`; diagnostic PDF plots written to `preprocessing/` and `clustering/` directories.
- Implementation notes: uses `Seurat` CCA-based integration (`FindIntegrationAnchors` + `IntegrateData`, *dims=1:30*, *anchor.features=2000*); the script is written for `R` (v3.6.3) and requires `Seurat`, `ggplot2`, `dplyr`, `cowplot`, and `future`.

### 1.3 Differential expression analysis of human tumor pseudobulks

[DEGs_pseudobulk_Met-vs-Primary.R](DEGs_pseudobulk_Met-vs-Primary.R) script finds differentially expressed genes (**DEGs**) between primary and metastatic tumor pseudobulks. The script loads per-sample Seurat objects; merges and renames cells by sample; aggregates counts into pseudobulks per sample; filters low-expressed genes; runs edgeR (normalization, dispersion, GLM) to test **metastasis** vs **primary**; extracts significant genes and prepares summary tables. The pseudobulk were generated using *aggregate.Matrix* function from the `Matrix.utils` R package (version 0.9.8), with sample as the grouping variable and *sum* as the aggregation function. For differential expression computation, standard `edgeR` pipeline was followed. To ensure reproducibility, tumor cells used for this analysis were reported in the [tumor subsets metadata](Tumor_subsets_metadata.csv).

- Input: list of all the 24 analyzed samples, as Seurat objects, containing only tumor cells and a custom label indicating if the sample is a **metastasis** or a **primary** tumor - `sample_list.rds`.
- Outputs: lists of DEGs (`edgeR-Metastasis_vs_Primary.txt` and `edgeR-Metastasis_vs_Primary-SIG.txt`) and PDF plots (violin plots, heatmap, volcano plot)
- Implementation notes: uses `Seurat` (v3.1.5), `Matrix.utils` (v0.9.8) for aggregation, `edgeR` (v3.28.1) for DEGs computation, and `pheatmap`/`ggplot2` for figures; written for `R` (v3.6.3).

### 1.4 Analysis of 2D morphology in TNBC H&E-derived tumor areas

[2D_morphology_TNBC_tumors.R](2D_morphology_TNBC_tumors.R) script quantifies 2D tumor morphology from QuPath-exported GeoJSON files. The script reads per-cell annotations, extracts tumor-labeled objects, unions/simplifies and buffers polygons (tested with different buffer sizes, but the final buffer size chosen in the paper is 50 pixels), and computes morphology metrics such as area, perimeter, polygon counts, polygons per area/cell, and perimeter per area. It also joins metrics with clinical annotations for comparisons.

- Inputs: per-sample json objects (.geojson) containing Stardist segmentation masks and cell labels identifying tumor cells.
- Outputs: per-sample polygon images (PNG/PDF), `polygon_metrics.csv`, and `polygon_metrics_boxplots.pdf`.
- Implementation: written for `R` (v4.4.3) using `sf`, `ggplot2`, `dplyr`, `data.table`, `future` (multisession), and related packages.

### 1.5 Metadata of scRNA-seq data

#### 1.5.1 Metadata of human tumor subsets

[tumor subsets metadata](Tumor_subsets_metadata.csv) is a simple `csv` two-column table (`sample`, cell `barcode`). Each row lists a cell barcode assigned to a tumor subset for a given sample (e.g., *BBM5*, *AAACCCAAGACATCCT*). Used to reproduce tumor-subset selections and downstream analyses (pseudobulks, DE tests, plots).

#### 1.5.2 Metadata of human integrated normal mammary gland

[hMG metadata](hMG_epi_metadata.csv) is a simple `csv` three-column table (`sample`, cell `barcode`, `epithelial annotation`). Each row lists, for a given epithelial cell, the sample of origin, the barcode assigned and the epithelial class **L1**, **L2**, and **Basal**. Used to reproduce the inputs for [hMG_reference_integration.R](hMG_reference_integration.R) and the differential analysis involving normal epithelial cell types.

## 2. Code reproducibility

2D morphology analysis in TNBC H&E-derived tumor areas was performed using `R` (v4.4.3) and `sf` package (v1.0-21). Tumor areas were identified using a Stardist-based weakly supervised machine learning object classifier in `QuPath` (v0.5.0).  

All the other scripts were performed in `R` (v3.6.3), using `Seurat` (v3.1.5), `popsicleR` (v0.2.1), `scrublet` (v0.2.1), `SingleR` (v1.0.6), `infercnv` (v1.2.1).  

For pseudobulk analysis `Matrix.utils` R package (v0.9.8) was used to aggregate scRNA-seq data in pseudobulk, and `edgeR` (v3.28.1) was used to compute differential expression analysis  

For more details refer to methods section of the paper.

## 3. Citation

`Caire et al., 2026` paper (in press):
**A 3D Morphogenetic Program for Metastatic Outgrowth in Breast Cancer**.

Robin Caire, Roberta Bordo, Francesca Zanconato, Tito Panciera, Paolo Contessotto, Mikaela Fakiola, Ramona Bason, Oriana Romano, Ambela Suli, Giusy Battilana, Matteo Marchionni, Mattia Forcato, Estelle Audoux, Sara Donzelli, Maria Vittoria Dieci, Gaia Griguolo, Maria AntoniaCarosi, Matteo Fassan, Vincenza Guzzardo, Angelo Paolo Dei Tos, Silvia Marsoni, Pei-Hsun Wu, Denis Wirtz, Shanshan He, Cecilia Casali, Francesco Volpin, Giovanni Blandino, Claudio Tripodo, Silvio Bicciato, Valentina Guarneri, Massimiliano Pagani, Michelangelo Cordenonsi and Stefano
Piccolo
