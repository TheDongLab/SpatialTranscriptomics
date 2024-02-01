# Clustering and analysis pipeline for Visium 10X spatial transcriptomic data
Jie Yuan 1/8/2024

Code applied to Visium 10X spatial transcriptomic data in the Middle Temporal Gyrus (MTG), part of the ASAP PD5D Parkinson's Cell Atlas project. Analysis steps are split into numbered folders.


### 1_sample_QC
- st_sample_qc.R: Flag sample outliers for removal by aggregating per-sample count matrices and clustering.
- post_qc_all_data_to_seurat.R: Combine all sample count matrices and metadata into a single Seurat object.

### 2_clustering_pipeline
- spatial_clustering_pipeline.R: Runs the main clustering pipeline. This includes the following:
    - variable feature selection
    - clustering using GLM-PCA (2019)
    - batch correction using Harmony (2019)
    - Neighborhood detection using Leiden algorithm scanning across resolution values
- run_spatial_pipeline.R: Runs the clustering pipeline specified in the above script.

### 3_post_clustering_validation
- validate_best_layers_by_resolution.R: In middle temporal gyrus, compares cluster annotations across resolution values to cortical layer spatial data from Maynard et al., 2021. A permutation function identifies the best cluster-to-layer assignment for each resolution.
- validate_new_layer4.R: Layer 4 is extracted by sub-clustering Layer 5 and extracting one of its sub-clusters. This reruns the validation with the new Layer 4 annotation.

### 4_pseudobulk_DE
- pseudobulk_DE_DESeq2.R: Performs pseudobulk DE analysis using DESeq2 on each cortical layer independently. Different types of regression variables are demonstrated, such as quantitative variables, factors, and ordered factors.


## System requirements
- R (4.2.3)
- Seurat (4.3.0)
- SeuratWrappers (0.3.1)
- harmony (1.0.3)
- glmpca (0.2.0)
- DESeq2 (1.38.3)
- dplyr (1.1.2)
- ggplot2 (3.4.2)
- testit (0.13)
- leidenalg (0.10.1, use conda with reticulate)


## Maintenance and Contribution
This code is developed and maintained by Jie Yuan and Xianjun Dong. Please email jyuan15(AT)bwh.harvard,edu OR xdong(AT)bwh.harvard.edu for questions and comments.
