Clustering and analysis pipeline for Visium 10X spatial transcriptomic data
Jie Yuan 1/8/2024

Analysis steps are split into numbered folders.


1_sample_QC
- st_sample_qc.R: Flag sample outliers for removal by aggregating per-sample count matrices and clustering.
- post_qc_all_data_to_seurat.R: Combine all sample count matrices and metadata into a single Seurat object.

2_clustering_pipeline
- spatial_clustering_pipeline.R: Runs the main clustering pipeline. This includes the following:
    - variable feature selection
    - clustering using GLM-PCA (2019)
    - batch correction using Harmony (2019)
    - Neighborhood detection using Leiden algorithm scanning across resolution values
- run_spatial_pipeline.R: Runs the clustering pipeline specified in the above script.

3_post_clustering_validation
- validate_best_layers_by_resolution.R: In middle temporal gyrus, compares cluster annotations across resolution values to cortical layer spatial data from Maynard et al., 2021. A permutation function identifies the best cluster-to-layer assignment for each resolution.
- validate_new_layer4.R: Layer 4 is extracted by sub-clustering Layer 5 and extracting one of its sub-clusters. This reruns the validation with the new Layer 4 annotation.

4_pseudobulk_DE
- pseudobulk_DE_DESeq2.R: Performs pseudobulk DE analysis using DESeq2 on each cortical layer independently. Different types of regression variables are demonstrated, such as quantitative variables, factors, and ordered factors.
