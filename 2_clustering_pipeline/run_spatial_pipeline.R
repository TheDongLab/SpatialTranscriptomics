source("./spatial_clustering_pipeline.R")

merged_seurat <- "comb_seurat_postqc.rds"
combined <- readRDS(merged_seurat)
resolutions_list <- c(0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6)
harmony_vars <- c("orig.ident", "sex", "expired_age", "PMI", "RIN")
harmony_thetas <- c(0.4, 0.4, 0.4, 0.4, 0.4)

random_seed <- 101
combined_looprun <- combined
job_name <- sprintf("all_data_5covar0p4_seed%s", random_seed)
cluster_pipeline(combined_looprun, job_name,
                 harmony_vars = harmony_vars,
                 harmony_thetas = harmony_thetas,
                 resolutions_list = resolutions_list,
                 random_seed = random_seed)

