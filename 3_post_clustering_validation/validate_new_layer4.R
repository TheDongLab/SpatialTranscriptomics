# Redo the validation, but now having identified cortical Layer 4 as a subcluster of Layer 5.

library(Seurat)
library(SeuratWrappers)
library(ggplot2)
library(reshape2)
# library(gridExtra)
# load main clustering script
source("../2_clustering_pipeline/spatial_clustering_pipeline.R")

# updated for panda server
# seed 101, select res = 0.45
comb_seurat <- readRDS("comb_seurat_diet_out.rds")
comb_seurat@meta.data$num_label <- comb_seurat@meta.data$Spatial_snn_res.0.45
comb_seurat@meta.data$layer_label_v1 <- "Noise"
comb_seurat@meta.data[comb_seurat@meta.data$num_label == 7, "layer_label_v1"] <- "Layer 1"
comb_seurat@meta.data[comb_seurat@meta.data$num_label == 5, "layer_label_v1"] <- "Layer 2"
comb_seurat@meta.data[comb_seurat@meta.data$num_label == 2, "layer_label_v1"] <- "Layer 3"
comb_seurat@meta.data[comb_seurat@meta.data$num_label == 1, "layer_label_v1"] <- "Layer 5"
comb_seurat@meta.data[comb_seurat@meta.data$num_label == 4, "layer_label_v1"] <- "Layer 6"
comb_seurat@meta.data[comb_seurat@meta.data$num_label == 26, "layer_label_v1"] <- "WM_1"
comb_seurat@meta.data[comb_seurat@meta.data$num_label == 6, "layer_label_v1"] <- "WM_2"
comb_seurat@meta.data[comb_seurat@meta.data$num_label == 3, "layer_label_v1"] <- "WM_3"
comb_seurat@meta.data$num_label <- NULL
comb_seurat@meta.data$layer_label_v1 <- as.factor(comb_seurat@meta.data$layer_label_v1)
Idents(comb_seurat) <- comb_seurat@meta.data$layer_label_v1

layers_v1 <- sort(unique(comb_seurat@meta.data$layer_label_v1))
layers_v1 <- layers_v1[ layers_v1 != "Noise" ]
resolutions_list <- c(0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4)
harmony_vars <- c("orig.ident", "sex", "expired_age", "PMI", "RIN")
harmony_thetas <- c(0.4, 0.4, 0.4, 0.4, 0.4)
random_seed <- 101
for( layer in layers_v1 ) {
    print(paste("Subclustering", layer))
    comb_sub <- subset(comb_seurat, subset = layer_label_v1 == layer)
    job_name <- sprintf("seurat_subclust_%s_seed%s", gsub(' ', '', layer), random_seed)
    cluster_pipeline(comb_sub, job_name,
                     harmony_vars = harmony_vars,
                     harmony_thetas = harmony_thetas,
                     resolutions_list = resolutions_list,
                     random_seed = random_seed)
}


# extract Layer 4 from Layer 5; this is stored in layer_label_v2
layer <- "Layer 5"
subcluster_res <- 0.2
cluster_col <- paste0("Spatial_snn_res.", subcluster_res)
comb_sub <- readRDS(sprintf("comb_seurat_diet_out.rds", gsub(' ','',layer)))
Idents(comb_sub) <- comb_sub@meta.data[, cluster_col]
comb_sub@meta.data[, "layer_label_v2"] <- as.character(comb_sub@meta.data[, cluster_col])
comb_sub@meta.data[comb_sub@meta.data$layer_label_v2 == 3, "layer_label_v2"] <- "Layer 4"
comb_sub@meta.data[comb_sub@meta.data$layer_label_v2 == 1, "layer_label_v2"] <- "Layer 5"
comb_sub@meta.data[comb_sub@meta.data$layer_label_v2 == 2, "layer_label_v2"] <- "Layer 5"
comb_sub@meta.data[comb_sub@meta.data$layer_label_v2 == 4, "layer_label_v2"] <- "Layer 5"
comb_sub@meta.data[comb_sub@meta.data$layer_label_v2 == 5, "layer_label_v2"] <- "Layer 5"

comb_seurat@meta.data$Row.names <- rownames(comb_seurat@meta.data)
comb_sub@meta.data$Row.names <- rownames(comb_sub@meta.data)
temp <- merge(comb_seurat@meta.data, comb_sub@meta.data[, c("Row.names", "layer_label_v2")], by="Row.names", all.x=TRUE)
rownames(temp) <- temp$Row.names
temp$Row.names <- NULL
temp <- temp[rownames(comb_seurat@meta.data), ]
temp[is.na(temp$layer_label_v2), "layer_label_v2"] <- as.character(temp[is.na(temp$layer_label_v2), "layer_label_v1"])
temp$layer_label_v2 <- as.factor(temp$layer_label_v2)
comb_seurat@meta.data <- temp



comb_seurat@meta.data$sublayer_label <- NA
subcluster_res <- 0.1




# color palette from:
# https://colorswall.com/palette/142
group.colors <- c("Layer 1"="#f43545", "Layer 2"="#ff8901", "Layer 3"="#fad717",
                  "Layer 4"="#00ba71", "Layer 5"="#00c2de", "Layer 6"="#ba55d3",
                  "WM_1"="#3b2175", "WM_2"="#75215b", "WM_3"="#652175", "Noise"="#c8c8c8")
pdf("recolored_layers_plusnewL4.pdf")
for (img in names(comb_seurat@images)) {
    meta_row <- rownames(comb_seurat@images[[img]]@coordinates)[1]
    sample_id <- comb_seurat@meta.data[meta_row, "orig.ident"]
    diagnosis <- comb_seurat@meta.data[meta_row, "diagnosis"]

    sp_pt <- SpatialPlot(comb_seurat, group.by="layer_label_v2", image.alpha=0.3, images=img, stroke=NA) + scale_fill_manual(values=group.colors) + ggtitle(paste(sample_id, diagnosis))
    print(sp_pt)
}
dev.off()

saveRDS(comb_seurat, "comb_seurat_diet_out_plusL4.rds")


generate_layer_corr_table <- function(top_markers, maynard, maynard_t_stats) {
    corr_table <- as.data.frame(matrix(ncol=length(unique(top_markers$cluster)),
                                       nrow=length(maynard_t_stats)))
    rownames(corr_table)<-maynard_t_stats
    colnames(corr_table)<-unique(top_markers$cluster)

    corr_table_pvals <- as.data.frame(matrix(ncol=length(unique(top_markers$cluster)),
                                             nrow=length(maynard_t_stats)))
    rownames(corr_table_pvals)<-maynard_t_stats
    colnames(corr_table_pvals)<-unique(top_markers$cluster)

    for (i in unique(top_markers$cluster)) {
        sig_genes_layer <- subset(top_markers, cluster==i)
        for(j in maynard_t_stats) {
            merged_mat <- merge(sig_genes_layer, maynard[, c(j, "gene")], by="gene")
            # res <- cor.test(merged_mat$avg_log2FC, merged_mat[[j]], method="spearman")
            res <- cor.test(merged_mat$avg_log2FC, merged_mat[[j]], method="pearson")
            corr_table[j,i] <- res$estimate
            corr_table_pvals[j,i] <- res$p.value
        }
    }
    results <- list()
    results$pearson <- corr_table
    results$pvalues <- corr_table_pvals
    return(results)
}

permute_k_subset <- function(set, k) {
    permut_mat <- c()
    permute_k_helper <- function(chain, remaining, k) {
        if ( length(chain) == k | length(remaining) == 0 ) {
            # print(chain)
            permut_mat <<- rbind(permut_mat, chain)
        } else {
            for (i in 1:length(remaining)) {
                permute_k_helper(c(chain, remaining[i]), remaining[-i], k)
            }
        }
    }
    permute_k_helper(c(), set, k)
    rownames(permut_mat) <- NULL
    return(permut_mat)
}



maynard <- read.csv("maynard_layer_markers.csv", sep=',', header=TRUE)
maynard_t_stats <- c("t_stat_WM", "t_stat_Layer1", "t_stat_Layer2", "t_stat_Layer3",
                     "t_stat_Layer4", "t_stat_Layer5", "t_stat_Layer6")

set.seed(101)

# don't use every gene for correlation
Idents(comb_seurat_all) <- comb_seurat_all@meta.data[, "new_label"]
top_markers <- RunPrestoAll(comb_seurat_all, assay='Spatial', slot='data', logfc.threshold=0.1, return.thresh=0.1, min.pct = 0.1)

write.table(top_markers, "top_markers_newlayer4_all_data.csv", sep = ',', row.names = TRUE, col.names = NA)

corr_table <- generate_layer_corr_table(top_markers, maynard, maynard_t_stats)
write.table(corr_table$pearson, "maynard_validation_pearson_newlayer4_all_data.csv", sep = ',', row.names = TRUE, col.names = NA)
write.table(corr_table$pvalues, "maynard_validation_pvalues_newlayer4_all_data.csv", sep = ',', row.names = TRUE, col.names = NA)




# compare to the old labels (w/out Layer 4)
Idents(comb_seurat_all) <- comb_seurat_all@meta.data[, "layer_label"]
top_markers <- RunPrestoAll(comb_seurat_all, assay='Spatial', slot='data', logfc.threshold=0.1, return.thresh=0.1, min.pct = 0.1)

write.table(top_markers, "top_markers_oldlabels_all_data.csv", sep = ',', row.names = TRUE, col.names = NA)

corr_table <- generate_layer_corr_table(top_markers, maynard, maynard_t_stats)
write.table(corr_table$pearson, "maynard_validation_pearson_oldlabels_all_data.csv", sep = ',', row.names = TRUE, col.names = NA)
write.table(corr_table$pvalues, "maynard_validation_pvalues_oldlabels_all_data.csv", sep = ',', row.names = TRUE, col.names = NA)

