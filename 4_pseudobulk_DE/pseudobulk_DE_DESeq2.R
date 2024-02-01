# Pseudobulk spatial cortical samples by layer and perform DE analysis using DESeq2
# Though it is somewhat redundant, this script illustrates different types of regression
# variables: quantitative, case/control factor, and multi-level ordinal factor.

library(Seurat)
library(DESeq2)
library(stringr)
library(testit)

comb_seurat <- readRDS("comb_seurat.rds")

random_seed <- 101

# Map diagnosis to pseudotime: Control = 0, ILBD = 1, Case = 2
comb_seurat@meta.data$diagnosis.int <- NA
comb_seurat@meta.data$diagnosis.int[comb_seurat@meta.data$diagnosis == "Control"] <- 0
comb_seurat@meta.data$diagnosis.int[comb_seurat@meta.data$diagnosis == "ILBD"] <- 1
comb_seurat@meta.data$diagnosis.int[comb_seurat@meta.data$diagnosis == "Case"] <- 2
comb_seurat@meta.data$diagnosis.int <- as.numeric(comb_seurat@meta.data$diagnosis.int)

layer_annot <- "smoothed_label_s5"
layers <- sort(unique(comb_seurat@meta.data[, layer_annot]))
features_quantitative <- c("nctx_temporal", "brain_stem_sn", "last_mmse_test_score", "motor_updrs_score")
min_layer_spot_count <- 20

# loop over layers: create seurat object for each
for(layer in layers) {
    temp_layer <- FetchData(comb_seurat, vars = layer_annot)
    comb_sub <- comb_seurat[, which(temp_layer == layer)]
    Idents(object = comb_sub) <- "sample_name"
    
    # get countData and filter out lowly expressed genes
    countData <- AggregateExpression(comb_sub, assays="Spatial", slot="counts")
    countData <- countData$Spatial
    keep_genes <- rowSums(countData >= 5) >= 3
    countData <- countData[keep_genes, ]
    
    # filter out samples with fewer than N spots belonging to the given layer
    keep_samples <- names( which( table(comb_sub@meta.data$sample_name) > min_layer_spot_count ) )
    countData <- countData[, keep_samples]

    count_out_csv <- sprintf("agg_counts_filtered_%s_%s.csv", layer_annot, gsub(' ', '', layer))
    if(!file.exists(count_out_csv)) {
        write.csv(countData, file = count_out_csv)
    }
    
    # for each feature, create metadata table
    for(feat_quant in features_quantitative) {
        print(paste("Running", layer, feat_quant))
        
        sample_meta <- comb_sub@meta.data[, c("sample_name", "batch", "sex", "expired_age_num", "PMI_num", "RIN_num", feat_quant)]
        sample_meta <- sample_meta[!duplicated(sample_meta$sample_name), ]
        sample_meta <- sample_meta[sample_meta$sample_name %in% colnames(countData), ]
        # filter out samples whose feature label is NA
        sample_meta <- sample_meta[!is.na(sample_meta[, feat_quant]), ]
        rownames(sample_meta) <- sample_meta$sample_name
        countData_filt <- countData[, as.character(sample_meta$sample_name)]
        print(paste("Samples passing filtering", ncol(countData_filt)))
        
        assert(all(rownames(sample_meta) == colnames(countData_filt)))
        
        # DESeq2 likes to have all numeric variables centered and scaled
        sample_meta$sample_name <- as.factor(sample_meta$sample_name)
        sample_meta$batch <- as.factor(sample_meta$batch)
        sample_meta$sex <- as.factor(sample_meta$sex)
        sample_meta$PMI_num <- scale(sample_meta$PMI_num, center=TRUE, scale=TRUE)
        sample_meta$RIN_num <- scale(sample_meta$RIN_num, center=TRUE, scale=TRUE)
        sample_meta$expired_age_num <- scale(sample_meta$expired_age_num, center=TRUE, scale=TRUE)
        sample_meta[, feat_quant] <- scale(sample_meta[, feat_quant], center=TRUE, scale=TRUE)
        
        assert(is.factor(sample_meta$sample_name))
        assert(is.factor(sample_meta$batch))
        assert(is.factor(sample_meta$sex))
        assert(is.numeric(sample_meta$expired_age_num))
        assert(is.numeric(sample_meta$PMI_num))
        assert(is.numeric(sample_meta$RIN_num))
        assert(is.numeric(sample_meta[, feat_quant]))
        
        
        tryCatch(
            expr = {

                set.seed(random_seed)
                dds <- DESeqDataSetFromMatrix(countData=countData_filt, colData=sample_meta, design= model.matrix(~ sample_meta$batch + sample_meta$PMI_num + sample_meta$RIN_num + sample_meta$sex + sample_meta$expired_age_num + sample_meta[, feat_quant]), tidy=FALSE)
                saveRDS(dds, sprintf("pseudobulk_clinical_deseq2_dds_%s_%s.rds", feat_quant, gsub(' ', '', layer)))

                dds_lrt <- DESeq(dds, test="LRT", reduced = model.matrix(~ sample_meta$batch + sample_meta$PMI_num + sample_meta$RIN_num + sample_meta$sex + sample_meta$expired_age_num))
                saveRDS(dds_lrt, sprintf("pseudobulk_clinical_deseq2_dds_lrt_%s_%s.rds", feat_quant, gsub(' ', '', layer)))

                results_table <- results(dds_lrt)

                write.csv(results_table, sprintf("pseudobulk_clinical_deseq2_%s_%s.csv", feat_quant, gsub(' ', '', layer)))

            },
            error = function(e){
                # (Optional)
                # Do this if an error is caught...
                print(paste("Running", layer, feat_quant))
                message(e)
            }
        )

    }
}



# Control/ILBD and Control/Case comparisons

comb_seurat <- readRDS("comb_seurat.rds")

random_seed <- 101

# in nctx_temporal, there is only 1 sample with a score of 4: set this to 3 so as not to overly bias the regression
comb_seurat@meta.data$nctx_temporal_fix <- comb_seurat@meta.data$nctx_temporal
comb_seurat@meta.data$nctx_temporal_fix[comb_seurat@meta.data$nctx_temporal_fix > 3] <- 3

layer_annot <- "smoothed_label_s5"
layers <- sort(unique(comb_seurat@meta.data[, layer_annot]))
# features_quantitative <- c("nctx_temporal_fix", "brain_stem_sn", "last_mmse_test_score", "motor_updrs_score")
min_layer_spot_count <- 20

comb_seurat_original <- comb_seurat
# control vs ILBD
feat_diag <- "diagnosis_cont_ilbd"
comb_seurat <- comb_seurat[, which(comb_seurat@meta.data$diagnosis %in% c("Control", "ILBD"))]
comb_seurat@meta.data[, feat_diag] <- factor(comb_seurat@meta.data$diagnosis, levels=c("Control", "ILBD"))

# loop over layers: create seurat object for each
for(layer in layers) {
    temp_layer <- FetchData(comb_seurat, vars = layer_annot)
    comb_sub <- comb_seurat[, which(temp_layer == layer)]
    Idents(object = comb_sub) <- "sample_name"
    
    # get countData and filter out lowly expressed genes
    countData <- AggregateExpression(comb_sub, assays="Spatial", slot="counts")
    countData <- countData$Spatial
    keep_genes <- rowSums(countData >= 5) >= 3
    countData <- countData[keep_genes, ]
    
    # filter out samples with fewer than N spots belonging to the given layer
    keep_samples <- names( which( table(comb_sub@meta.data$sample_name) > min_layer_spot_count ) )
    countData <- countData[, keep_samples]

    # for each feature, create metadata table
    print(paste("Running", layer, feat_diag))
    
    sample_meta <- comb_sub@meta.data[, c("sample_name", "batch", "sex", "expired_age_num", "PMI_num", "RIN_num", feat_diag)]
    sample_meta <- sample_meta[!duplicated(sample_meta$sample_name), ]
    sample_meta <- sample_meta[sample_meta$sample_name %in% colnames(countData), ]
    # filter out samples whose feature label is NA
    sample_meta <- sample_meta[!is.na(sample_meta[, feat_diag]), ]
    rownames(sample_meta) <- sample_meta$sample_name
    countData_filt <- countData[, as.character(sample_meta$sample_name)]
    print(paste("Samples passing filtering", ncol(countData_filt)))
    
    assert(all(rownames(sample_meta) == colnames(countData_filt)))
    
    # DESeq2 likes to have all numeric variables centered and scaled
    sample_meta$sample_name <- as.factor(sample_meta$sample_name)
    sample_meta$batch <- as.factor(sample_meta$batch)
    sample_meta$sex <- as.factor(sample_meta$sex)
    sample_meta$PMI_num <- scale(sample_meta$PMI_num, center=TRUE, scale=TRUE)
    sample_meta$RIN_num <- scale(sample_meta$RIN_num, center=TRUE, scale=TRUE)
    sample_meta$expired_age_num <- scale(sample_meta$expired_age_num, center=TRUE, scale=TRUE)
    # sample_meta[, feat_quant] <- scale(sample_meta[, feat_quant], center=TRUE, scale=TRUE)
    
    assert(is.factor(sample_meta$sample_name))
    assert(is.factor(sample_meta$batch))
    assert(is.factor(sample_meta$sex))
    assert(is.numeric(sample_meta$expired_age_num))
    assert(is.numeric(sample_meta$PMI_num))
    assert(is.numeric(sample_meta$RIN_num))
    # assert(is.numeric(sample_meta[, feat_quant]))
    assert(is.factor(sample_meta[, feat_diag]))
    
    tryCatch(
        expr = {

            set.seed(random_seed)
            dds <- DESeqDataSetFromMatrix(countData=countData_filt, colData=sample_meta, design= model.matrix(~ sample_meta$batch + sample_meta$PMI_num + sample_meta$RIN_num + sample_meta$sex + sample_meta$expired_age_num + sample_meta[, feat_diag]), tidy=FALSE)
            saveRDS(dds, sprintf("pseudobulk_clinical_deseq2_dds_%s_%s.rds", feat_diag, gsub(' ', '', layer)))
            
            dds_lrt <- DESeq(dds, test="LRT", reduced = model.matrix(~ sample_meta$batch + sample_meta$PMI_num + sample_meta$RIN_num + sample_meta$sex + sample_meta$expired_age_num))
            saveRDS(dds_lrt, sprintf("pseudobulk_clinical_deseq2_dds_lrt_%s_%s.rds", feat_diag, gsub(' ', '', layer)))

            results_table <- results(dds_lrt)

            write.csv(results_table, sprintf("pseudobulk_clinical_deseq2_%s_%s.csv", feat_diag, gsub(' ', '', layer)))

        },
        error = function(e){
            # (Optional)
            # Do this if an error is caught...
            print(paste("Running", layer, feat_diag))
            message(e)
        }
    )

}




# diagnosis as ordered factor
comb_seurat <- readRDS("comb_seurat.rds")

random_seed <- 101

# map diagnosis to ordered factor
comb_seurat@meta.data$diagnosis <- factor(comb_seurat@meta.data$diagnosis, ordered=TRUE, levels=c("Control", "ILBD", "Case"))

layer_annot <- "smoothed_label_s5"
layers <- sort(unique(comb_seurat@meta.data[, layer_annot]))
features_ordered <- c("diagnosis")
min_layer_spot_count <- 20

# loop over layers: create seurat object for each
for(layer in layers) {
    temp_layer <- FetchData(comb_seurat, vars = layer_annot)
    comb_sub <- comb_seurat[, which(temp_layer == layer)]
    Idents(object = comb_sub) <- "sample_name"
    
    # get countData and filter out lowly expressed genes
    countData <- AggregateExpression(comb_sub, assays="Spatial", slot="counts")
    countData <- countData$Spatial
    keep_genes <- rowSums(countData >= 5) >= 3
    countData <- countData[keep_genes, ]
    
    # filter out samples with fewer than N spots belonging to the given layer
    keep_samples <- names( which( table(comb_sub@meta.data$sample_name) > min_layer_spot_count ) )
    countData <- countData[, keep_samples]

    count_out_csv <- sprintf("agg_counts_filtered_%s_%s.csv", layer_annot, gsub(' ', '', layer))
    if(!file.exists(count_out_csv)) {
        write.csv(countData, file = count_out_csv)
    }
    
    # for each feature, create metadata table
    for(feat_ordered in features_ordered) {
        print(paste("Running", layer, feat_ordered))
        
        sample_meta <- comb_sub@meta.data[, c("sample_name", "batch", "sex", "expired_age_num", "PMI_num", "RIN_num", feat_ordered)]
        sample_meta <- sample_meta[!duplicated(sample_meta$sample_name), ]
        sample_meta <- sample_meta[sample_meta$sample_name %in% colnames(countData), ]
        # filter out samples whose feature label is NA
        sample_meta <- sample_meta[!is.na(sample_meta[, feat_ordered]), ]
        rownames(sample_meta) <- sample_meta$sample_name
        countData_filt <- countData[, as.character(sample_meta$sample_name)]
        print(paste("Samples passing filtering", ncol(countData_filt)))
        
        assert(all(rownames(sample_meta) == colnames(countData_filt)))
        
        # DESeq2 likes to have all numeric variables centered and scaled
        sample_meta$sample_name <- as.factor(sample_meta$sample_name)
        sample_meta$batch <- as.factor(sample_meta$batch)
        sample_meta$sex <- as.factor(sample_meta$sex)
        sample_meta$PMI_num <- scale(sample_meta$PMI_num, center=TRUE, scale=TRUE)
        sample_meta$RIN_num <- scale(sample_meta$RIN_num, center=TRUE, scale=TRUE)
        sample_meta$expired_age_num <- scale(sample_meta$expired_age_num, center=TRUE, scale=TRUE)
        # sample_meta[, feat_quant] <- scale(sample_meta[, feat_quant], center=TRUE, scale=TRUE)
        
        
        assert(is.factor(sample_meta$sample_name))
        assert(is.factor(sample_meta$batch))
        assert(is.factor(sample_meta$sex))
        assert(is.numeric(sample_meta$expired_age_num))
        assert(is.numeric(sample_meta$PMI_num))
        assert(is.numeric(sample_meta$RIN_num))
        # assert(is.numeric(sample_meta[, feat_quant]))
        assert(is.factor(sample_meta[, feat_ordered]))
        assert(is.ordered(sample_meta[, feat_ordered]))
        
        tryCatch(
            expr = {

                set.seed(random_seed)
                dds <- DESeqDataSetFromMatrix(countData=countData_filt, colData=sample_meta, design= model.matrix(~ sample_meta$batch + sample_meta$PMI_num + sample_meta$RIN_num + sample_meta$sex + sample_meta$expired_age_num + sample_meta[, feat_ordered]), tidy=FALSE)
                saveRDS(dds, sprintf("pseudobulk_clinical_deseq2_dds_%s_%s.rds", feat_ordered, gsub(' ', '', layer)))

                dds_lrt <- DESeq(dds, test="LRT", reduced = model.matrix(~ sample_meta$batch + sample_meta$PMI_num + sample_meta$RIN_num + sample_meta$sex + sample_meta$expired_age_num))
                saveRDS(dds_lrt, sprintf("pseudobulk_clinical_deseq2_dds_lrt_%s_%s.rds", feat_ordered, gsub(' ', '', layer)))

                results_table <- results(dds_lrt)

                write.csv(results_table, sprintf("pseudobulk_clinical_deseq2_%s_%s.csv", feat_ordered, gsub(' ', '', layer)))

            },
            error = function(e){
                # (Optional)
                # Do this if an error is caught...
                print(paste("Running", layer, feat_ordered))
                message(e)
            }
        )

    }
}

