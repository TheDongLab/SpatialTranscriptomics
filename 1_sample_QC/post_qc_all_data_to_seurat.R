# Helper script to create a merged Seurat object of all samples.

library(Seurat)
library(data.table)
library(reshape2)
library(dplyr)
library(ggplot2)
library(stringr)
library(testit)

set.seed(101)

# initial spatial data samples file
input_file_list <- ""
clinical_features_list <- ""
output_data_file <- "st_files_postqc.csv"
output_postqc_seurat_dir <- "output_postqc_seurat_objs"
output_seurat_name <- "postqc"

st_input_files <- read.csv(input_file_list, comment.char = '#')
# remove QC-filtered samples according to dedicated column, then remove QC columns
st_input_files <- st_input_files[st_input_files$qc_final == "keep", ]
st_input_files$qc_final <- NULL
st_input_files$qc_threshold <- NULL
st_input_files$qc_worse_replicate <- NULL

st_input_files$orig.ident <- st_input_files$sample_name
st_input_files$sample_name <- str_split_fixed(st_input_files$orig.ident, '_', n=2)[, 1]
st_input_files$sample_name <- as.factor(st_input_files$sample_name)
st_input_files$batch <- as.factor(st_input_files$batch)
st_input_files$sex <- as.factor(st_input_files$sex)
st_input_files$diagnosis <- factor(st_input_files$diagnosis, levels=c("Control", "ILBD", "Case"))

clinical_data <- read.csv(clinical_features_list)
clinical_data$sample_name <- clinical_data$CaseID
clinical_data$CaseID <- NULL
# combine UPDRS score on/off medication
clinical_data$motor_updrs_score <- clinical_data$motor_updrs_score_off
clinical_data$motor_updrs_score[is.na(clinical_data$motor_updrs_score_off)] <- clinical_data$motor_updrs_score_on[is.na(clinical_data$motor_updrs_score_off)]
clinical_feature_set <- c("sample_name", "last_mmse_test_score",
                       "Unified.LB.Stage",
                       "PlaqueT", "TangleT",
                       "motor_updrs_score", "nctx_temporal", "brain_stem_sn")
clinical_data <- clinical_data[, clinical_feature_set]
rownames(clinical_data) <- clinical_data$CaseID
clinical_data$sample_name <- as.factor(clinical_data$sample_name)
clinical_data$Unified.LB.Stage <- str_split_fixed(clinical_data$Unified.LB.Stage, '\\.', n=2)[, 1]
clinical_data$Unified.LB.Stage[clinical_data$Unified.LB.Stage == "lla"] <- "ll"
clinical_data$Unified.LB.Stage[clinical_data$Unified.LB.Stage == "llb"] <- NA
# Unified.LB.Stage as numeric {0,1,2,3,4}
clinical_data$Unified.LB.Stage.int <- clinical_data$Unified.LB.Stage
clinical_data$Unified.LB.Stage.int[clinical_data$Unified.LB.Stage == ""] <- NA
clinical_data$Unified.LB.Stage.int[clinical_data$Unified.LB.Stage == "0"] <- 0
clinical_data$Unified.LB.Stage.int[clinical_data$Unified.LB.Stage == "l"] <- 1
clinical_data$Unified.LB.Stage.int[clinical_data$Unified.LB.Stage == "ll"] <- 2
clinical_data$Unified.LB.Stage.int[clinical_data$Unified.LB.Stage == "lll"] <- 3
clinical_data$Unified.LB.Stage.int[clinical_data$Unified.LB.Stage == "lV"] <- 4
clinical_data$Unified.LB.Stage.int <- as.numeric(clinical_data$Unified.LB.Stage.int)

temp <- merge(x=st_input_files, y=clinical_data, by="sample_name", all.x=TRUE)
assert(nrow(temp) == nrow(st_input_files))
st_input_files <- temp

write.table(st_input_files, file=output_data_file, sep='\t', row.names=TRUE)


ii <- 1
for (i in rownames(st_input_files)) {
   row_vals = st_input_files[i, ]
   if (row_vals$fold_removed) {
       print(paste(ii, row_vals$sample_name, "fold removed"), sep=', ')
       st_data <- readRDS(row_vals$dir)
       st_data <- subset(st_data, subset=selected_spot == FALSE)
   } else {
       print(paste(ii, row_vals$sample_name, "no fold removed"), sep=', ')
       st_data <- Load10X_Spatial(row_vals$dir)
   }
   st_data@meta.data["orig.ident"] <- row_vals$orig.ident
   st_data@meta.data["sample_name"] <- row_vals$sample_name
   st_data@meta.data["batch"] <- row_vals$batch
   st_data@meta.data["PMI"] <- row_vals$PMI
   st_data@meta.data["RIN"] <- row_vals$RIN
   st_data@meta.data["sex"] <- row_vals$sex
   st_data@meta.data["expired_age"] <- row_vals$expired_age
   st_data@meta.data["diagnosis"] <- row_vals$diagnosis
   st_data@meta.data["last_mmse_test_score"] <- row_vals$last_mmse_test_score
   st_data@meta.data["motor_updrs_score"] <- row_vals$motor_updrs_score
   st_data@meta.data["nctx_temporal"] <- row_vals$nctx_temporal
   st_data@meta.data["brain_stem_sn"] <- row_vals$brain_stem_sn

   ii <- ii + 1
   
   saveRDS(st_data, file.path(output_postqc_seurat_dir, paste0(row_vals$sample_name, "_", output_seurat_name, ".rds")))
}

data.list <- list()
sample.list <- c()
ii <- 1
files <- list.files(path="/postqc/seurat/obj/dir", pattern="*.rds", full.names=TRUE, recursive=FALSE)
for(st_file in files) {
    st_data <- readRDS(st_file)
    data.list[[ii]] = st_data
    ii <- ii + 1
    print(st_file)
    sample_name <- as.character(st_data@meta.data$sample_name[1])
    expired_age <- st_data@meta.data$expired_age[1]
    sample.list <- c(sample.list, sample_name)
    print(paste(sample_name, expired_age))
}
combined <- merge(data.list[[1]], y = data.list[2:length(data.list)], add.cell.ids = sample.list)
combined@meta.data$orig.ident <- as.factor(combined@meta.data$orig.ident)
combined@meta.data$sample_name <- as.factor(combined@meta.data$sample_name)
combined@meta.data$batch <- as.factor(combined@meta.data$batch)
combined@meta.data$PMI_num <- combined@meta.data$PMI
combined@meta.data$PMI <- droplevels(cut(combined@meta.data$PMI, seq(0,10,by=1)))
combined@meta.data$RIN_num <- combined@meta.data$RIN
combined@meta.data$RIN <- droplevels(cut(combined@meta.data$RIN, seq(0,11,by=1)))
combined@meta.data$sex <- as.factor(combined@meta.data$sex)
age_batch_size = 10
combined@meta.data$expired_age_num <- combined@meta.data$expired_age
combined@meta.data$expired_age <- droplevels(cut(combined@meta.data$expired_age, seq(0,130,by=age_batch_size)))
combined@meta.data$diagnosis <- as.factor(combined@meta.data$diagnosis)
assert(is.numeric(combined@meta.data$PMI_num))
assert(is.numeric(combined@meta.data$RIN_num))
assert(is.numeric(combined@meta.data$expired_age_num))
assert(is.numeric(combined@meta.data$last_mmse_test_score))
assert(is.numeric(combined@meta.data$motor_updrs_score))
assert(is.numeric(combined@meta.data$nctx_temporal))
assert(is.numeric(combined@meta.data$brain_stem_sn))
saveRDS(combined, "comb_seurat_postqc.rds")
