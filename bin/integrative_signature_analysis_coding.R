#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(ggpubr))

args <- commandArgs()

working_dir <- args[6]
sbs_cosmic2 <- args[7]
sbs_cosmic3_2 <- args [8]

#1) Load cosine similarity values

path_folder_matrix1 <- list.files(paste(working_dir, "Cosine_Similarity", "matrix1", sep="/"), full.names=TRUE)
name_files_matrix1 <- list.files(paste(working_dir, "Cosine_Similarity", "matrix1", sep="/"))

cos_sim_matrix1 <- lapply(path_folder_matrix1, function(x) read.table(file=x, row.names=1))
names(cos_sim_matrix1) <- str_sub(name_files_matrix1, start = 16, end = -5)
cos_sim_matrix1 <- lapply(cos_sim_matrix1, function(x) {colnames(x) <- "cos_sim"; x})

path_folder_matrix2 <- list.files(paste(working_dir, "Cosine_Similarity", "matrix2", sep="/"), full.names=TRUE)
name_files_matrix2 <- list.files(paste(working_dir, "Cosine_Similarity", "matrix2", sep="/"))

cos_sim_matrix2 <- lapply(path_folder_matrix2, function(x) read.table(file=x, row.names=1))
names(cos_sim_matrix2) <- str_sub(name_files_matrix2, start = 16, end = -5)
cos_sim_matrix2 <- lapply(cos_sim_matrix2, function(x) {colnames(x) <- "cos_sim"; x})

#coding mutations only

path_folder_matrix1_coding <- list.files(paste(working_dir, "Cosine_Similarity", "matrix1_coding", sep="/"), full.names=TRUE)
name_files_matrix1_coding <- list.files(paste(working_dir, "Cosine_Similarity", "matrix1_coding", sep="/"))

cos_sim_matrix1_coding <- lapply(path_folder_matrix1_coding, function(x) read.table(file=x, row.names=1))
names(cos_sim_matrix1_coding) <- str_sub(name_files_matrix1_coding, start = 16, end = -5)
cos_sim_matrix1_coding <- lapply(cos_sim_matrix1_coding, function(x) {colnames(x) <- "cos_sim"; x})

path_folder_matrix2_coding <- list.files(paste(working_dir, "Cosine_Similarity", "matrix2_coding", sep="/"), full.names=TRUE)
name_files_matrix2_coding <- list.files(paste(working_dir, "Cosine_Similarity", "matrix2_coding", sep="/"))

cos_sim_matrix2_coding <- lapply(path_folder_matrix2_coding, function(x) read.table(file=x, row.names=1))
names(cos_sim_matrix2_coding) <- str_sub(name_files_matrix2_coding, start = 16, end = -5)
cos_sim_matrix2_coding <- lapply(cos_sim_matrix2_coding, function(x) {colnames(x) <- "cos_sim"; x})


#non-coding mutations only

path_folder_matrix1_non_coding <- list.files(paste(working_dir, "Cosine_Similarity", "matrix1_non_coding", sep="/"), full.names=TRUE)
name_files_matrix1_non_coding <- list.files(paste(working_dir, "Cosine_Similarity", "matrix1_non_coding", sep="/"))

cos_sim_matrix1_non_coding <- lapply(path_folder_matrix1_non_coding, function(x) read.table(file=x, row.names=1))
names(cos_sim_matrix1_non_coding) <- str_sub(name_files_matrix1_non_coding, start = 16, end = -5)
cos_sim_matrix1_non_coding <- lapply(cos_sim_matrix1_non_coding, function(x) {colnames(x) <- "cos_sim"; x})

path_folder_matrix2_non_coding <- list.files(paste(working_dir, "Cosine_Similarity", "matrix2_non_coding", sep="/"), full.names=TRUE)
name_files_matrix2_non_coding <- list.files(paste(working_dir, "Cosine_Similarity", "matrix2_non_coding", sep="/"))

cos_sim_matrix2_non_coding <- lapply(path_folder_matrix2_non_coding, function(x) read.table(file=x, row.names=1))
names(cos_sim_matrix2_non_coding) <- str_sub(name_files_matrix2_non_coding, start = 16, end = -5)
cos_sim_matrix2_non_coding <- lapply(cos_sim_matrix2_non_coding, function(x) {colnames(x) <- "cos_sim"; x})


# Load SBS signatures values

path_folder_matrix1 <- list.files(paste(working_dir, "SBS_signature_contributions", "matrix1", sep="/"), full.names=TRUE)
name_files_matrix1 <- list.files(paste(working_dir, "SBS_signature_contributions", "matrix1", sep="/"))

sbs_matrix1 <- lapply(path_folder_matrix1, function(x) read.table(file=x, row.names=1, header=T))
names(sbs_matrix1) <- str_sub(name_files_matrix1, start = 9, end = -5)

path_folder_matrix2 <- list.files(paste(working_dir, "SBS_signature_contributions", "matrix2", sep="/"), full.names=TRUE)
name_files_matrix2 <- list.files(paste(working_dir, "SBS_signature_contributions", "matrix2", sep="/"))

sbs_matrix2 <- lapply(path_folder_matrix2, function(x) read.table(file=x, row.names=1, header=T))
names(sbs_matrix2) <- str_sub(name_files_matrix2, start = 9, end = -5) 	

#coding mutations only

path_folder_matrix1_coding <- list.files(paste(working_dir, "SBS_signature_contributions", "matrix1_coding", sep="/"), full.names=TRUE)
name_files_matrix1_coding <- list.files(paste(working_dir, "SBS_signature_contributions", "matrix1_coding", sep="/"))

sbs_matrix1_coding <- lapply(path_folder_matrix1_coding, function(x) read.table(file=x, row.names=1, header=T))
names(sbs_matrix1_coding) <- str_sub(name_files_matrix1_coding, start = 9, end = -5)

path_folder_matrix2_coding <- list.files(paste(working_dir, "SBS_signature_contributions", "matrix2_coding", sep="/"), full.names=TRUE)
name_files_matrix2_coding <- list.files(paste(working_dir, "SBS_signature_contributions", "matrix2_coding", sep="/"))

sbs_matrix2_coding <- lapply(path_folder_matrix2_coding, function(x) read.table(file=x, row.names=1, header=T))
names(sbs_matrix2_coding) <- str_sub(name_files_matrix2_coding, start = 9, end = -5)

#non-coding mutations only

path_folder_matrix1_non_coding <- list.files(paste(working_dir, "SBS_signature_contributions", "matrix1_non_coding", sep="/"), full.names=TRUE)
name_files_matrix1_non_coding <- list.files(paste(working_dir, "SBS_signature_contributions", "matrix1_non_coding", sep="/"))

sbs_matrix1_non_coding <- lapply(path_folder_matrix1_non_coding, function(x) read.table(file=x, row.names=1, header=T))
names(sbs_matrix1_non_coding) <- str_sub(name_files_matrix1_non_coding, start = 9, end = -5)

path_folder_matrix2_non_coding <- list.files(paste(working_dir, "SBS_signature_contributions", "matrix2_non_coding", sep="/"), full.names=TRUE)
name_files_matrix2_non_coding <- list.files(paste(working_dir, "SBS_signature_contributions", "matrix2_non_coding", sep="/"))

sbs_matrix2_non_coding <- lapply(path_folder_matrix2_non_coding, function(x) read.table(file=x, row.names=1, header=T))
names(sbs_matrix2_non_coding) <- str_sub(name_files_matrix2_non_coding, start = 9, end = -5)

#2) load SBS to evaluate

sbs_cosmic2 <- unlist(strsplit(sbs_cosmic2, ","))
sbs_cosmic3_2 <- unlist(strsplit(sbs_cosmic3_2, ","))

#3) calc contribution

C2_sbs_cosmic2_MP_m1 <- (apply(sbs_matrix1$C2_MutationalPatterns[sbs_cosmic2,], 2, sum))
C2_sbs_cosmic2_MP_m1_df <- data.frame(C2_sbs_cosmic2_MP_m1, "MutationalPatterns", "matrix1")
colnames(C2_sbs_cosmic2_MP_m1_df) <- c("SBS_contr", "Tool", "Condition")

C2_sbs_cosmic2_MP_m2 <- (apply(sbs_matrix2$C2_MutationalPatterns[sbs_cosmic2,], 2, sum))
C2_sbs_cosmic2_MP_m2_df <- data.frame(C2_sbs_cosmic2_MP_m2, "MutationalPatterns", "matrix2")
colnames(C2_sbs_cosmic2_MP_m2_df) <- c("SBS_contr", "Tool", "Condition")

C2_sbs_cosmic2_DS_m1 <- (apply(sbs_matrix1$C2_deconstructsigs[sbs_cosmic2,], 2, sum))
C2_sbs_cosmic2_DS_m1_df <- data.frame(C2_sbs_cosmic2_DS_m1, "deconstructsigs", "matrix1")
colnames(C2_sbs_cosmic2_DS_m1_df) <- c("SBS_contr", "Tool", "Condition")

C2_sbs_cosmic2_DS_m2 <- (apply(sbs_matrix2$C2_deconstructsigs[sbs_cosmic2,], 2, sum))
C2_sbs_cosmic2_DS_m2_df <- data.frame(C2_sbs_cosmic2_DS_m2, "deconstructsigs", "matrix2")
colnames(C2_sbs_cosmic2_DS_m2_df) <- c("SBS_contr", "Tool", "Condition")

C2_sbs_cosmic2_STL_m1 <- (apply(sbs_matrix1$C2_signature.tools.lib[sbs_cosmic2,], 2, sum))
C2_sbs_cosmic2_STL_m1_df <- data.frame(C2_sbs_cosmic2_STL_m1, "signature.tools.lib", "matrix1")
colnames(C2_sbs_cosmic2_STL_m1_df) <- c("SBS_contr", "Tool", "Condition")

C2_sbs_cosmic2_STL_m2 <- (apply(sbs_matrix2$C2_signature.tools.lib[sbs_cosmic2,], 2, sum))
C2_sbs_cosmic2_STL_m2_df <- data.frame(C2_sbs_cosmic2_STL_m2, "signature.tools.lib", "matrix2")
colnames(C2_sbs_cosmic2_STL_m2_df) <- c("SBS_contr", "Tool", "Condition")

C2_sbs_cosmic2_SA_m1 <- (apply(sbs_matrix1$C2_sigprofilerassignment[sbs_cosmic2,], 2, sum))
C2_sbs_cosmic2_SA_m1_df <- data.frame(C2_sbs_cosmic2_SA_m1, "sigprofilerassignment", "matrix1")
colnames(C2_sbs_cosmic2_SA_m1_df) <- c("SBS_contr", "Tool", "Condition")

C2_sbs_cosmic2_SA_m2 <- (apply(sbs_matrix2$C2_sigprofilerassignment[sbs_cosmic2,], 2, sum))
C2_sbs_cosmic2_SA_m2_df <- data.frame(C2_sbs_cosmic2_SA_m2, "sigprofilerassignment", "matrix2")
colnames(C2_sbs_cosmic2_SA_m2_df) <- c("SBS_contr", "Tool", "Condition")

SBS_C2 <- rbind(C2_sbs_cosmic2_MP_m1_df, C2_sbs_cosmic2_MP_m2_df, C2_sbs_cosmic2_DS_m1_df, C2_sbs_cosmic2_DS_m2_df, C2_sbs_cosmic2_STL_m1_df, C2_sbs_cosmic2_STL_m2_df, C2_sbs_cosmic2_SA_m1_df, C2_sbs_cosmic2_SA_m2_df) ### TO PLOT ###


#coding mutations only


C2_sbs_cosmic2_MP_m1_coding <- (apply(sbs_matrix1_coding$C2_MutationalPatterns[sbs_cosmic2,], 2, sum))
C2_sbs_cosmic2_MP_m1_coding_df <- data.frame(C2_sbs_cosmic2_MP_m1_coding, "MutationalPatterns", "matrix1")
colnames(C2_sbs_cosmic2_MP_m1_coding_df) <- c("SBS_contr", "Tool", "Condition")

C2_sbs_cosmic2_MP_m2_coding <- (apply(sbs_matrix2_coding$C2_MutationalPatterns[sbs_cosmic2,], 2, sum))
C2_sbs_cosmic2_MP_m2_coding_df <- data.frame(C2_sbs_cosmic2_MP_m2_coding, "MutationalPatterns", "matrix2")
colnames(C2_sbs_cosmic2_MP_m2_coding_df) <- c("SBS_contr", "Tool", "Condition")

C2_sbs_cosmic2_DS_m1_coding <- (apply(sbs_matrix1_coding$C2_deconstructsigs[sbs_cosmic2,], 2, sum))
C2_sbs_cosmic2_DS_m1_coding_df <- data.frame(C2_sbs_cosmic2_DS_m1_coding, "deconstructsigs", "matrix1")
colnames(C2_sbs_cosmic2_DS_m1_coding_df) <- c("SBS_contr", "Tool", "Condition")

C2_sbs_cosmic2_DS_m2_coding <- (apply(sbs_matrix2_coding$C2_deconstructsigs[sbs_cosmic2,], 2, sum))
C2_sbs_cosmic2_DS_m2_coding_df <- data.frame(C2_sbs_cosmic2_DS_m2_coding, "deconstructsigs", "matrix2")
colnames(C2_sbs_cosmic2_DS_m2_coding_df) <- c("SBS_contr", "Tool", "Condition")

C2_sbs_cosmic2_STL_m1_coding <- (apply(sbs_matrix1_coding$C2_signature.tools.lib[sbs_cosmic2,], 2, sum))
C2_sbs_cosmic2_STL_m1_coding_df <- data.frame(C2_sbs_cosmic2_STL_m1_coding, "signature.tools.lib", "matrix1")
colnames(C2_sbs_cosmic2_STL_m1_coding_df) <- c("SBS_contr", "Tool", "Condition")

C2_sbs_cosmic2_STL_m2_coding <- (apply(sbs_matrix2_coding$C2_signature.tools.lib[sbs_cosmic2,], 2, sum))
C2_sbs_cosmic2_STL_m2_coding_df <- data.frame(C2_sbs_cosmic2_STL_m2_coding, "signature.tools.lib", "matrix2")
colnames(C2_sbs_cosmic2_STL_m2_coding_df) <- c("SBS_contr", "Tool", "Condition")

C2_sbs_cosmic2_SA_m1_coding <- (apply(sbs_matrix1_coding$C2_sigprofilerassignment[sbs_cosmic2,], 2, sum))
C2_sbs_cosmic2_SA_m1_coding_df <- data.frame(C2_sbs_cosmic2_SA_m1_coding, "sigprofilerassignment", "matrix1")
colnames(C2_sbs_cosmic2_SA_m1_coding_df) <- c("SBS_contr", "Tool", "Condition")

C2_sbs_cosmic2_SA_m2_coding <- (apply(sbs_matrix2_coding$C2_sigprofilerassignment[sbs_cosmic2,], 2, sum))
C2_sbs_cosmic2_SA_m2_coding_df <- data.frame(C2_sbs_cosmic2_SA_m2_coding, "sigprofilerassignment", "matrix2")
colnames(C2_sbs_cosmic2_SA_m2_coding_df) <- c("SBS_contr", "Tool", "Condition")

SBS_C2_coding <- rbind(C2_sbs_cosmic2_MP_m1_coding_df, C2_sbs_cosmic2_MP_m2_coding_df, C2_sbs_cosmic2_DS_m1_coding_df, C2_sbs_cosmic2_DS_m2_coding_df, C2_sbs_cosmic2_STL_m1_coding_df, C2_sbs_cosmic2_STL_m2_coding_df, C2_sbs_cosmic2_SA_m1_coding_df, C2_sbs_cosmic2_SA_m2_coding_df) ### TO PLOT ###


#non-coding mutations only

C2_sbs_cosmic2_MP_m1_non_coding <- (apply(sbs_matrix1_non_coding$C2_MutationalPatterns[sbs_cosmic2,], 2, sum))
C2_sbs_cosmic2_MP_m1_non_coding_df <- data.frame(C2_sbs_cosmic2_MP_m1_non_coding, "MutationalPatterns", "matrix1")
colnames(C2_sbs_cosmic2_MP_m1_non_coding_df) <- c("SBS_contr", "Tool", "Condition")

C2_sbs_cosmic2_MP_m2_non_coding <- (apply(sbs_matrix2_non_coding$C2_MutationalPatterns[sbs_cosmic2,], 2, sum))
C2_sbs_cosmic2_MP_m2_non_coding_df <- data.frame(C2_sbs_cosmic2_MP_m2_non_coding, "MutationalPatterns", "matrix2")
colnames(C2_sbs_cosmic2_MP_m2_non_coding_df) <- c("SBS_contr", "Tool", "Condition")

C2_sbs_cosmic2_DS_m1_non_coding <- (apply(sbs_matrix1_non_coding$C2_deconstructsigs[sbs_cosmic2,], 2, sum))
C2_sbs_cosmic2_DS_m1_non_coding_df <- data.frame(C2_sbs_cosmic2_DS_m1_non_coding, "deconstructsigs", "matrix1")
colnames(C2_sbs_cosmic2_DS_m1_non_coding_df) <- c("SBS_contr", "Tool", "Condition")

C2_sbs_cosmic2_DS_m2_non_coding <- (apply(sbs_matrix2_non_coding$C2_deconstructsigs[sbs_cosmic2,], 2, sum))
C2_sbs_cosmic2_DS_m2_non_coding_df <- data.frame(C2_sbs_cosmic2_DS_m2_non_coding, "deconstructsigs", "matrix2")
colnames(C2_sbs_cosmic2_DS_m2_non_coding_df) <- c("SBS_contr", "Tool", "Condition")

C2_sbs_cosmic2_STL_m1_non_coding <- (apply(sbs_matrix1_non_coding$C2_signature.tools.lib[sbs_cosmic2,], 2, sum))
C2_sbs_cosmic2_STL_m1_non_coding_df <- data.frame(C2_sbs_cosmic2_STL_m1_non_coding, "signature.tools.lib", "matrix1")
colnames(C2_sbs_cosmic2_STL_m1_non_coding_df) <- c("SBS_contr", "Tool", "Condition")

C2_sbs_cosmic2_STL_m2_non_coding <- (apply(sbs_matrix2_non_coding$C2_signature.tools.lib[sbs_cosmic2,], 2, sum))
C2_sbs_cosmic2_STL_m2_non_coding_df <- data.frame(C2_sbs_cosmic2_STL_m2_non_coding, "signature.tools.lib", "matrix2")
colnames(C2_sbs_cosmic2_STL_m2_non_coding_df) <- c("SBS_contr", "Tool", "Condition")

C2_sbs_cosmic2_SA_m1_non_coding <- (apply(sbs_matrix1_non_coding$C2_sigprofilerassignment[sbs_cosmic2,], 2, sum))
C2_sbs_cosmic2_SA_m1_non_coding_df <- data.frame(C2_sbs_cosmic2_SA_m1_non_coding, "sigprofilerassignment", "matrix1")
colnames(C2_sbs_cosmic2_SA_m1_non_coding_df) <- c("SBS_contr", "Tool", "Condition")

C2_sbs_cosmic2_SA_m2_non_coding <- (apply(sbs_matrix2_non_coding$C2_sigprofilerassignment[sbs_cosmic2,], 2, sum))
C2_sbs_cosmic2_SA_m2_non_coding_df <- data.frame(C2_sbs_cosmic2_SA_m2_non_coding, "sigprofilerassignment", "matrix2")
colnames(C2_sbs_cosmic2_SA_m2_non_coding_df) <- c("SBS_contr", "Tool", "Condition")

SBS_C2_non_coding <- rbind(C2_sbs_cosmic2_MP_m1_non_coding_df, C2_sbs_cosmic2_MP_m2_non_coding_df, C2_sbs_cosmic2_DS_m1_non_coding_df, C2_sbs_cosmic2_DS_m2_non_coding_df, C2_sbs_cosmic2_STL_m1_non_coding_df, C2_sbs_cosmic2_STL_m2_non_coding_df, C2_sbs_cosmic2_SA_m1_non_coding_df, C2_sbs_cosmic2_SA_m2_non_coding_df) ### TO PLOT ###


#SBS contribution cosmic v3.2

C3_2_sbs_cosmic3_2_MP_m1 <- (apply(sbs_matrix1$C3.2_MutationalPatterns[sbs_cosmic3_2,], 2, sum))
C3_2_sbs_cosmic3_2_MP_m1_df <- data.frame(C3_2_sbs_cosmic3_2_MP_m1, "MutationalPatterns", "matrix1")
colnames(C3_2_sbs_cosmic3_2_MP_m1_df) <- c("SBS_contr", "Tool", "Condition")

C3_2_sbs_cosmic3_2_MP_m2 <- (apply(sbs_matrix2$C3.2_MutationalPatterns[sbs_cosmic3_2,], 2, sum))
C3_2_sbs_cosmic3_2_MP_m2_df <- data.frame(C3_2_sbs_cosmic3_2_MP_m2, "MutationalPatterns", "matrix2")
colnames(C3_2_sbs_cosmic3_2_MP_m2_df) <- c("SBS_contr", "Tool", "Condition")

C3_2_sbs_cosmic3_2_DS_m1 <- (apply(sbs_matrix1$C3.2_deconstructsigs[sbs_cosmic3_2,], 2, sum))
C3_2_sbs_cosmic3_2_DS_m1_df <- data.frame(C3_2_sbs_cosmic3_2_DS_m1, "deconstructsigs", "matrix1")
colnames(C3_2_sbs_cosmic3_2_DS_m1_df) <- c("SBS_contr", "Tool", "Condition")

C3_2_sbs_cosmic3_2_DS_m2 <- (apply(sbs_matrix2$C3.2_deconstructsigs[sbs_cosmic3_2,], 2, sum))
C3_2_sbs_cosmic3_2_DS_m2_df <- data.frame(C3_2_sbs_cosmic3_2_DS_m2, "deconstructsigs", "matrix2")
colnames(C3_2_sbs_cosmic3_2_DS_m2_df) <- c("SBS_contr", "Tool", "Condition")

C3_2_sbs_cosmic3_2_STL_m1 <- (apply(sbs_matrix1$C3.2_signature.tools.lib[sbs_cosmic3_2,], 2, sum))
C3_2_sbs_cosmic3_2_STL_m1_df <- data.frame(C3_2_sbs_cosmic3_2_STL_m1, "signature.tools.lib", "matrix1")
colnames(C3_2_sbs_cosmic3_2_STL_m1_df) <- c("SBS_contr", "Tool", "Condition")

C3_2_sbs_cosmic3_2_STL_m2 <- (apply(sbs_matrix2$C3.2_signature.tools.lib[sbs_cosmic3_2,], 2, sum))
C3_2_sbs_cosmic3_2_STL_m2_df <- data.frame(C3_2_sbs_cosmic3_2_STL_m2, "signature.tools.lib", "matrix2")
colnames(C3_2_sbs_cosmic3_2_STL_m2_df) <- c("SBS_contr", "Tool", "Condition")

C3_2_sbs_cosmic3_2_SA_m1 <- (apply(sbs_matrix1$C3.2_sigprofilerassignment[sbs_cosmic3_2,], 2, sum))
C3_2_sbs_cosmic3_2_SA_m1_df <- data.frame(C3_2_sbs_cosmic3_2_SA_m1, "sigprofilerassignment", "matrix1")
colnames(C3_2_sbs_cosmic3_2_SA_m1_df) <- c("SBS_contr", "Tool", "Condition")

C3_2_sbs_cosmic3_2_SA_m2 <- (apply(sbs_matrix2$C3.2_sigprofilerassignment[sbs_cosmic3_2,], 2, sum))
C3_2_sbs_cosmic3_2_SA_m2_df <- data.frame(C3_2_sbs_cosmic3_2_SA_m2, "sigprofilerassignment", "matrix2")
colnames(C3_2_sbs_cosmic3_2_SA_m2_df) <- c("SBS_contr", "Tool", "Condition")

SBS_C3_2 <- rbind(C3_2_sbs_cosmic3_2_MP_m1_df, C3_2_sbs_cosmic3_2_MP_m2_df, C3_2_sbs_cosmic3_2_DS_m1_df, C3_2_sbs_cosmic3_2_DS_m2_df, C3_2_sbs_cosmic3_2_STL_m1_df, C3_2_sbs_cosmic3_2_STL_m2_df, C3_2_sbs_cosmic3_2_SA_m1_df, C3_2_sbs_cosmic3_2_SA_m2_df) ### TO PLOT ###


#coding mutations only

C3_2_sbs_cosmic3_2_MP_m1_coding <- (apply(sbs_matrix1_coding$C3.2_MutationalPatterns[sbs_cosmic3_2,], 2, sum))
C3_2_sbs_cosmic3_2_MP_m1_coding_df <- data.frame(C3_2_sbs_cosmic3_2_MP_m1_coding, "MutationalPatterns", "matrix1")
colnames(C3_2_sbs_cosmic3_2_MP_m1_coding_df) <- c("SBS_contr", "Tool", "Condition")

C3_2_sbs_cosmic3_2_MP_m2_coding <- (apply(sbs_matrix2_coding$C3.2_MutationalPatterns[sbs_cosmic3_2,], 2, sum))
C3_2_sbs_cosmic3_2_MP_m2_coding_df <- data.frame(C3_2_sbs_cosmic3_2_MP_m2_coding, "MutationalPatterns", "matrix2")
colnames(C3_2_sbs_cosmic3_2_MP_m2_coding_df) <- c("SBS_contr", "Tool", "Condition")

C3_2_sbs_cosmic3_2_DS_m1_coding <- (apply(sbs_matrix1_coding$C3.2_deconstructsigs[sbs_cosmic3_2,], 2, sum))
C3_2_sbs_cosmic3_2_DS_m1_coding_df <- data.frame(C3_2_sbs_cosmic3_2_DS_m1_coding, "deconstructsigs", "matrix1")
colnames(C3_2_sbs_cosmic3_2_DS_m1_coding_df) <- c("SBS_contr", "Tool", "Condition")

C3_2_sbs_cosmic3_2_DS_m2_coding <- (apply(sbs_matrix2_coding$C3.2_deconstructsigs[sbs_cosmic3_2,], 2, sum))
C3_2_sbs_cosmic3_2_DS_m2_coding_df <- data.frame(C3_2_sbs_cosmic3_2_DS_m2_coding, "deconstructsigs", "matrix2")
colnames(C3_2_sbs_cosmic3_2_DS_m2_coding_df) <- c("SBS_contr", "Tool", "Condition")

C3_2_sbs_cosmic3_2_STL_m1_coding <- (apply(sbs_matrix1_coding$C3.2_signature.tools.lib[sbs_cosmic3_2,], 2, sum))
C3_2_sbs_cosmic3_2_STL_m1_coding_df <- data.frame(C3_2_sbs_cosmic3_2_STL_m1_coding, "signature.tools.lib", "matrix1")
colnames(C3_2_sbs_cosmic3_2_STL_m1_coding_df) <- c("SBS_contr", "Tool", "Condition")

C3_2_sbs_cosmic3_2_STL_m2_coding <- (apply(sbs_matrix2_coding$C3.2_signature.tools.lib[sbs_cosmic3_2,], 2, sum))
C3_2_sbs_cosmic3_2_STL_m2_coding_df <- data.frame(C3_2_sbs_cosmic3_2_STL_m2_coding, "signature.tools.lib", "matrix2")
colnames(C3_2_sbs_cosmic3_2_STL_m2_coding_df) <- c("SBS_contr", "Tool", "Condition")

C3_2_sbs_cosmic3_2_SA_m1_coding <- (apply(sbs_matrix1_coding$C3.2_sigprofilerassignment[sbs_cosmic3_2,], 2, sum))
C3_2_sbs_cosmic3_2_SA_m1_coding_df <- data.frame(C3_2_sbs_cosmic3_2_SA_m1_coding, "sigprofilerassignment", "matrix1")
colnames(C3_2_sbs_cosmic3_2_SA_m1_coding_df) <- c("SBS_contr", "Tool", "Condition")

C3_2_sbs_cosmic3_2_SA_m2_coding <- (apply(sbs_matrix2_coding$C3.2_sigprofilerassignment[sbs_cosmic3_2,], 2, sum))
C3_2_sbs_cosmic3_2_SA_m2_coding_df <- data.frame(C3_2_sbs_cosmic3_2_SA_m2_coding, "sigprofilerassignment", "matrix2")
colnames(C3_2_sbs_cosmic3_2_SA_m2_coding_df) <- c("SBS_contr", "Tool", "Condition")

SBS_C3_2_coding <- rbind(C3_2_sbs_cosmic3_2_MP_m1_coding_df, C3_2_sbs_cosmic3_2_MP_m2_coding_df, C3_2_sbs_cosmic3_2_DS_m1_coding_df, C3_2_sbs_cosmic3_2_DS_m2_coding_df, C3_2_sbs_cosmic3_2_STL_m1_coding_df, C3_2_sbs_cosmic3_2_STL_m2_coding_df, C3_2_sbs_cosmic3_2_SA_m1_coding_df, C3_2_sbs_cosmic3_2_SA_m2_coding_df) ### TO PLOT ###


#non-coding mutations only

C3_2_sbs_cosmic3_2_MP_m1_non_coding <- (apply(sbs_matrix1_non_coding$C3.2_MutationalPatterns[sbs_cosmic3_2,], 2, sum))
C3_2_sbs_cosmic3_2_MP_m1_non_coding_df <- data.frame(C3_2_sbs_cosmic3_2_MP_m1_non_coding, "MutationalPatterns", "matrix1")
colnames(C3_2_sbs_cosmic3_2_MP_m1_non_coding_df) <- c("SBS_contr", "Tool", "Condition")

C3_2_sbs_cosmic3_2_MP_m2_non_coding <- (apply(sbs_matrix2_non_coding$C3.2_MutationalPatterns[sbs_cosmic3_2,], 2, sum))
C3_2_sbs_cosmic3_2_MP_m2_non_coding_df <- data.frame(C3_2_sbs_cosmic3_2_MP_m2_non_coding, "MutationalPatterns", "matrix2")
colnames(C3_2_sbs_cosmic3_2_MP_m2_non_coding_df) <- c("SBS_contr", "Tool", "Condition")
 
C3_2_sbs_cosmic3_2_DS_m1_non_coding <- (apply(sbs_matrix1_non_coding$C3.2_deconstructsigs[sbs_cosmic3_2,], 2, sum))
C3_2_sbs_cosmic3_2_DS_m1_non_coding_df <- data.frame(C3_2_sbs_cosmic3_2_DS_m1_non_coding, "deconstructsigs", "matrix1")
colnames(C3_2_sbs_cosmic3_2_DS_m1_non_coding_df) <- c("SBS_contr", "Tool", "Condition")

C3_2_sbs_cosmic3_2_DS_m2_non_coding <- (apply(sbs_matrix2_non_coding$C3.2_deconstructsigs[sbs_cosmic3_2,], 2, sum))
C3_2_sbs_cosmic3_2_DS_m2_non_coding_df <- data.frame(C3_2_sbs_cosmic3_2_DS_m2_non_coding, "deconstructsigs", "matrix2")
colnames(C3_2_sbs_cosmic3_2_DS_m2_non_coding_df) <- c("SBS_contr", "Tool", "Condition")

C3_2_sbs_cosmic3_2_STL_m1_non_coding <- (apply(sbs_matrix1_non_coding$C3.2_signature.tools.lib[sbs_cosmic3_2,], 2, sum))
C3_2_sbs_cosmic3_2_STL_m1_non_coding_df <- data.frame(C3_2_sbs_cosmic3_2_STL_m1_non_coding, "signature.tools.lib", "matrix1")
colnames(C3_2_sbs_cosmic3_2_STL_m1_non_coding_df) <- c("SBS_contr", "Tool", "Condition")

C3_2_sbs_cosmic3_2_STL_m2_non_coding <- (apply(sbs_matrix2_non_coding$C3.2_signature.tools.lib[sbs_cosmic3_2,], 2, sum))
C3_2_sbs_cosmic3_2_STL_m2_non_coding_df <- data.frame(C3_2_sbs_cosmic3_2_STL_m2_non_coding, "signature.tools.lib", "matrix2")
colnames(C3_2_sbs_cosmic3_2_STL_m2_non_coding_df) <- c("SBS_contr", "Tool", "Condition")

C3_2_sbs_cosmic3_2_SA_m1_non_coding <- (apply(sbs_matrix1_non_coding$C3.2_sigprofilerassignment[sbs_cosmic3_2,], 2, sum))
C3_2_sbs_cosmic3_2_SA_m1_non_coding_df <- data.frame(C3_2_sbs_cosmic3_2_SA_m1_non_coding, "sigprofilerassignment", "matrix1")
colnames(C3_2_sbs_cosmic3_2_SA_m1_non_coding_df) <- c("SBS_contr", "Tool", "Condition")

C3_2_sbs_cosmic3_2_SA_m2_non_coding <- (apply(sbs_matrix2_non_coding$C3.2_sigprofilerassignment[sbs_cosmic3_2,], 2, sum))
C3_2_sbs_cosmic3_2_SA_m2_non_coding_df <- data.frame(C3_2_sbs_cosmic3_2_SA_m2_non_coding, "sigprofilerassignment", "matrix2")
colnames(C3_2_sbs_cosmic3_2_SA_m2_non_coding_df) <- c("SBS_contr", "Tool", "Condition")

SBS_C3_2_non_coding <- rbind(C3_2_sbs_cosmic3_2_MP_m1_non_coding_df, C3_2_sbs_cosmic3_2_MP_m2_non_coding_df, C3_2_sbs_cosmic3_2_DS_m1_non_coding_df, C3_2_sbs_cosmic3_2_DS_m2_non_coding_df, C3_2_sbs_cosmic3_2_STL_m1_non_coding_df, C3_2_sbs_cosmic3_2_STL_m2_non_coding_df, C3_2_sbs_cosmic3_2_SA_m1_non_coding_df, C3_2_sbs_cosmic3_2_SA_m2_non_coding_df) ### TO PLOT ###


#COSINE SIMILARITY
#Cosmic2
C2_MP_m1 <- data.frame(cos_sim_matrix1$C2_MutationalPatterns, "MutationalPatterns")
colnames(C2_MP_m1) <- c("cosine_similarity", "Tool")

C2_DS_m1 <- data.frame(cos_sim_matrix1$C2_deconstructsigs, "deconstructsig")
colnames(C2_DS_m1) <- c("cosine_similarity", "Tool")

C2_STL_m1 <- data.frame(cos_sim_matrix1$C2_signature.tools.lib, "signature.tools.lib")
colnames(C2_STL_m1) <- c("cosine_similarity", "Tool")

C2_SA_m1 <- data.frame(cos_sim_matrix1$C2_sigprofilerassignment, "sigprofilerassignment")
colnames(C2_SA_m1) <- c("cosine_similarity", "Tool")

cos_sim_cosmic2_matrix1 <- rbind(C2_MP_m1,C2_DS_m1,C2_STL_m1,C2_SA_m1) ### TO PLOT ###

C2_MP_m2 <- data.frame(cos_sim_matrix2$C2_MutationalPatterns, "MutationalPatterns")
colnames(C2_MP_m2) <- c("cosine_similarity", "Tool")

C2_DS_m2 <- data.frame(cos_sim_matrix2$C2_deconstructsigs, "deconstructsig")
colnames(C2_DS_m2) <- c("cosine_similarity", "Tool")

C2_STL_m2 <- data.frame(cos_sim_matrix2$C2_signature.tools.lib, "signature.tools.lib")
colnames(C2_STL_m2) <- c("cosine_similarity", "Tool")

C2_SA_m2 <- data.frame(cos_sim_matrix2$C2_sigprofilerassignment, "sigprofilerassignment")
colnames(C2_SA_m2) <- c("cosine_similarity", "Tool")

cos_sim_cosmic2_matrix2 <- rbind(C2_MP_m2,C2_DS_m2,C2_STL_m2,C2_SA_m2) ### TO PLOT ###

#Cosmic3.2

C3.2_MP_m1 <- data.frame(cos_sim_matrix1$C3.2_MutationalPatterns, "MutationalPatterns")
colnames(C3.2_MP_m1) <- c("cosine_similarity", "Tool")

C3.2_DS_m1 <- data.frame(cos_sim_matrix1$C3.2_deconstructsigs, "deconstructsig")
colnames(C3.2_DS_m1) <- c("cosine_similarity", "Tool")

C3.2_STL_m1 <- data.frame(cos_sim_matrix1$C3.2_signature.tools.lib, "signature.tools.lib")
colnames(C3.2_STL_m1) <- c("cosine_similarity", "Tool")

C3.2_SA_m1 <- data.frame(cos_sim_matrix1$C3.2_sigprofilerassignment, "sigprofilerassignment")
colnames(C3.2_SA_m1) <- c("cosine_similarity", "Tool")

cos_sim_cosmic3_2_matrix1 <- rbind(C3.2_MP_m1,C3.2_DS_m1,C3.2_STL_m1,C3.2_SA_m1) ### TO PLOT ###

C3.2_MP_m2 <- data.frame(cos_sim_matrix2$C3.2_MutationalPatterns, "MutationalPatterns")
colnames(C3.2_MP_m2) <- c("cosine_similarity", "Tool")

C3.2_DS_m2 <- data.frame(cos_sim_matrix2$C3.2_deconstructsigs, "deconstructsig")
colnames(C3.2_DS_m2) <- c("cosine_similarity", "Tool")

C3.2_STL_m2 <- data.frame(cos_sim_matrix2$C3.2_signature.tools.lib, "signature.tools.lib")
colnames(C3.2_STL_m2) <- c("cosine_similarity", "Tool")

C3.2_SA_m2 <- data.frame(cos_sim_matrix2$C3.2_sigprofilerassignment, "sigprofilerassignment")
colnames(C3.2_SA_m2) <- c("cosine_similarity", "Tool")

cos_sim_cosmic3_2_matrix2 <- rbind(C3.2_MP_m2,C3.2_DS_m2,C3.2_STL_m2,C3.2_SA_m2) ### TO PLOT ###

#Coding mutations only

C2_MP_m1_coding <- data.frame(cos_sim_matrix1_coding$C2_MutationalPatterns, "MutationalPatterns")
colnames(C2_MP_m1_coding) <- c("cosine_similarity", "Tool")

C2_DS_m1_coding <- data.frame(cos_sim_matrix1_coding$C2_deconstructsigs, "deconstructsig")
colnames(C2_DS_m1_coding) <- c("cosine_similarity", "Tool")

C2_STL_m1_coding <- data.frame(cos_sim_matrix1_coding$C2_signature.tools.lib, "signature.tools.lib")
colnames(C2_STL_m1_coding) <- c("cosine_similarity", "Tool")

C2_SA_m1_coding <- data.frame(cos_sim_matrix1_coding$C2_sigprofilerassignment, "sigprofilerassignment")
colnames(C2_SA_m1_coding) <- c("cosine_similarity", "Tool")

cos_sim_cosmic2_matrix1_coding <- rbind(C2_MP_m1_coding,C2_DS_m1_coding,C2_STL_m1_coding,C2_SA_m1_coding) ### TO PLOT ###

C2_MP_m2_coding <- data.frame(cos_sim_matrix2_coding$C2_MutationalPatterns, "MutationalPatterns")
colnames(C2_MP_m2_coding) <- c("cosine_similarity", "Tool")

C2_DS_m2_coding <- data.frame(cos_sim_matrix2_coding$C2_deconstructsigs, "deconstructsig")
colnames(C2_DS_m2_coding) <- c("cosine_similarity", "Tool")

C2_STL_m2_coding <- data.frame(cos_sim_matrix2_coding$C2_signature.tools.lib, "signature.tools.lib")
colnames(C2_STL_m2_coding) <- c("cosine_similarity", "Tool")

C2_SA_m2_coding <- data.frame(cos_sim_matrix2_coding$C2_sigprofilerassignment, "sigprofilerassignment")
colnames(C2_SA_m2_coding) <- c("cosine_similarity", "Tool")

cos_sim_cosmic2_matrix2_coding <- rbind(C2_MP_m2_coding,C2_DS_m2_coding,C2_STL_m2_coding,C2_SA_m2_coding) ### TO PLOT ###

#Cosmic3.2

C3.2_MP_m1_coding <- data.frame(cos_sim_matrix1_coding$C3.2_MutationalPatterns, "MutationalPatterns")
colnames(C3.2_MP_m1_coding) <- c("cosine_similarity", "Tool")

C3.2_DS_m1_coding <- data.frame(cos_sim_matrix1_coding$C3.2_deconstructsigs, "deconstructsig")
colnames(C3.2_DS_m1_coding) <- c("cosine_similarity", "Tool")

C3.2_STL_m1_coding <- data.frame(cos_sim_matrix1_coding$C3.2_signature.tools.lib, "signature.tools.lib")
colnames(C3.2_STL_m1_coding) <- c("cosine_similarity", "Tool")

C3.2_SA_m1_coding <- data.frame(cos_sim_matrix1_coding$C3.2_sigprofilerassignment, "sigprofilerassignment")
colnames(C3.2_SA_m1_coding) <- c("cosine_similarity", "Tool")

cos_sim_cosmic3_2_matrix1_coding <- rbind(C3.2_MP_m1_coding,C3.2_DS_m1_coding,C3.2_STL_m1_coding,C3.2_SA_m1_coding) ### TO PLOT ###

C3.2_MP_m2_coding <- data.frame(cos_sim_matrix2_coding$C3.2_MutationalPatterns, "MutationalPatterns")
colnames(C3.2_MP_m2_coding) <- c("cosine_similarity", "Tool")

C3.2_DS_m2_coding <- data.frame(cos_sim_matrix2_coding$C3.2_deconstructsigs, "deconstructsig")
colnames(C3.2_DS_m2_coding) <- c("cosine_similarity", "Tool")

C3.2_STL_m2_coding <- data.frame(cos_sim_matrix2_coding$C3.2_signature.tools.lib, "signature.tools.lib")
colnames(C3.2_STL_m2_coding) <- c("cosine_similarity", "Tool")

C3.2_SA_m2_coding <- data.frame(cos_sim_matrix2_coding$C3.2_sigprofilerassignment, "sigprofilerassignment")
colnames(C3.2_SA_m2_coding) <- c("cosine_similarity", "Tool")

cos_sim_cosmic3_2_matrix2_coding <- rbind(C3.2_MP_m2_coding,C3.2_DS_m2_coding,C3.2_STL_m2_coding,C3.2_SA_m2_coding) ### TO PLOT ###

#non_coding mutations only

C2_MP_m1_non_coding <- data.frame(cos_sim_matrix1_coding$C2_MutationalPatterns, "MutationalPatterns")
colnames(C2_MP_m1_non_coding) <- c("cosine_similarity", "Tool")

C2_DS_m1_non_coding <- data.frame(cos_sim_matrix1_coding$C2_deconstructsigs, "deconstructsig")
colnames(C2_DS_m1_non_coding) <- c("cosine_similarity", "Tool")

C2_STL_m1_non_coding <- data.frame(cos_sim_matrix1_coding$C2_signature.tools.lib, "signature.tools.lib")
colnames(C2_STL_m1_non_coding) <- c("cosine_similarity", "Tool")

C2_SA_m1_non_coding <- data.frame(cos_sim_matrix1_coding$C2_sigprofilerassignment, "sigprofilerassignment")
colnames(C2_SA_m1_non_coding) <- c("cosine_similarity", "Tool")

cos_sim_cosmic2_matrix1_non_coding <- rbind(C2_MP_m1_coding,C2_DS_m1_coding,C2_STL_m1_coding,C2_SA_m1_coding) ### TO PLOT ###

C2_MP_m2_non_coding <- data.frame(cos_sim_matrix2_coding$C2_MutationalPatterns, "MutationalPatterns")
colnames(C2_MP_m2_non_coding) <- c("cosine_similarity", "Tool")

C2_DS_m2_non_coding <- data.frame(cos_sim_matrix2_coding$C2_deconstructsigs, "deconstructsig")
colnames(C2_DS_m2_non_coding) <- c("cosine_similarity", "Tool")

C2_STL_m2_non_coding <- data.frame(cos_sim_matrix2_coding$C2_signature.tools.lib, "signature.tools.lib")
colnames(C2_STL_m2_non_coding) <- c("cosine_similarity", "Tool")

C2_SA_m2_non_coding <- data.frame(cos_sim_matrix2_coding$C2_sigprofilerassignment, "sigprofilerassignment")
colnames(C2_SA_m2_non_coding) <- c("cosine_similarity", "Tool")

cos_sim_cosmic2_matrix2_non_coding <- rbind(C2_MP_m2_coding,C2_DS_m2_coding,C2_STL_m2_coding,C2_SA_m2_coding) ### TO PLOT ###

#Cosmic3.2

C3.2_MP_m1_non_coding <- data.frame(cos_sim_matrix1_non_coding$C3.2_MutationalPatterns, "MutationalPatterns")
colnames(C3.2_MP_m1_non_coding) <- c("cosine_similarity", "Tool")

C3.2_DS_m1_non_coding <- data.frame(cos_sim_matrix1_non_coding$C3.2_deconstructsigs, "deconstructsig")
colnames(C3.2_DS_m1_non_coding) <- c("cosine_similarity", "Tool")

C3.2_STL_m1_non_coding <- data.frame(cos_sim_matrix1_non_coding$C3.2_signature.tools.lib, "signature.tools.lib")
colnames(C3.2_STL_m1_non_coding) <- c("cosine_similarity", "Tool")

C3.2_SA_m1_non_coding <- data.frame(cos_sim_matrix1_non_coding$C3.2_sigprofilerassignment, "sigprofilerassignment")
colnames(C3.2_SA_m1_non_coding) <- c("cosine_similarity", "Tool")

cos_sim_cosmic3_2_matrix1_non_coding <- rbind(C3.2_MP_m1_non_coding,C3.2_DS_m1_non_coding,C3.2_STL_m1_non_coding,C3.2_SA_m1_non_coding) ### TO PLOT ###

C3.2_MP_m2_non_coding <- data.frame(cos_sim_matrix2_non_coding$C3.2_MutationalPatterns, "MutationalPatterns")
colnames(C3.2_MP_m2_non_coding) <- c("cosine_similarity", "Tool")

C3.2_DS_m2_non_coding <- data.frame(cos_sim_matrix2_non_coding$C3.2_deconstructsigs, "deconstructsig")
colnames(C3.2_DS_m2_non_coding) <- c("cosine_similarity", "Tool")

C3.2_STL_m2_non_coding <- data.frame(cos_sim_matrix2_non_coding$C3.2_signature.tools.lib, "signature.tools.lib")
colnames(C3.2_STL_m2_non_coding) <- c("cosine_similarity", "Tool")

C3.2_SA_m2_non_coding <- data.frame(cos_sim_matrix2_non_coding$C3.2_sigprofilerassignment, "sigprofilerassignment")
colnames(C3.2_SA_m2_non_coding) <- c("cosine_similarity", "Tool")

cos_sim_cosmic3_2_matrix2_non_coding <- rbind(C3.2_MP_m2_non_coding,C3.2_DS_m2_non_coding,C3.2_STL_m2_non_coding,C3.2_SA_m2_non_coding) ### TO PLOT ###

####COSMIC 2####


#plot Cosine Similarity

cos_sim_matrix1_C2_plot <- ggplot(cos_sim_cosmic2_matrix1, aes(x = Tool, y = cosine_similarity)) +
	geom_boxplot(color="black", fill="grey") +
	ggtitle("All Mutations - cosmic v2 - Matrix1") +
	labs(x= "Mutational Signature Tools", y="Cosine Similarity") +
	scale_y_continuous(breaks = seq(0, 1, 0.1),
			   limits = c(0,1),
			   expand = c(0, 0))+
	geom_hline(yintercept=0.9, linetype="dashed") +
	theme_bw() +
	theme(axis.title.x = element_text(size = 5),
	      axis.title.y = element_text(size = 5),
	      axis.text = element_text(size = 5),
	      plot.title = element_text(hjust = 0.5, size = 6.5, face = "bold"),
	      panel.border = element_rect(color = "black", linewidth = 0.4))



cos_sim_matrix2_C2_plot <- ggplot(cos_sim_cosmic2_matrix2, aes(x = Tool, y = cosine_similarity)) +
	geom_boxplot(color="black", fill="grey") +
	ggtitle("All Mutations - cosmic v2 - Matrix2") +
	labs(x= "Mutational Signature Tools", y="Cosine Similarity") +
	scale_y_continuous(breaks = seq(0, 1, 0.1),
			   limits = c(0,1),
			   expand = c(0, 0))+ 
	geom_hline(yintercept=0.9, linetype="dashed") +
	theme_bw() +
	theme(axis.title.x = element_text(size = 5),
	      axis.title.y = element_text(size = 5), 
	      axis.text = element_text(size = 5),
	      plot.title = element_text(hjust = 0.5, size = 6.5, face = "bold"), 
	      panel.border = element_rect(color = "black", linewidth = 0.4))


# coding

cos_sim_matrix1_C2_coding_plot <- ggplot(cos_sim_cosmic2_matrix1_coding, aes(x = Tool, y = cosine_similarity)) +
	geom_boxplot(color="black", fill="grey") +
	ggtitle("Coding Mutations - cosmic v2 - Matrix1") +
	labs(x= "Mutational Signature Tools", y="Cosine Similarity") +
	scale_y_continuous(breaks = seq(0, 1, 0.1),
			   limits = c(0,1),
		           expand = c(0, 0))+
	geom_hline(yintercept=0.9, linetype="dashed") +
        theme_bw() +
	theme(axis.title.x = element_text(size = 5),
	      axis.title.y = element_text(size = 5),
	      axis.text = element_text(size = 5),
	      plot.title = element_text(hjust = 0.5, size = 6.5, face = "bold"),
	      panel.border = element_rect(color = "black", linewidth = 0.4))


cos_sim_matrix2_C2_coding_plot <- ggplot(cos_sim_cosmic2_matrix2_coding, aes(x = Tool, y = cosine_similarity)) +
	geom_boxplot(color="black", fill="grey") +
	ggtitle("Coding Mutations - cosmic v2 - Matrix2") +
	labs(x= "Mutational Signature Tools", y="Cosine Similarity") +
 	scale_y_continuous(breaks = seq(0, 1, 0.1),
   			    limits = c(0,1),
   			    expand = c(0, 0))+
	geom_hline(yintercept=0.9, linetype="dashed") +
	theme_bw() +
	theme(axis.title.x = element_text(size = 5),
	      axis.title.y = element_text(size = 5),
	      axis.text = element_text(size = 5),
	      plot.title = element_text(hjust = 0.5, size = 6.5, face = "bold"),
	      panel.border = element_rect(color = "black", linewidth = 0.4))

#non coding

cos_sim_matrix1_C2_non_coding_plot <- ggplot(cos_sim_cosmic2_matrix1_non_coding, aes(x = Tool, y = cosine_similarity)) +
	geom_boxplot(color="black", fill="grey") +
	ggtitle("Non-Coding Mutations - cosmic v2 - Matrix1") +
	labs(x= "Mutational Signature Tools", y="Cosine Similarity") +
 	scale_y_continuous(breaks = seq(0, 1, 0.1),
  			    limits = c(0,1),  
			    expand = c(0, 0))+
	geom_hline(yintercept=0.9, linetype="dashed") +
	theme_bw() +
	theme(axis.title.x = element_text(size = 5),
	      axis.title.y = element_text(size = 5),
	      axis.text = element_text(size = 5),
	      plot.title = element_text(hjust = 0.5, size = 6.5, face = "bold"),
	      panel.border = element_rect(color = "black", linewidth = 0.4))

cos_sim_matrix2_C2_non_coding_plot <- ggplot(cos_sim_cosmic2_matrix2_non_coding, aes(x = Tool, y = cosine_similarity)) +
	geom_boxplot(color="black", fill="grey") +
	ggtitle("Non-Coding Mutations - cosmic v2 - Matrix2") +
	labs(x= "Mutational Signature Tools", y="Cosine Similarity") +
 	scale_y_continuous(breaks = seq(0, 1, 0.1),
  			    limits = c(0,1),
  			    expand = c(0, 0))+
	geom_hline(yintercept=0.9, linetype="dashed") +
	theme_bw() +
	theme(axis.title.x = element_text(size = 5),
	      axis.title.y = element_text(size = 5),
	      axis.text = element_text(size = 5),
	      plot.title = element_text(hjust = 0.5, size = 6.5, face = "bold"),
	      panel.border = element_rect(color = "black", linewidth = 0.4))


ggarrange(cos_sim_matrix1_C2_plot, cos_sim_matrix2_C2_plot, cos_sim_matrix1_C2_coding_plot, cos_sim_matrix2_C2_coding_plot, cos_sim_matrix1_C2_non_coding_plot, cos_sim_matrix2_C2_non_coding_plot,
	  labels = c("A", "B", "C", "D", "E", "F"),
	  ncol = 2, nrow = 3,
	  font.label = list(size = 10, color = "black", face = "bold"))



ggsave(paste(working_dir, "Figures", "cosine_similarity_cosmic2.png", sep="/"), dpi=300)



####COSMIC 3.2####


#plot Cosine Similarity

cos_sim_matrix1_C3_2_plot <- ggplot(cos_sim_cosmic3_2_matrix1, aes(x = Tool, y = cosine_similarity)) +
        geom_boxplot(color="black", fill="grey") +
	ggtitle("All Mutations - cosmic v3.2 - Matrix1") +
        labs(x= "Mutational Signature Tools", y="Cosine Similarity") +
        scale_y_continuous(breaks = seq(0, 1, 0.1),
                           limits = c(0,1),
                           expand = c(0, 0))+
	geom_hline(yintercept=0.9, linetype="dashed") +
        theme_bw() +
        theme(axis.title.x = element_text(size = 5),
              axis.title.y = element_text(size = 5),
              axis.text = element_text(size = 5),
	      plot.title = element_text(hjust = 0.5, size = 6.5, face = "bold"),
              panel.border = element_rect(color = "black", linewidth = 0.4))



cos_sim_matrix2_C3_2_plot <- ggplot(cos_sim_cosmic3_2_matrix2, aes(x = Tool, y = cosine_similarity)) +
        geom_boxplot(color="black", fill="grey") +
        ggtitle("All Mutations - cosmic v3.2 - Matrix2") +
	labs(x= "Mutational Signature Tools", y="Cosine Similarity") +
 	scale_y_continuous(breaks = seq(0, 1, 0.1),
  			    limits = c(0,1),
   			    expand = c(0, 0))+
	geom_hline(yintercept=0.9, linetype="dashed") +
        theme_bw() +
        theme(axis.title.x = element_text(size = 5),
              axis.title.y = element_text(size = 5),
              axis.text = element_text(size = 5),
	      plot.title = element_text(hjust = 0.5, size = 6.5, face = "bold"),
              panel.border = element_rect(color = "black", linewidth = 0.4))

# coding

cos_sim_matrix1_C3_2_coding_plot <- ggplot(cos_sim_cosmic3_2_matrix1_coding, aes(x = Tool, y = cosine_similarity)) +
        geom_boxplot(color="black", fill="grey") +
	ggtitle("Coding Mutations - cosmic v3.2 - Matrix1") +
        labs(x= "Mutational Signature Tools", y="Cosine Similarity") +
  	scale_y_continuous(breaks = seq(0, 1, 0.1),
    			     limits = c(0,1),
      			     expand = c(0, 0))+
	geom_hline(yintercept=0.9, linetype="dashed") +
        theme_bw() +
        theme(axis.title.x = element_text(size = 5),
              axis.title.y = element_text(size = 5),
              axis.text = element_text(size = 5),
	      plot.title = element_text(hjust = 0.5, size = 6.5, face = "bold"),
              panel.border = element_rect(color = "black", linewidth = 0.4))


cos_sim_matrix2_C3_2_coding_plot <- ggplot(cos_sim_cosmic3_2_matrix2_coding, aes(x = Tool, y = cosine_similarity)) +
        geom_boxplot(color="black", fill="grey") +
        ggtitle("Coding Mutations - cosmic v3.2 - Matrix2") +
	labs(x= "Mutational Signature Tools", y="Cosine Similarity") +
 	scale_y_continuous(breaks = seq(0, 1, 0.1),
   			    limits = c(0,1),
  			    expand = c(0, 0))+
	geom_hline(yintercept=0.9, linetype="dashed") +
        theme_bw() +
        theme(axis.title.x = element_text(size = 5),
              axis.title.y = element_text(size = 5),
              axis.text = element_text(size = 5),
	      plot.title = element_text(hjust = 0.5, size = 6.5, face = "bold"),
              panel.border = element_rect(color = "black", linewidth = 0.4))

#non coding

cos_sim_matrix1_C3_2_non_coding_plot <- ggplot(cos_sim_cosmic3_2_matrix1_non_coding, aes(x = Tool, y = cosine_similarity)) +
        geom_boxplot(color="black", fill="grey") +
	ggtitle("Non-Coding Mutations - cosmic v3.2 - Matrix1") +
        labs(x= "Mutational Signature Tools", y="Cosine Similarity") +
 	scale_y_continuous(breaks = seq(0, 1, 0.1),
  			    limits = c(0,1),
  			    expand = c(0, 0))+
	geom_hline(yintercept=0.9, linetype="dashed") +
        theme_bw() +
        theme(axis.title.x = element_text(size = 5),
              axis.title.y = element_text(size = 5),
              axis.text = element_text(size = 5),
	      plot.title = element_text(hjust = 0.5, size = 6.5, face = "bold"),
              panel.border = element_rect(color = "black", linewidth = 0.4))


cos_sim_matrix2_C3_2_non_coding_plot <- ggplot(cos_sim_cosmic3_2_matrix2_non_coding, aes(x = Tool, y = cosine_similarity)) +
        geom_boxplot(color="black", fill="grey") +
        ggtitle("Non-Coding Mutations - cosmic v3.2 - Matrix2") +
	labs(x= "Mutational Signature Tools", y="Cosine Similarity") +
  	scale_y_continuous(breaks = seq(0, 1, 0.1),
   			     limits = c(0,1),
   			     expand = c(0, 0))+
	geom_hline(yintercept=0.9, linetype="dashed") +
        theme_bw() +
        theme(axis.title.x = element_text(size = 5),
              axis.title.y = element_text(size = 5),
              axis.text = element_text(size = 5),
	      plot.title = element_text(hjust = 0.5, size = 6.5, face = "bold"),
              panel.border = element_rect(color = "black", linewidth = 0.4))



ggarrange(cos_sim_matrix1_C3_2_plot, cos_sim_matrix2_C3_2_plot, cos_sim_matrix1_C3_2_coding_plot, cos_sim_matrix2_C3_2_coding_plot, cos_sim_matrix1_C3_2_non_coding_plot, cos_sim_matrix2_C3_2_non_coding_plot,
          labels = c("A", "B", "C", "D", "E", "F"),
          ncol = 2, nrow = 3,
          font.label = list(size = 10, color = "black", face = "bold"))


ggsave(paste(working_dir, "Figures", "cosine_similarity_cosmic3.2.png", sep="/"), dpi=300)


#SBS analysis

SBS_C2_plot <- ggplot(SBS_C2, aes(x = Tool, y = SBS_contr, color = factor(Condition))) +
        geom_boxplot() +
	ggtitle("All Mutations - cosmic v2") +
        labs(x= "Mutational Signature Tools", y="Signature Contribution", color="Conditions") +
        scale_y_continuous(breaks = seq(0, 1, 0.1),
                           limits = c(0,1),
                           expand = c(0, 0))+
        theme_bw() +
        scale_color_manual(values=c("goldenrod1", "firebrick4")) +
	theme(axis.title.x = element_text(size = 5),
              axis.title.y = element_text(size = 5),
              axis.text = element_text(size = 5),
	      plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
              panel.border = element_rect(color = "black", linewidth = 0.4))



SBS_C2_coding_plot <- ggplot(SBS_C2_coding, aes(x = Tool, y = SBS_contr, color = factor(Condition))) +
        geom_boxplot() +
	ggtitle("Coding Mutations - cosmic v2") +
        labs(x= "Mutational Signature Tools", y="Signature Contribution", color="Conditions") +
        scale_y_continuous(breaks = seq(0, 1, 0.1),
                           limits = c(0,1),
                           expand = c(0, 0))+
        theme_bw() +
	scale_color_manual(values=c("goldenrod1", "firebrick4")) +
        theme(axis.title.x = element_text(size = 5),
              axis.title.y = element_text(size = 5),
              axis.text = element_text(size = 5),
	      plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
              panel.border = element_rect(color = "black", linewidth = 0.4))



SBS_C2_non_coding_plot <- ggplot(SBS_C2_non_coding, aes(x = Tool, y = SBS_contr, color = factor(Condition))) +
        geom_boxplot() +
	ggtitle("Non-Coding Mutations - cosmic v2") +
        labs(x= "Mutational Signature Tools", y="Signature Contribution", color="Conditions") +
        scale_y_continuous(breaks = seq(0, 1, 0.1),
                           limits = c(0,1),
                           expand = c(0, 0))+
        theme_bw() +
        scale_color_manual(values=c("goldenrod1", "firebrick4")) +	
	theme(axis.title.x = element_text(size = 5),
              axis.title.y = element_text(size = 5),
              axis.text = element_text(size = 5),
	      plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
              panel.border = element_rect(color = "black", linewidth = 0.4))


ggarrange(SBS_C2_plot, SBS_C2_coding_plot, SBS_C2_non_coding_plot,
          labels = c("A", "B", "C"),
          ncol = 1, nrow = 3,
          font.label = list(size = 10, color = "black", face = "bold"))


ggsave(paste(working_dir, "Figures", "SBS_contribution_cosmic2.png", sep="/"), dpi=300)


#SBS analysis Cosmic v3.2

SBS_C3_2_plot <- ggplot(SBS_C3_2, aes(x = Tool, y = SBS_contr, color = factor(Condition))) +
        geom_boxplot() +
	ggtitle("All Mutations - cosmic v3.2") +
        labs(x= "Mutational Signature Tools", y="Signature Contribution", color="Conditions") +
        scale_y_continuous(breaks = seq(0, 1, 0.1),
                           limits = c(0,1),
                           expand = c(0, 0))+
        theme_bw() +
	scale_color_manual(values=c("goldenrod1", "firebrick4")) +
        theme(axis.title.x = element_text(size = 5),
              axis.title.y = element_text(size = 5),
              axis.text = element_text(size = 5),
	      plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
              panel.border = element_rect(color = "black", linewidth = 0.4))



SBS_C3_2_coding_plot <- ggplot(SBS_C3_2_coding, aes(x = Tool, y = SBS_contr, color = factor(Condition))) +
        geom_boxplot() +
	ggtitle("Coding Mutations - cosmic v3.2") +
        labs(x= "Mutational Signature Tools", y="Signature Contribution", color="Conditions") +
        scale_y_continuous(breaks = seq(0, 1, 0.1),
                           limits = c(0,1),
                           expand = c(0, 0))+
        theme_bw() +
        scale_color_manual(values=c("goldenrod1", "firebrick4")) +	
	theme(axis.title.x = element_text(size = 5),
              axis.title.y = element_text(size = 5),
              axis.text = element_text(size = 5),
	      plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
              panel.border = element_rect(color = "black", linewidth = 0.4))



SBS_C3_2_non_coding_plot <- ggplot(SBS_C3_2_non_coding, aes(x = Tool, y = SBS_contr, color = factor(Condition))) +
        geom_boxplot() +
	ggtitle("Non-Coding Mutations - cosmic v3.2") +
        labs(x= "Mutational Signature Tools", y="Signature Contribution", color="Conditions") +
        scale_y_continuous(breaks = seq(0, 1, 0.1),
                           limits = c(0,1),
                           expand = c(0, 0))+
        theme_bw() +
        scale_color_manual(values=c("goldenrod1", "firebrick4")) +
        theme(axis.title.x = element_text(size = 5),
              axis.title.y = element_text(size = 5),
              axis.text = element_text(size = 5),
	      plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
              panel.border = element_rect(color = "black", linewidth = 0.4))


ggarrange(SBS_C3_2_plot, SBS_C3_2_coding_plot, SBS_C3_2_non_coding_plot,
          labels = c("A", "B", "C"),
          ncol = 1, nrow = 3,
          font.label = list(size = 10, color = "black", face = "bold"))


ggsave(paste(working_dir, "Figures", "SBS_contribution_cosmic3.2.png", sep="/"), dpi=300)


