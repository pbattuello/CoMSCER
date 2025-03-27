#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(deconstructSigs))
suppressPackageStartupMessages(library(BSgenome))
suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg38))
suppressPackageStartupMessages(library(MutationalPatterns))
args <- commandArgs()

mut_matrix_file <- args[6]
script_path <- args[7]
working_dir <- args[8]
condition <- args[9]

COSMIC2 <- read.table(paste(script_path, "cosmic_ref/COSMIC_v2_SBS_GRCh38.txt", sep = "/"), header = T, row.names = 1)
COSMIC2 <- as.data.frame(COSMIC2)

COSMIC3.2 <- read.table(paste(script_path, "cosmic_ref/COSMIC_v3.2_SBS_GRCh38.txt", sep = "/"), header = T, row.names = 1)
COSMIC3.2 <- as.data.frame(COSMIC3.2)

mut_matrix <- read.table(mut_matrix_file, header = T, row.names = 1)

weigth_matrix_C2 <- read.table(paste(working_dir, "TMP" , condition,"sbs_sig_C2_sigprofilerassignment.txt", sep = "/"), header = T, row.names = 1)
weigth_matrix_C3.2 <- read.table(paste(working_dir, "TMP", condition, "sbs_sig_C3.2_sigprofilerassignment.txt", sep = "/"), header = T, row.names = 1)

weigth_matrix_C2 <- t(weigth_matrix_C2)
weigth_matrix_C3.2 <- t(weigth_matrix_C3.2)

weight_matrix_C2_norm <- weigth_matrix_C2/colSums(weigth_matrix_C2)[col(weigth_matrix_C2)]
weight_matrix_C3.2_norm <- weigth_matrix_C3.2/colSums(weigth_matrix_C3.2)[col(weigth_matrix_C3.2)]

#Reconstructed CRC cosmic2

reconstructed_Cosmic2 <- data.frame(matrix(NA, nrow = length(seq(1:96)), ncol = length(seq(1:length(colnames(mut_matrix))))))
rownames(reconstructed_Cosmic2) <- rownames(mut_matrix)
colnames(reconstructed_Cosmic2) <- colnames(mut_matrix)

results_values_Cosmic2 <- data.frame(matrix(NA, nrow = length(seq(1:96)), ncol = length(seq(1:30))))
rownames(results_values_Cosmic2) <- rownames(mut_matrix)
colnames(results_values_Cosmic2) <- colnames(COSMIC2)

for (j in 1:length(colnames(mut_matrix))){
	for (i in 1:30){
  		results_values_Cosmic2[,i] <- weight_matrix_C2_norm[i,j]*(COSMIC2[,i])
  	}
	reconstructed_Cosmic2[,j] <- apply(results_values_Cosmic2, 1, sum)
}

#Cosine Similarity CRC cosmic2

cos_sim_mut_matrix_COSMIC2 <- diag(cos_sim_matrix(mut_matrix, reconstructed_Cosmic2))
names(cos_sim_mut_matrix_COSMIC2) <- colnames(mut_matrix)

cos_sim_mut_matrix_COSMIC2 <- as.data.frame(cos_sim_mut_matrix_COSMIC2)
write.table(cos_sim_mut_matrix_COSMIC2, paste(working_dir, "Cosine_Similarity", condition, "cos_similarity_C2_sigprofilerassignment.txt", sep = "/"), sep = "\t", quote = F, col.names = F)


#Reconstructed CRC cosmic3.2

reconstructed_Cosmic3.2 <- data.frame(matrix(NA, nrow = length(seq(1:96)), ncol = length(seq(1:length(colnames(mut_matrix))))))
rownames(reconstructed_Cosmic3.2) <- rownames(mut_matrix)
colnames(reconstructed_Cosmic3.2) <- colnames(mut_matrix)

results_values_Cosmic3.2 <- data.frame(matrix(NA, nrow = length(seq(1:96)), ncol = length(seq(1:78))))
rownames(results_values_Cosmic3.2) <- rownames(mut_matrix)
colnames(results_values_Cosmic3.2) <- colnames(COSMIC3.2)


for (j in 1:length(colnames(mut_matrix))){
	for (i in 1:78){
		results_values_Cosmic3.2[,i] <- weight_matrix_C3.2_norm[i,j]*(COSMIC3.2[,i])
	}
	reconstructed_Cosmic3.2[,j] <- apply(results_values_Cosmic3.2, 1, sum)
}

#Cosine Similarity CRC cosmic3.2

cos_sim_mut_matrix_COSMIC3.2 <- diag(cos_sim_matrix(mut_matrix, reconstructed_Cosmic3.2))
names(cos_sim_mut_matrix_COSMIC3.2) <- colnames(mut_matrix)

cos_sim_mut_matrix_COSMIC3.2 <- as.data.frame(cos_sim_mut_matrix_COSMIC3.2)

write.table(cos_sim_mut_matrix_COSMIC3.2, paste(working_dir, "Cosine_Similarity", condition, "cos_similarity_C3.2_sigprofilerassignment.txt", sep = "/"), sep = "\t", quote = F, col.names = F)

#Write weights

write.table(weight_matrix_C2_norm, paste(working_dir, "SBS_signature_contributions", condition, "sbs_sig_C2_sigprofilerassignment.txt", sep = "/"), sep = "\t", quote = F, col.names = NA)
write.table(weight_matrix_C3.2_norm, paste(working_dir, "SBS_signature_contributions", condition, "sbs_sig_C3.2_sigprofilerassignment.txt", sep = "/"), sep = "\t", quote = F, col.names = NA)

