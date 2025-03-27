#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(signature.tools.lib))
suppressPackageStartupMessages(library(MutationalPatterns))

args <- commandArgs()

mut_matrix_file <- args[6]
script_path <- args[7]
working_dir <- args[8]
condition <- args[9]
genome <- args[10]

cosmic2_file <- paste("COSMIC_v2_SBS_", genome, ".txt", sep="")
cosmic3.2_file <- paste("COSMIC_v3.2_SBS_", genome, ".txt", sep="")

COSMIC2 <- read.table(paste(script_path, "cosmic_ref", cosmic2_file, sep = "/"), header = T, row.names = 1)
COSMIC2 <- as.matrix(sapply(COSMIC2, as.numeric))

COSMIC3.2 <- read.table(paste(script_path, "cosmic_ref", cosmic3.2_file, sep = "/"), header = T, row.names = 1)
COSMIC3.2 <- as.matrix(sapply(COSMIC3.2, as.numeric))


#COSMIC2 <- read.table(paste(script_path, "cosmic_ref/COSMIC_v2_SBS_GRCh38.txt", sep = "/"), header = T, row.names = 1)
#COSMIC2 <- as.matrix(sapply(COSMIC2, as.numeric))

#COSMIC3.2 <- read.table(paste(script_path, "cosmic_ref/COSMIC_v3.2_SBS_GRCh38.txt", sep = "/"), header = T, row.names = 1)
#COSMIC3.2 <- as.matrix(sapply(COSMIC3.2, as.numeric))

mut_matrix <- read.table(mut_matrix_file, header = T, row.names = 1)

#Fitting Analysis

Fitting_mut_matrix_Cosmic2 <- SignatureFit(mut_matrix, COSMIC2)
Fitting_mut_matrix_Cosmic3.2 <- SignatureFit(mut_matrix, COSMIC3.2)

#Reconstructed CRC

reconstructed_Cosmic2_STL <- data.frame(matrix(NA, nrow = length(seq(1:96)), ncol = length(seq(1:length(colnames(mut_matrix))))))
rownames(reconstructed_Cosmic2_STL) <- rownames(mut_matrix)
colnames(reconstructed_Cosmic2_STL) <- colnames(mut_matrix)

results_values_Cosmic2_STL <- data.frame(matrix(NA, nrow = length(seq(1:96)), ncol = length(seq(1:30))))
rownames(results_values_Cosmic2_STL) <- rownames(mut_matrix)
colnames(results_values_Cosmic2_STL) <- colnames(COSMIC2)

for (j in 1:length(colnames(mut_matrix))){
	for (i in 1:30){
		results_values_Cosmic2_STL[,i] <- Fitting_mut_matrix_Cosmic2[i,j]*(COSMIC2[,i])
	}
	reconstructed_Cosmic2_STL[,j] <- apply(results_values_Cosmic2_STL, 1, sum)
}



#Cosine Similarity CRC

cos_sim_mut_matrix_COSMIC2 <- diag(cos_sim_matrix(mut_matrix, reconstructed_Cosmic2_STL))
names(cos_sim_mut_matrix_COSMIC2) <- colnames(mut_matrix)
write.table(cos_sim_mut_matrix_COSMIC2, paste(working_dir, "Cosine_Similarity", condition, "cos_similarity_C2_signature.tools.lib.txt", sep = "/"), sep = "\t", quote = F, col.names = F)

#Reconstructed CRC

reconstructed_Cosmic3.2_STL <- data.frame(matrix(NA, nrow = length(seq(1:96)), ncol = length(seq(1:length(colnames(mut_matrix))))))
rownames(reconstructed_Cosmic3.2_STL) <- rownames(mut_matrix)
colnames(reconstructed_Cosmic3.2_STL) <- colnames(mut_matrix)

results_values_Cosmic3.2_STL <- data.frame(matrix(NA, nrow = length(seq(1:96)), ncol = length(seq(1:78))))
rownames(results_values_Cosmic3.2_STL) <- rownames(mut_matrix)
colnames(results_values_Cosmic3.2_STL) <- colnames(COSMIC3.2)

for (j in 1:length(colnames(mut_matrix))){
	for (i in 1:78){
		results_values_Cosmic3.2_STL[,i] <- Fitting_mut_matrix_Cosmic3.2[i,j]*(COSMIC3.2[,i])
	}
	reconstructed_Cosmic3.2_STL[,j] <- apply(results_values_Cosmic3.2_STL, 1, sum)
}



#Cosine Similarity CRC

cos_sim_mut_matrix_COSMIC3.2 <- diag(cos_sim_matrix(mut_matrix, reconstructed_Cosmic3.2_STL))
names(cos_sim_mut_matrix_COSMIC3.2) <- colnames(mut_matrix)
write.table(cos_sim_mut_matrix_COSMIC3.2, paste(working_dir, "Cosine_Similarity", condition, "cos_similarity_C3.2_signature.tools.lib.txt", sep = "/"), sep = "\t", quote = F, col.names = F)

#Write Normalized Fitting Matrix

Fitting_mut_matrix_Cosmic2_norm <- Fitting_mut_matrix_Cosmic2/colSums(Fitting_mut_matrix_Cosmic2)[col(Fitting_mut_matrix_Cosmic2)]
Fitting_mut_matrix_Cosmic3.2_norm <- Fitting_mut_matrix_Cosmic3.2/colSums(Fitting_mut_matrix_Cosmic3.2)[col(Fitting_mut_matrix_Cosmic3.2)]

Fitting_mut_matrix_Cosmic2_norm <- as.data.frame(Fitting_mut_matrix_Cosmic2_norm)
Fitting_mut_matrix_Cosmic3.2_norm <- as.data.frame(Fitting_mut_matrix_Cosmic3.2_norm)


write.table(Fitting_mut_matrix_Cosmic2_norm, paste(working_dir, "SBS_signature_contributions", condition, "sbs_sig_C2_signature.tools.lib.txt", sep = "/"), sep = "\t", quote = F, col.names = NA, row.names = T)
write.table(Fitting_mut_matrix_Cosmic3.2_norm, paste(working_dir, "SBS_signature_contributions", condition, "sbs_sig_C3.2_signature.tools.lib.txt", sep = "/"), sep = "\t", quote = F, col.names = NA, row.names = T)


