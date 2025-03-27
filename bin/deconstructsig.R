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
genome <- args[10]

cosmic2_file <- paste("COSMIC_v2_SBS_", genome, ".txt", sep="")
cosmic3.2_file <- paste("COSMIC_v3.2_SBS_", genome, ".txt", sep="")

COSMIC2 <- read.table(paste(script_path, "cosmic_ref", cosmic2_file, sep = "/"), header = T, row.names = 1)
#COSMIC2 <- as.matrix(sapply(COSMIC2, as.numeric))

COSMIC3.2 <- read.table(paste(script_path, "cosmic_ref", cosmic3.2_file, sep = "/"), header = T, row.names = 1)
#COSMIC3.2 <- as.matrix(sapply(COSMIC3.2, as.numeric))


#COSMIC2 <- read.table(paste(script_path, "cosmic_ref/COSMIC_v2_SBS_GRCh38.txt", sep = "/"), header = T, row.names = 1)
COSMIC2 <- t(COSMIC2)
COSMIC2 <- as.data.frame(COSMIC2)

#COSMIC3.2 <- read.table(paste(script_path, "cosmic_ref/COSMIC_v3.2_SBS_GRCh38.txt", sep = "/"), header = T, row.names = 1)
COSMIC3.2 <- t(COSMIC3.2)
COSMIC3.2 <- as.data.frame(COSMIC3.2)

mut_matrix <- read.table(mut_matrix_file, header = T, row.names = 1)

sample_list <- colnames(mut_matrix)
mut_matrix <- t(mut_matrix)

mut_matrix <- mut_matrix[, colnames(randomly.generated.tumors)]

results_cosmic2 <- vector(mode = "list", length = length(colnames(mut_matrix)))
results_cosmic3.2 <- vector(mode = "list", length = length(colnames(mut_matrix)))

mut_matrix <- as.data.frame(mut_matrix)

for (i in 1:length(rownames(mut_matrix))){results_cosmic2[[i]] <- whichSignatures(tumor.ref = mut_matrix, sample.id = rownames(mut_matrix)[i], signatures.ref = COSMIC2[,colnames(randomly.generated.tumors)], contexts.needed = TRUE)}
for (i in 1:length(rownames(mut_matrix))){results_cosmic3.2[[i]] <- whichSignatures(tumor.ref = mut_matrix, sample.id = rownames(mut_matrix)[i], signatures.ref = COSMIC3.2[,colnames(randomly.generated.tumors)], contexts.needed = TRUE)}


results_cosmic2_weight <- vector(mode = "list", length = length(colnames(mut_matrix)))
results_cosmic3.2_weight <- vector(mode = "list", length = length(colnames(mut_matrix)))

for (i in 1:length(rownames(mut_matrix))){results_cosmic2_weight[[i]] <- results_cosmic2[[i]]$weights}
for (i in 1:length(rownames(mut_matrix))){results_cosmic3.2_weight[[i]] <- results_cosmic3.2[[i]]$weights}


results_cosmic2_weight_df <- results_cosmic2_weight[[1]]
results_cosmic3.2_weight_df <- results_cosmic3.2_weight[[1]]


for (i in 2:length(rownames(mut_matrix))){results_cosmic2_weight_df <- rbind(results_cosmic2_weight_df, results_cosmic2_weight[[i]])}
for (i in 2:length(rownames(mut_matrix))){results_cosmic3.2_weight_df <- rbind(results_cosmic3.2_weight_df, results_cosmic3.2_weight[[i]])}

results_cosmic2_weight_df <- t(results_cosmic2_weight_df)
results_cosmic3.2_weight_df <- t(results_cosmic3.2_weight_df)

mut_matrix <- read.table(mut_matrix_file, header = T, row.names = 1)
COSMIC2 <- read.table(paste(script_path, "cosmic_ref/COSMIC_v2_SBS_GRCh38.txt", sep = "/"), header = T, row.names = 1)
COSMIC2 <- as.data.frame(COSMIC2)
COSMIC3.2 <- read.table(paste(script_path, "cosmic_ref/COSMIC_v3.2_SBS_GRCh38.txt", sep = "/"), header = T, row.names = 1)
COSMIC3.2 <- as.data.frame(COSMIC3.2)


#Reconstructed CRC cosmic2

reconstructed_Cosmic2 <- data.frame(matrix(NA, nrow = length(seq(1:96)), ncol = length(seq(1:length(colnames(mut_matrix))))))
rownames(reconstructed_Cosmic2) <- rownames(mut_matrix)
colnames(reconstructed_Cosmic2) <- colnames(mut_matrix)

results_values_Cosmic2 <- data.frame(matrix(NA, nrow = length(seq(1:96)), ncol = length(seq(1:30))))
rownames(results_values_Cosmic2) <- rownames(mut_matrix)
colnames(results_values_Cosmic2) <- colnames(COSMIC2)


for (j in 1:length(colnames(mut_matrix))){
	for (i in 1:30){
  		results_values_Cosmic2[,i] <- results_cosmic2_weight_df[i,j]*(COSMIC2[,i])
  	}
	reconstructed_Cosmic2[,j] <- apply(results_values_Cosmic2, 1, sum)
}

#Cosine Similarity CRC cosmic2

cos_sim_mut_matrix_COSMIC2 <- diag(cos_sim_matrix(mut_matrix, reconstructed_Cosmic2))
names(cos_sim_mut_matrix_COSMIC2) <- colnames(mut_matrix)

cos_sim_mut_matrix_COSMIC2 <- as.data.frame(cos_sim_mut_matrix_COSMIC2)

write.table(cos_sim_mut_matrix_COSMIC2, paste(working_dir, "Cosine_Similarity", condition, "cos_similarity_C2_deconstructsigs.txt", sep = "/"), sep = "\t", quote = F, col.names = F)


#Reconstructed CRC cosmic3.2

reconstructed_Cosmic3.2 <- data.frame(matrix(NA, nrow = length(seq(1:96)), ncol = length(seq(1:length(colnames(mut_matrix))))))
rownames(reconstructed_Cosmic3.2) <- rownames(mut_matrix)
colnames(reconstructed_Cosmic3.2) <- colnames(mut_matrix)

results_values_Cosmic3.2 <- data.frame(matrix(NA, nrow = length(seq(1:96)), ncol = length(seq(1:78))))
rownames(results_values_Cosmic3.2) <- rownames(mut_matrix)
colnames(results_values_Cosmic3.2) <- colnames(COSMIC3.2)


for (j in 1:length(colnames(mut_matrix))){
	for (i in 1:78){
		results_values_Cosmic3.2[,i] <- results_cosmic3.2_weight_df[i,j]*(COSMIC3.2[,i])
	}
	reconstructed_Cosmic3.2[,j] <- apply(results_values_Cosmic3.2, 1, sum)
}

#Cosine Similarity CRC cosmic3.2

cos_sim_mut_matrix_COSMIC3.2 <- diag(cos_sim_matrix(mut_matrix, reconstructed_Cosmic3.2))
names(cos_sim_mut_matrix_COSMIC3.2) <- colnames(mut_matrix)

cos_sim_mut_matrix_COSMIC3.2 <- as.data.frame(cos_sim_mut_matrix_COSMIC3.2)

write.table(cos_sim_mut_matrix_COSMIC3.2, paste(working_dir, "Cosine_Similarity", condition, "cos_similarity_C3.2_deconstructsigs.txt", sep = "/"), sep = "\t", quote = F, col.names = F)


#Write weights

write.table(results_cosmic2_weight_df, paste(working_dir, "SBS_signature_contributions", condition, "sbs_sig_C2_deconstructsigs.txt", sep = "/"), sep = "\t", quote = F, col.names = NA)
write.table(results_cosmic3.2_weight_df, paste(working_dir, "SBS_signature_contributions", condition, "sbs_sig_C3.2_deconstructsigs.txt", sep = "/"), sep = "\t", quote = F, col.names = NA)

