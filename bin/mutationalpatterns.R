#!/usr/bin/env Rscript

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

mut_matrix <- read.table(mut_matrix_file, header = T, row.names = 1)

#Fitting Analysis

fitting_mut_matrix_C2 <- fit_to_signatures(mut_matrix, COSMIC2)
fitting_mut_matrix_C3.2 <- fit_to_signatures(mut_matrix, COSMIC3.2) 

sbs_matrix_C2 <- fitting_mut_matrix_C2$contribution
sbs_matrix_C3.2 <- fitting_mut_matrix_C3.2$contribution

sbs_matrix_C2_norm <- sbs_matrix_C2/colSums(sbs_matrix_C2)[col(sbs_matrix_C2)]
sbs_matrix_C3.2_norm <- sbs_matrix_C3.2/colSums(sbs_matrix_C3.2)[col(sbs_matrix_C3.2)]

sbs_matrix_C2_norm <- as.data.frame(sbs_matrix_C2_norm)
sbs_matrix_C3.2_norm <- as.data.frame(sbs_matrix_C3.2_norm)

#Cosine Similarity

cos_sim_sbs_matrix_C2 <- diag(cos_sim_matrix(mut_matrix, fitting_mut_matrix_C2$reconstructed))
cos_sim_sbs_matrix_C3.2 <- diag(cos_sim_matrix(mut_matrix, fitting_mut_matrix_C3.2$reconstructed))

write.table(cos_sim_sbs_matrix_C2, paste(working_dir, "Cosine_Similarity", condition, "cos_similarity_C2_MutationalPatterns.txt", sep = "/"), sep = "\t", quote = F, col.names = F)
write.table(cos_sim_sbs_matrix_C3.2, paste(working_dir, "Cosine_Similarity", condition, "cos_similarity_C3.2_MutationalPatterns.txt", sep = "/"), sep = "\t", quote = F, col.names = F)


#Signatures Contribution

write.table(sbs_matrix_C2_norm, paste(working_dir, "SBS_signature_contributions", condition, "sbs_sig_C2_MutationalPatterns.txt", sep = "/"), sep = "\t", quote = F, col.names = NA)
write.table(sbs_matrix_C3.2_norm, paste(working_dir, "SBS_signature_contributions", condition, "sbs_sig_C3.2_MutationalPatterns.txt", sep = "/"), sep = "\t", quote = F, col.names = NA)

