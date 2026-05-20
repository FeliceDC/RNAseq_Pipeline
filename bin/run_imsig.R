#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
counts_file = args[1]

if (!requireNamespace("imsig", quietly = TRUE)) {
    install.packages("imsig") }
library(imsig)

cat("Reading the counts file...\n")
counts_data <- read.table(counts_file, header=TRUE, row.names=1, sep="\t", check.names=FALSE)

if("Length" %in% colnames(counts_data)){
    exp_matrix <- counts_data[, 6:ncol(counts_data)]
} else {
    exp_matrix <- counts_data
}

cat("ImSig execution...\n")
imsig_res <- imsig(exp = exp_matrix, r = 0.7)

cat("Saving results...\n")
write.csv(imsig_res, "ImSig_results.csv", row.names=TRUE)
