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



library(ggplot2)
library(tidyr)


imsig_perc <- sweep(imsig_res, 1, rowSums(imsig_res), FUN="/") * 100


imsig_perc$Sample <- rownames(imsig_perc)

df_long <- pivot_longer(imsig_perc, cols = -Sample, names_to = "CellType", values_to = "Percentage")


p <- ggplot(df_long, aes(x = Sample, y = Percentage, fill = CellType)) +
    geom_bar(stat = "identity", position = "fill") +
    theme_minimal() +

    theme(
        axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
        title = element_text(face = "bold")
    ) +
    labs(
        title = "Relative Composition of Immune Infiltrate (ImSig)",
        x = "Samples",
        y = "Relative Proportion",
        fill = "Cell type"
    ) +
    scale_y_continuous(labels = scales::percent_format())

ggsave("ImSig_plot.pdf", plot = p, width = 10, height = 7)
