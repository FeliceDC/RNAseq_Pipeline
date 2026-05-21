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

imsig_perc <- sweep(imsig_res, 2, colSums(imsig_res), FUN="/") * 100
imsig_perc$CellType <- rownames(imsig_perc)
df_long <- pivot_longer(imsig_perc, cols = -CellType, names_to = "Sample", values_to = "Percentage")

p <- ggplot(df_long, aes(x = Sample, y = Percentage, fill = CellType)) +
    geom_bar(stat = "identity", position = "fill") + 
    coord_flip() + # Barra orizzontale, perfetta per nomi lunghi SRR
    theme_minimal() +
    theme(
        title = element_text(face = "bold"),
        legend.position = "right" # Mettiamo la legenda a destra
    ) +
    labs(
        title = "Relative Composition of Immune Infiltrate (ImSig)",
        subtitle = "Visualization of the composition for each sample",
        x = "Samples", # Etichetta corretta per asse Y (dopo coord_flip)
        y = "Relative Proportion", # Etichetta corretta per asse X (dopo coord_flip)
        fill = "Cell type"
    ) +
    scale_y_continuous(labels = scales::percent_format())

ggsave("ImSig_plot.pdf", plot = p, width = 12, height = 8)
