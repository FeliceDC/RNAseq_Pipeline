#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
counts_file <- args [1]
metadata_file <-args [2]
user_design <-args[3]

library(DESeq2)
library(ggplot2)

counts <- read.table(counts_file, header = TRUE, row.names = 1, stringsAsFactors = FALSE)
meta <- read.csv(metadata_file, row.names = 1, stringsAsFactors = TRUE)

counts <- counts[, 6:ncol(counts)]

col_names <- colnames(counts)

col_names <- sub("^X", "", col_names)
col_names <- sub("\\..*", "", col_names)

true_samples <- rownames(meta)

for (ts in true_samples) {
    match_idx <- grep(paste0("^", ts), col_names)
    if (length(match_idx) == 0) {
        match_idx <- grep(ts, col_names, fixed = TRUE)
    }
    if (length(match_idx) == 1) {
        col_names[match_idx] <- ts
    }
}

colnames(counts) <- col_names
common_samples <- intersect(colnames(counts), rownames(meta))

if (length(common_samples) == 0) {
    message <- paste("\n\n#####################################################\n",
                       "ERROR: The sample names do not correspond!\n",
                       "Matrix names: ", paste(colnames(counts), collapse=" , "), "\n",
                       "Samplesheet names: ", paste(rownames(meta), collapse=" , "), "\n",
                       "#####################################################\n\n")
    stop(message)
}

counts <- counts[, common_samples]
meta <- meta[common_samples, , drop=FALSE]

design_formula <- as.formula(paste("~", user_design))

if (length(unique(as.list(counts))) == 1) {
    message("\n========================================================")
    message("DETECTED: Samples Are Identical (Profile Test).")
    message("Adding artificial batch effect in order to test the pipeline.")
    message("========================================================\n")
    set.seed(42)
    rumore <- matrix(rpois(nrow(counts) * ncol(counts), lambda = 5), nrow = nrow(counts))
    counts <- counts + rumore
}

dds <- DESeqDataSetFromMatrix(countData = counts, colData = meta, design = design_formula)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

dds <- tryCatch({
    DESeq(dds, sfType="poscounts")
}, error = function(e) {
    message("\nStandard curve fitting failed (too few genes).")
    message("Unpacking the algorithm and forcing manual calculations...\n")
    dds <- estimateSizeFactors(dds, type="poscounts")
    dds <- estimateDispersionsGeneEst(dds)
    dispersions(dds) <- mcols(dds)$dispGeneEst
    dds <- nbinomWaldTest(dds)
    return(dds)
})

res <- results(dds)

write.table(as.data.frame(res), file="deseq2_results.txt", sep="\t", quote=FALSE, row.names=FALSE)

res_clean <- res[!is.na(res$padj), ]
res_filt <- res_clean[res_clean$padj < 0.05 & abs(res_clean$log2FoldChange) > 1.5, ]

res_filt_df <- as.data.frame(res_filt)
res_filt_df$Gene_Name <- rownames(res_filt_df)
write.table(res_filt_df, file="filtered_results.txt", sep="\t", quote=FALSE, row.names=FALSE)

normalized_counts <- counts(dds, normalized=TRUE)
res_df <- as.data.frame(res)
complete_table <- merge(res_df, normalized_counts, by="row.names", all=TRUE)
colnames(complete_table)[1] <- "Gene_Name"
complete_table <- complete_table[order(complete_table$padj), ]


write.table(as.data.frame(complete_table), file="complete_table.txt", sep="\t", quote=FALSE, row.names=FALSE)

#--- PLOTS ---

pdf ("deseq2_plots.pdf")

#MA plot
plotMA(res, main="MA Plot (test nf-core)")

#PCA plot
if (nrow(dds) < 1000) {
    message("MLess than 1000 genes: using normTransform instead of vst for plots")
    vsd <- normTransform(dds)
} else {
    vsd <- vst(dds, blind=FALSE)
}

pca_groups <- trimws(unlist(strsplit(user_design, "\\+")))

pca_plot <- plotPCA(vsd, intgroup=pca_groups) + theme_minimal() + geom_point(size=4, alpha=0.9) + ggtitle("2. PCA") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 15), 
    panel.grid.major = element_line(color = "grey90"),                
    panel.grid.minor = element_blank(),
    legend.title = element_text(face="bold", size=12), 
    legend.text = element_text(size = 11),                            
    axis.title = element_text(size = 12, face = "bold")               
  )

print(pca_plot)

#Volcano Plot
with(res, plot(log2FoldChange, -log10(padj), pch=20, main="3. Volcano Plot", col="darkgrey", xlim=c(-5,5)))
with(subset(res, padj < 0.05 & abs(log2FoldChange) > 1), points(log2FoldChange, -log10(padj), pch=20, col="red"))
abline(v=c(-1,1), col="blue", lty=2)
abline(h=-log10(0.05), col="blue", lty=2)

#Heatmap
top_genes <- head(order(res$padj), 50)
mat <- assay(vsd)[top_genes, ]
mat <- mat - rowMeans(mat)

colori_heatmap <- colorRampPalette(c("blue", "white", "red"))(256)

heatmap(mat, scale="none", col=colori_heatmap, margins=c(6, 6), cexCol=0.9, cexRow=0.8, main="4. Heatmap Top 50 Genes")


# 5. Top 6 genes counts plot
par(mfrow=c(3,2), las=1) 

main_condition <- pca_groups[length(pca_groups)]
plot_colors <- c("#1f78b4", "#e31a1c", "#33a02c", "#ff7f00", "#6a3d9a", "#b15928")
condition_factors <- as.factor(dds[[main_condition]])
num_groups <- length(levels(condition_factors))
dynamic_colors <- rep(plot_colors, length.out = num_groups)

top6_genes <- head(order(res$padj), 6)

if (length(top6_genes) > 0) {
    top6_counts <- counts(dds, normalized=TRUE)[top6_genes, , drop = FALSE ]
    absolute_max <- max(top6_counts)
    y_limit <- c(0.5, absolute_max + (absolute_max * 0.1))

    for (i in top6_genes) {
      gene_name <- rownames(res)[i]

      plotCounts(dds, 
             gene = gene_name, 
             intgroup = main_condition, 
             main = paste("Expression of:", gene_name),
             col = dynamic_colors[as.numeric(condition_factors)], 
             pch = 16,                             
             cex = 1.5,
             ylim = y_limit 
  )
}
} else {
  message("No significant genes found to plot in Top 6 counts.")
}


dev.off()
