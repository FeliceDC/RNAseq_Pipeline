#!/usr/bin/env Rscript

if (!requireNamespace("enrichR", quietly = TRUE)) {
    install.packages(c("enrichR", "optparse"), repos="https://cloud.r-project.org")
}

if (!requireNamespace("optparse", quietly = TRUE)) {
    install.packages("optparse", repos="http://cran.us.r-project.org")
}

suppressPackageStartupMessages({
    library(enrichR)
    library(optparse)
})

option_list <- list(
    make_option(c("-i", "--input"), type="character", default=NULL,
                help="TXT file of the filtered DESeq2 results"),
    make_option(c("-d", "--databases"), type="character", default="KEGG_2023_Human",
                help="Comma-separated Enrichr databases"),
    make_option(c("-o", "--outdir"), type="character", default=".",
                help="Outut directory")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input)){
    print_help(opt_parser)
    stop("Wait: You must provide an input file using: --input", call.=FALSE)
}

cat("Reading the DEGs gene file...\n")
signif_genes <- read.table(opt$input, header=TRUE, sep="\t", stringsAsFactors=FALSE)

gene_list <- as.character(signif_genes$Gene_Name)
gene_list <- gene_list[!grepl("^ENSG", gene_list)]

if (length(gene_list) == 0) {
    stop("No genes found in the input file")
}

dbs <- unlist(strsplit(opt$databases, ","))
setEnrichrSite("Enrichr") 

cat(paste("Querying Enrichr for", length(gene_list), "genes...\n"))
enriched <- enrichr(gene_list, dbs)

for (db in dbs) {
    write.csv(enriched[[db]], file=file.path(opt$outdir, paste0("Enrichr_", db, ".csv")), row.names=FALSE)
    
    if (nrow(enriched[[db]]) > 0) {
        pdf(file=file.path(opt$outdir, paste0("Barplot_", db, ".pdf")), width=10, height=6)
        print(plotEnrich(enriched[[db]], showTerms = 15, numChar = 50, y = "Count", orderBy = "P.value", title = db))
        dev.off()
    }
}
cat("Pathways analysis completed succesfully\n")
