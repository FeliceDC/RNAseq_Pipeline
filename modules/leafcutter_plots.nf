process LEAFCUTTER_PLOT {
    tag "Plotting LeafCutter"
    label 'process_low'
    
    container 'rocker/geospatial:4.3.1'

    input:
    path input_files // Accetta qualsiasi cosa (singolo file o cartella)

    output:
    path "*_mqc.png", emit: multiqc_png
    path "*.pdf"    , emit: plots

    script:
    """
    cat << 'EOF' > plot_leafcutter.R
    library(ggplot2)

    # 1. Ricerca intelligente del file di input
    # Cerca prima il file annotato, altrimenti ripiega su quello normale
    files <- list.files(pattern = "Annotated_LeafCutter_Significance.tsv", full.names = TRUE)
    if(length(files) == 0) {
        files <- list.files(pattern = "cluster_significance.txt", full.names = TRUE)
    }
    
    if(length(files) == 0) {
        stop("Errore: Nessun file di significatività trovato per il plot!")
    }
    
    sig_file <- files[1]
    cat("Leggendo il file:", sig_file, "\\n")
    
    # Legge il file trovato
    sig_data <- read.delim(sig_file, header=TRUE, stringsAsFactors=FALSE)
    
    # Filtra i dati validi
    sig_data <- sig_data[!is.na(sig_data\$p.adjust), ]
    
    if(nrow(sig_data) > 0) {
        # 2. PREPARAZIONE DATI
        sig_data\$Significance <- ifelse(sig_data\$p.adjust < 0.05, "Significant", "Not Significant")
        sig_data\$log10_fdr <- -log10(sig_data\$p.adjust + 1e-300) # Evita log(0)
        
        n_sig <- sum(sig_data\$Significance == "Significant")
        
        # 3. ISTOGRAMMA P-VALUE
        p_hist <- ggplot(sig_data, aes(x=p.adjust, fill=Significance)) +
            geom_histogram(breaks=seq(0, 1, by=0.05), color="black", alpha=0.8) +
            scale_fill_manual(values=c("Not Significant" = "grey70", "Significant" = "firebrick")) +
            theme_minimal(base_size = 14) +
            labs(title=paste0("LeafCutter: Adjusted P-value Distribution\\n(Significant clusters: ", n_sig, ")"),
                 x="FDR (Adjusted P-value)", y="Number of Clusters") +
            theme(legend.position="bottom")
        
        ggsave("LeafCutter_Pvalue_Histogram.pdf", plot=p_hist, width=8, height=6)
        
        # 4. VOLCANO PLOT
        # Controlla se la colonna Gene_Name esiste (se è stato annotato)
        has_annotation <- "Gene_Name" %in% colnames(sig_data)
        
        p_volc <- ggplot(sig_data, aes(x=loglr, y=log10_fdr, color=Significance)) +
            geom_point(alpha=0.6, size=1.5) +
            scale_color_manual(values=c("Not Significant" = "grey60", "Significant" = "dodgerblue4")) +
            geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "darkred", alpha = 0.5) +
            geom_vline(xintercept = 2.0, linetype = "dashed", color = "darkred", alpha = 0.5) +
            theme_minimal(base_size = 14) +
            labs(title="LeafCutter: Differential Splicing Volcano Plot",
                 x="Log-Likelihood Ratio (Magnitude of Change)",
                 y="-log10(FDR)") +
            theme(legend.position="bottom")
            
        # Aggiungiamo le etichette solo se i geni sono stati annotati
        if(has_annotation) {
            sig_genes_only <- sig_data[sig_data\$Significance == "Significant" & 
                                       sig_data\$Gene_Name != "Intergenic" & 
                                       sig_data\$Gene_Name != "Unknown", ]
            if(nrow(sig_genes_only) > 0) {
                top_genes <- head(sig_genes_only[order(sig_genes_only\$p.adjust), ], 15)
                p_volc <- p_volc + geom_text(data=top_genes, aes(label=Gene_Name), 
                                             vjust=-1, color="black", size=3, check_overlap=TRUE)
            }
        }
        
        ggsave("LeafCutter_Volcano_Plot.pdf", plot=p_volc, width=8, height=6)
        ggsave("LeafCutter_Volcano_mqc.png", plot=p_volc, width=10, height=6, dpi=300)
    }
EOF
    
    Rscript plot_leafcutter.R
    """
}
