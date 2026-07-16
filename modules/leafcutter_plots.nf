process LEAFCUTTER_PLOT {
    tag "Plotting LeafCutter"
    label 'process_low'
    
    container 'rocker/geospatial:4.3.1'

    input:
    // Ora in input prende direttamente il file annotato!
    path annotated_sig_file

    output:
    path "*_mqc.png", emit: multiqc_png
    path "*.pdf"    , emit: plots

    script:
    """
    cat << 'EOF' > plot_leafcutter.R
    library(ggplot2)

    # Legge il file annotato passato da Nextflow
    sig_data <- read.delim("${annotated_sig_file}", header=TRUE, stringsAsFactors=FALSE)
    
    # Filtra i dati validi
    sig_data <- sig_data[!is.na(sig_data\$p.adjust), ]
    
    if(nrow(sig_data) > 0) {
        # 1. PREPARAZIONE DATI
        sig_data\$Significance <- ifelse(sig_data\$p.adjust < 0.05, "Significant", "Not Significant")
        sig_data\$log10_fdr <- -log10(sig_data\$p.adjust + 1e-300) # Evita log(0)
        
        n_sig <- sum(sig_data\$Significance == "Significant")
        
        # 2. ISTOGRAMMA P-VALUE (Come prima)
        p_hist <- ggplot(sig_data, aes(x=p.adjust, fill=Significance)) +
            geom_histogram(breaks=seq(0, 1, by=0.05), color="black", alpha=0.8) +
            scale_fill_manual(values=c("Not Significant" = "grey70", "Significant" = "firebrick")) +
            theme_minimal(base_size = 14) +
            labs(title=paste0("LeafCutter: Adjusted P-value Distribution\\n(Significant clusters: ", n_sig, ")"),
                 x="FDR (Adjusted P-value)", y="Number of Clusters") +
            theme(legend.position="bottom")
        
        ggsave("LeafCutter_Pvalue_Histogram.pdf", plot=p_hist, width=8, height=6)
        
        # 3. VOLCANO PLOT (La magia dell'annotazione)
        # Troviamo i Top 15 geni più significativi da etichettare (escludendo Intergenic/Unknown)
        sig_genes_only <- sig_data[sig_data\$Significance == "Significant" & 
                                   sig_data\$Gene_Name != "Intergenic" & 
                                   sig_data\$Gene_Name != "Unknown", ]
        top_genes <- head(sig_genes_only[order(sig_genes_only\$p.adjust), ], 15)
        
        p_volc <- ggplot(sig_data, aes(x=loglr, y=log10_fdr, color=Significance)) +
            geom_point(alpha=0.6, size=1.5) +
            scale_color_manual(values=c("Not Significant" = "grey60", "Significant" = "dodgerblue4")) +
            theme_minimal(base_size = 14) +
            labs(title="LeafCutter: Differential Splicing Volcano Plot",
                 x="Log-Likelihood Ratio (Magnitude of Change)",
                 y="-log10(FDR)") +
            theme(legend.position="bottom")
            
        # Aggiungiamo le etichette con i Nomi dei Geni per i top hits
        if(nrow(top_genes) > 0) {
            p_volc <- p_volc + geom_text(data=top_genes, aes(label=Gene_Name), 
                                         vjust=-1, color="black", size=3, check_overlap=TRUE)
        }
        
        ggsave("LeafCutter_Volcano_Plot.pdf", plot=p_volc, width=8, height=6)
        
        # Salviamo la versione per MultiQC (usiamo il Volcano che è più d'impatto)
        ggsave("LeafCutter_Volcano_mqc.png", plot=p_volc, width=10, height=6, dpi=300)
    }
EOF
    
    Rscript plot_leafcutter.R
    """
}
