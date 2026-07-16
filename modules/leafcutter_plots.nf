process LEAFCUTTER_PLOT {
    tag "Plotting LeafCutter"
    label 'process_low'

    container 'rocker/geospatial:4.3.1'

    input:
    path leafcutter_results

    output:
    path "*_mqc.png", emit: multiqc_png
    path "*.pdf"    , emit: plots

    script:
    """
    cat << 'EOF' > plot_leafcutter.R
    library(ggplot2)

    sig_file <- list.files(pattern = "cluster_significance.txt", recursive = TRUE, full.names = TRUE)[1]
    
    if(!is.na(sig_file)) {
        sig_data <- read.delim(sig_file, header=TRUE, stringsAsFactors=FALSE)
        
        # Filtra i dati non analizzabili
        sig_data <- sig_data[!is.na(sig_data\$p.adjust), ]
        
        if(nrow(sig_data) > 0) {
            
            sig_data\$Significance <- ifelse(sig_data\$p.adjust < 0.05, "Significant (FDR < 0.05)", "Not Significant")
            
            n_sig <- sum(sig_data\$Significance == "Significant (FDR < 0.05)")
            
            p <- ggplot(sig_data, aes(x=p.adjust, fill=Significance)) +
                geom_histogram(breaks=seq(0, 1, by=0.05), color="black", alpha=0.8) +
                scale_fill_manual(values=c("Not Significant" = "grey70", "Significant (FDR < 0.05)" = "firebrick")) +
                theme_minimal(base_size = 14) +
                labs(title=paste0("LeafCutter: Adjusted P-value Distribution\\n(Significant clusters: ", n_sig, ")"),
                     x="FDR (Adjusted P-value)", 
                     y="Number of Clusters") +
                theme(legend.position="bottom")
            
            ggsave("LeafCutter_Pvalue_Distribution.pdf", plot=p, width=8, height=6)
            ggsave("LeafCutter_Pvalue_Distribution_mqc.png", plot=p, width=10, height=6, dpi=300)
        }
    }
    EOF
    
    Rscript plot_leafcutter.R
    """
}
