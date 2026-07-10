process RMATS {
    tag "Alternative Splicing"
    label 'process_high'

    container 'xinglab/rmats:v4.3.0'

    input:
    path bams           
    path samplesheet    
    path gtf

    output:
    path "rmats_out/*", emit: splicing_results

script:
    """
    cat << 'EOF' > plot_splicing.R
    library(ggplot2)
    prefix <- "${label}"

    # Cerca il file ovunque si trovi (inclusa la sottocartella rmats_out)
    se_file <- list.files(pattern = "SE.MATS.JC.txt", recursive = TRUE, full.names = TRUE)

    # 1. VOLCANO PLOT
    if(length(se_file) > 0) {
        # Usa se_file[1] che contiene il percorso completo (es: "rmats_out/SE.MATS.JC.txt")
        data <- read.delim(se_file[1], header=TRUE, stringsAsFactors=FALSE)
        data <- data[!is.na(data\$FDR) & !is.na(data\$IncLevelDifference), ]
        data\$Significance <- ifelse(data\$FDR < 0.05 & abs(data\$IncLevelDifference) > 0.1, "Significant", "Not Significant")

        p <- ggplot(data, aes(x=IncLevelDifference, y=-log10(FDR), color=Significance)) +
              geom_point(alpha=0.7, size=1.5) +
              scale_color_manual(values=c("Not Significant" = "grey80", "Significant" = "firebrick")) +
              geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", alpha = 0.6) +
              geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed", color = "black", alpha = 0.6) +
              theme_minimal(base_size = 14) +
              labs(title=paste0(prefix, ": Splicing Volcano Plot (Skipped Exons)"),
                   x="Inclusion Level Difference (delta PSI)",
                   y="-log10(FDR)") +
              theme(legend.position="bottom")
              
        ggsave(paste0(prefix, "_Volcano_SE.pdf"), plot=p, width=8, height=6)
        ggsave(paste0(prefix, "_Volcano_SE_mqc.png"), plot=p, width=14, height=8, dpi=300)
    }

    # 2. BAR PLOT
    events <- c("SE"="SE.MATS.JC.txt", "RI"="RI.MATS.JC.txt", "MXE"="MXE.MATS.JC.txt", "A5SS"="A5SS.MATS.JC.txt", "A3SS"="A3SS.MATS.JC.txt")
    counts <- data.frame(Event=character(), Significant=numeric(), stringsAsFactors=FALSE)

    for (ev in names(events)) {
        f_name <- events[ev]
        # Cerchiamo anche questi file ricorsivamente
        f_path <- list.files(pattern = f_name, recursive = TRUE, full.names = TRUE)
        
        if(length(f_path) > 0) {
            d <- read.delim(f_path[1], header=TRUE, stringsAsFactors=FALSE)
            sig_count <- sum(d\$FDR < 0.05 & abs(d\$IncLevelDifference) > 0.1, na.rm=TRUE)
            counts <- rbind(counts, data.frame(Event=ev, Significant=sig_count))
        }
    }
    
    if(nrow(counts) > 0) {
        p2 <- ggplot(counts, aes(x=Event, y=Significant, fill=Event)) +
              geom_bar(stat="identity", color="black", alpha=0.8) +
              theme_minimal(base_size = 14) +
              labs(title=paste0(prefix, ": Significant Alternative Splicing Events"),
                   x="Event Type", y="Number of Significant Events") +
              theme(legend.position="none")

        ggsave(paste0(prefix, "_Summary_BarChart.pdf"), plot=p2, width=8, height=6)
        ggsave(paste0(prefix, "_Summary_BarChart_mqc.png"), plot=p2, width=14, height=8, dpi=300)
    }
    EOF
    Rscript plot_splicing.R
    """
}
