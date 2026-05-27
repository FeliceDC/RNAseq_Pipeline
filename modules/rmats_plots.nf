process RMATS_PLOT {
    tag "rMATS Volcano Plot"
    label 'process_low'

    container 'rocker/geospatial:4.3.1'

    input:
    path rmats_dir

    output:
    path "*.pdf", emit: plots

    script:
    """
    cat << 'EOF' > plot_rmats.R
    library(ggplot2)

    data <- read.table(file.path("${rmats_dir}", "SE.MATS.JC.txt"), header=TRUE, sep="\t", stringsAsFactors=FALSE)
    
    data <- data[!is.na(data\$FDR) & !is.na(data\$IncLevelDifference), ]
    
    data\$Significance <- ifelse(data\$FDR < 0.05 & abs(data\$IncLevelDifference) > 0.1, "Significant", "Not Significant")

    p <- ggplot(data, aes(x=IncLevelDifference, y=-log10(FDR), color=Significance)) +
         geom_point(alpha=0.7, size=1.5) +
         scale_color_manual(values=c("Not Significant" = "grey80", "Significant" = "firebrick")) +
         theme_minimal(base_size = 14) +
         labs(title="Alternative Splicing Volcano Plot (Skipped Exons)",
              subtitle="Comparison of Inclusion Level vs Significance",
              x="Inclusion Level Difference (delta PSI)",
              y="-log10(FDR)") +
         theme(legend.position="bottom")


    files <- c("SE"="SE.MATS.JC.txt", "RI"="RI.MATS.JC.txt", "MXE"="MXE.MATS.JC.txt", "A5SS"="A5SS.MATS.JC.txt", "A3SS"="A3SS.MATS.JC.txt")
    counts <- data.frame(Event=character(), Significant=numeric(), stringsAsFactors=FALSE)

    for (ev in names(files)) {
        f_path <- file.path("${rmats_dir}", files[ev])
        if(file.exists(f_path)) {
            d <- read.table(f_path, header=TRUE, sep="\t", stringsAsFactors=FALSE)
            # Conta solo quelli con FDR < 0.05 e differenza maggiore del 10%
            sig_count <- sum(d\$FDR < 0.05 & abs(d\$IncLevelDifference) > 0.1, na.rm=TRUE)
            counts <- rbind(counts, data.frame(Event=ev, Significant=sig_count))
        }
    }

    p2 <- ggplot(counts, aes(x=Event, y=Significant, fill=Event)) +
          geom_bar(stat="identity", color="black", alpha=0.8) +
          theme_minimal(base_size = 14) +
          labs(title="Significant Alternative Splicing Events",
               x="Event Type", y="Number of Significant Events") +
          theme(legend.position="none")

    ggsave("rMATS_Summary_BarChart.pdf", plot=p2, width=8, height=6)

    ggsave("rMATS_Volcano_SE.pdf", plot=p, width=8, height=6)
    EOF

    Rscript plot_rmats.R
    """
}
