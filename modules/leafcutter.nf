process LEAFCUTTER {
    tag "LeafCutter (Annotation-free Splicing)"
    label 'process_high'
    container 'lifebitai/leafcutter:latest'

    input:
    path bams
    path samplesheet

    output:
    path "leafcutter_out/*", emit: leafcutter_results

    script:
    """
    # 1. Creiamo lo script Python separatamente per BLOCCARE gli errori di Groovy
    cat << 'EOF' > generate_groups.py
import csv, glob, sys
    
samplesheet_file = sys.argv[1]
design_col = sys.argv[2]
    
bams = glob.glob('*.bam')
    
with open(samplesheet_file, 'r') as f, open('groups_file.txt', 'w') as out:
    reader = csv.DictReader(f, skipinitialspace=True)
    for row in reader:
        sample = row['sample']
        cond = row[design_col]
        for b in bams:
            if b.startswith(sample) and not b[len(sample):len(sample)+1].isdigit():
                prefix = b.replace('.bam', '')
                out.write("{}\\t{}\\n".format(prefix, cond))
EOF

    python generate_groups.py ${samplesheet} ${params.design}

    mkdir -p leafcutter_out

    # 3. Localizzazione DINAMICA sicura degli script nel container
    BAM2JUNC=\$(find / -name "bam2junc.sh" -type f 2>/dev/null | head -n 1)
    CLUSTER_PY=\$(find / -name "leafcutter_cluster.py" -type f 2>/dev/null | head -n 1)
    DS_R=\$(find / -name "leafcutter_ds.R" -type f 2>/dev/null | head -n 1)

    if [ -z "\$BAM2JUNC" ]; then echo "Errore critico: bam2junc.sh non trovato"; exit 1; fi
    if [ -z "\$CLUSTER_PY" ]; then echo "Errore critico: leafcutter_cluster.py non trovato"; exit 1; fi
    if [ -z "\$DS_R" ]; then echo "Errore critico: leafcutter_ds.R non trovato"; exit 1; fi

    # 4. Estrazione giunzioni
    touch juncfiles.txt
    for bam in *.bam; do
        prefix=\${bam%.bam}
        echo "Estraendo giunzioni da \$bam..."
        sh "\$BAM2JUNC" "\$bam" "\${prefix}.junc"
        echo "\${prefix}.junc" >> juncfiles.txt
    done

    # 5. Clustering
    python "\$CLUSTER_PY" -j juncfiles.txt -m 50 -o leafcutter_out/leafcutter -l 500000

    # 6. Differential Splicing (con i parametri corretti per piccoli gruppi)
    Rscript "\$DS_R" \\
        --num_threads ${task.cpus} \\
        -i 2 \\
        -g 2 \\
        -o leafcutter_out/leafcutter_ds \\
        leafcutter_out/leafcutter_perind_numers.counts.gz \\
        groups_file.txt
    """
}
