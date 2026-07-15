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
    # 1. Creiamo lo script Python separatamente per evitare QUALSIASI errore di sintassi Groovy
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

    # 2. Eseguiamo lo script Python passando i parametri di Nextflow in modo sicuro
    python generate_groups.py ${samplesheet} ${params.design}

    touch juncfiles.txt
    mkdir -p leafcutter_out

    # 3. Ricerca rapida degli script originali nel container (nella cartella /opt)
    BAM2JUNC=\$(find /opt -name "bam2junc.sh" -type f 2>/dev/null | head -n 1)
    CLUSTER_PY=\$(find /opt -name "leafcutter_cluster.py" -type f 2>/dev/null | head -n 1)
    DS_R=\$(find /opt -name "leafcutter_ds.R" -type f 2>/dev/null | head -n 1)

    # 4. Estrazione giunzioni
    for bam in *.bam; do
        prefix=\${bam%.bam}
        echo "Estraendo giunzioni da \$bam..."
        
        sh \$BAM2JUNC \$bam \${prefix}.junc
        
        echo \${prefix}.junc >> juncfiles.txt
    done

    # 5. Clustering (messo su una riga singola per evitare errori di escape)
    python \$CLUSTER_PY -j juncfiles.txt -m 50 -o leafcutter_out/fornax -l 500000

    # 6. Differential Splicing
    Rscript \$DS_R --num_threads ${task.cpus} leafcutter_out/fornax_perind_numers.counts.gz groups_file.txt -o leafcutter_out/fornax_ds
    """
}
