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
    # 1. Generazione file gruppi
    python -c "
import csv, glob
bams = glob.glob('*.bam')
with open('${samplesheet}', 'r') as f, open('groups_file.txt', 'w') as out:
    reader = csv.DictReader(f, skipinitialspace=True)
    for row in reader:
        sample = row['sample']
        cond = row['${params.design}']
        for b in bams:
            if b.startswith(sample) and not b[len(sample):len(sample)+1].isdigit():
                prefix = b.replace('.bam', '')
                out.write('{}\\t{}\\n'.format(prefix, cond))
"

    mkdir -p leafcutter_out

    # 2. Localizzazione DINAMICA (cerchiamo in tutto il filesystem)
    # Rimuoviamo il fallback sbagliato, se non trova il file deve andare in errore 1
    BAM2JUNC=\$(find / -name "bam2junc.sh" -type f 2>/dev/null | head -n 1)
    CLUSTER_PY=\$(find / -name "leafcutter_cluster.py" -type f 2>/dev/null | head -n 1)
    DS_R=\$(find / -name "leafcutter_ds.R" -type f 2>/dev/null | head -n 1)

    # Verifica immediata: se uno è vuoto, ti dice chiaramente quale manca
    if [ -z "\$BAM2JUNC" ]; then echo "Errore critico: bam2junc.sh non trovato nel container"; exit 1; fi
    if [ -z "\$CLUSTER_PY" ]; then echo "Errore critico: leafcutter_cluster.py non trovato nel container"; exit 1; fi
    if [ -z "\$DS_R" ]; then echo "Errore critico: leafcutter_ds.R non trovato nel container"; exit 1; fi

    # 3. Estrazione giunzioni
    touch juncfiles.txt
    for bam in *.bam; do
        prefix=\${bam%.bam}
        echo "Estraendo giunzioni da \$bam..."
        sh "\$BAM2JUNC" "\$bam" "\${prefix}.junc"
        echo "\${prefix}.junc" >> juncfiles.txt
    done

    # 4. Clustering
    python "\$CLUSTER_PY" -j juncfiles.txt -m 50 -o leafcutter_out/fornax -l 500000

    # 5. Differential Splicing
    Rscript "\$DS_R" \
        --num_threads ${task.cpus} \
        -i 2 \
        -m 5 \
        leafcutter_out/fornax_perind_numers.counts.gz \
        groups_file.txt \
        -o leafcutter_out/fornax_ds
    """
}
