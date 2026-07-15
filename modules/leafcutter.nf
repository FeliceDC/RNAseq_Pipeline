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
    # 1. Generazione del file dei gruppi (groups_file.txt)
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
                # Usiamo .format() con DOPPIO backslash per l'escape
                out.write('{}\\t{}\\n'.format(prefix, cond))
"

    touch juncfiles.txt
    mkdir -p leafcutter_out

    # --- RICERCA DINAMICA DEGLI SCRIPT ---
    # Troviamo esattamente dove il container ha nascosto questi file
    BAM2JUNC=\$(find / -name "bam2junc.sh" -type f 2>/dev/null | head -n 1)
    CLUSTER_PY=\$(find / -name "leafcutter_cluster.py" -type f 2>/dev/null | head -n 1)
    DS_R=\$(find / -name "leafcutter_ds.R" -type f 2>/dev/null | head -n 1)

    # 2. Estrazione delle giunzioni
    for bam in *.bam; do
        prefix=\${bam%.bam}
        echo "Estraendo giunzioni da \$bam..."
        
        # Usiamo la variabile trovata da Linux
        sh \$BAM2JUNC \$bam \${prefix}.junc
        
        echo \${prefix}.junc >> juncfiles.txt
    done

    # 3. Clustering delle giunzioni
    python \$CLUSTER_PY \\
        -j juncfiles.txt \\
        -m 50 \\
        -o leafcutter_out/fornax \\
        -l 500000

    # 4. Differential Splicing
    Rscript \$DS_R \\
        --num_threads ${task.cpus} \\
        leafcutter_out/fornax_perind_numers.counts.gz \\
        groups_file.txt \\
        -o leafcutter_out/fornax_ds
    """
