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
    # 1. Creazione file gruppi (corretto per evitare NameError)
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

    touch juncfiles.txt
    mkdir -p leafcutter_out

    # 2. Localizzazione sicura degli script
    # 'which' cerca nei percorsi standard del sistema, molto più affidabile di 'find'
    BAM2JUNC=\$(which bam2junc.sh || echo '/opt/software/leafcutter/scripts/bam2junc.sh')
    CLUSTER_PY=\$(which leafcutter_cluster.py || echo '/opt/software/leafcutter/clustering/leafcutter_cluster.py')
    DS_R=\$(which leafcutter_ds.R || echo '/opt/software/leafcutter/scripts/leafcutter_ds.R')

    # Controllo di sicurezza: se non troviamo gli script, blocchiamo subito
    if [ ! -f "\$BAM2JUNC" ]; then echo "Errore: bam2junc.sh non trovato!"; exit 1; fi

    # 3. Estrazione giunzioni
    for bam in *.bam; do
        prefix=\${bam%.bam}
        echo "Estraendo giunzioni da \$bam..."
        
        # Esecuzione sicura
        sh \$BAM2JUNC \$bam \${prefix}.junc
        
        echo \${prefix}.junc >> juncfiles.txt
    done

    # 4. Clustering
    python \$CLUSTER_PY -j juncfiles.txt -m 50 -o leafcutter_out/fornax -l 500000

    # 5. Differential Splicing
    Rscript \$DS_R --num_threads ${task.cpus} leafcutter_out/fornax_perind_numers.counts.gz groups_file.txt -o leafcutter_out/fornax_ds
    """
}
