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
                # Usiamo .format() con DOPPIO backslash per l'escape corretto in Nextflow
                out.write('{}\\t{}\\n'.format(prefix, cond))
"

    touch juncfiles.txt
    mkdir -p leafcutter_out

    # 2. Estrazione delle giunzioni usando lo script nativo di LeafCutter
    for bam in *.bam; do
        # ESCAPE AGGIUNTO QUI: \${bam%.bam}
        prefix=\${bam%.bam}
        
        # ESCAPE AGGIUNTO QUI: \$bam
        echo "Estraendo giunzioni da \$bam..."
        
        # ESCAPE AGGIUNTO QUI: \$bam e \${prefix}
        sh /opt/software/leafcutter/scripts/bam2junc.sh \$bam \${prefix}.junc
        
        # ESCAPE AGGIUNTO QUI: \${prefix}
        echo \${prefix}.junc >> juncfiles.txt
    done

    # 3. Clustering delle giunzioni
    python /opt/software/leafcutter/clustering/leafcutter_cluster.py \\
        -j juncfiles.txt \\
        -m 50 \\
        -o leafcutter_out/fornax \\
        -l 500000

    # 4. Differential Splicing
    /opt/software/leafcutter/scripts/leafcutter_ds.R \\
        --num_threads ${task.cpus} \\
        leafcutter_out/fornax_perind_numers.counts.gz \\
        groups_file.txt \\
        -o leafcutter_out/fornax_ds
    """
}
