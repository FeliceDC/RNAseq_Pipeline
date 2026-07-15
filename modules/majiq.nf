process MAJIQ {
    tag "MAJIQ (Local Splicing Variations)"
    label 'process_high'

container 'quay.io/biocontainers/majiq:2.5.3--py310h30b5030_0'

    input:
    path bams
    path samplesheet
    path gtf

    output:
    path "majiq_out/*", emit: majiq_results
    path "majiq_build_out/*", emit: majiq_build_results, optional: true

    script:
    """
    python -c "
import csv, os, glob

bams = glob.glob('*.bam')
groups = {}

with open('${samplesheet}', 'r') as f:
    reader = csv.DictReader(f, skipinitialspace=True)
    for row in reader:
        sample = row['sample']
        cond = row['${params.design}']
        if cond not in groups:
            groups[cond] = []
            
        for b in bams:
            if b.startswith(sample) and not b[len(sample):len(sample)+1].isdigit():
                groups[cond].append(b.replace('.bam', ''))

# Scrittura del file INI
with open('majiq_config.ini', 'w') as config:
    config.write('[info]\\n')
    config.write('readlen=${params.splicing_read_length}\\n')
    config.write('bamdirs=.\\n')
    config.write('genome=hg38\\n') 
    config.write('[experiments]\\n')
    
    for cond, samples_list in groups.items():
        config.write(cond + '=' + ','.join(samples_list) + '\\n')
"

    mkdir -p majiq_build_out
    mkdir -p majiq_out

    # 2. MAJIQ Build
    majiq build \\
        ${gtf} \\
        -c majiq_config.ini \\
        -j ${task.cpus} \\
        -o majiq_build_out

    # Estrazione nomi condizioni in variabili Bash - CORRETTO L'ESCAPE
    COND1=\$(awk -F'=' '/^\\[experiments\\]/{flag=1; next} flag && NF {print \$1; exit}' majiq_config.ini)
    COND2=\$(awk -F'=' '/^\\[experiments\\]/{flag=1; next} flag && NF {NR++; if(NR==2) print \$1; exit}' majiq_config.ini)

    # 3. MAJIQ DeltaPSI
    majiq deltapsi \\
        -grp1 majiq_build_out/\${COND1}.majiq \\
        -grp2 majiq_build_out/\${COND2}.majiq \\
        -j ${task.cpus} \\
        -o majiq_out \\
        -n \${COND1} \${COND2}
    """
}
