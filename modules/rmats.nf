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
    python -c "
    import csv, os, glob
    bams = glob.glob('*.bam')
    groups = {}
    with open('${samplesheet}', 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            sample = row['sample']
            cond = row['${params.design}']
            if cond not in groups:
                groups[cond] = []
            for b in bams:
                if b.startswith(sample + "."):
    
    conds = list(groups.keys())
    with open('b1.txt', 'w') as f1:
        f1.write(','.join(groups[conds[0]]))
    with open('b2.txt', 'w') as f2:
        f2.write(','.join(groups[conds[1]]))
    "

    python /rmats/rmats.py \\
        --b1 b1.txt \\
        --b2 b2.txt \\
        --gtf ${gtf} \\
        -t ${params.single_end ? 'single' : 'paired'} \ \\
       --readLength ${params.splicing_read_length} \ \\
        --nthread ${task.cpus} \\
        --od rmats_out \\
        --tmp rmats_tmp
    """
}
