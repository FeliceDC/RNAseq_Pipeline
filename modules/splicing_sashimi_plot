process SASHIMI_PLOT {
    tag "Sashimi Plot"
    label 'process_medium'

    container 'xinglab/rmats2sashimiplot:v3.0.0'

    input:
    path bams
    path samplesheet
    path rmats_files 

    output:
    path "sashimi_out/Sashimi_plot/*.pdf", emit: plots

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
            cond = row['condition']
            if cond not in groups:
                groups[cond] = []
            for b in bams:
                if sample in b:
                    groups[cond].append(os.path.abspath(b))
    
    conds = list(groups.keys())
    with open('b1.txt', 'w') as f1:
        f1.write(','.join(groups[conds[0]]))
    with open('b2.txt', 'w') as f2:
        f2.write(','.join(groups[conds[1]]))
    "

    head -n 6 SE.MATS.JC.txt > top5_SE.txt

    rmats2sashimiplot \\
        --b1 b1.txt \\
        --b2 b2.txt \\
        --event-type SE \\
        -e top5_SE.txt \\
        --l1 Control --l2 Treated \\
        --exon_s 1 --intron_s 5 \\
        -o sashimi_out
    """
}
