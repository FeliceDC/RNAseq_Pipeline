process LEAFCUTTER_ANNOTATE {
    tag "Annotating LeafCutter"
    label 'process_low'
    
    container 'lifebitai/leafcutter:latest'

    input:
    path leafcutter_results
    path gtf

    output:
    path "Annotated_LeafCutter_Significance.tsv", emit: annotated_results

    script:
    """
    cat << 'EOF' > annotate_clusters.py
import sys

gtf_file = sys.argv[1]
sig_file = sys.argv[2]
effect_file = sys.argv[3]

print("1. Parsing GTF...")
genes = {}
with open(gtf_file, 'r') as f:
    for line in f:
        if line.startswith('#'): continue
        parts = line.split('\\t')
        if len(parts) >= 9 and parts[2] == 'gene':
            chrom = parts[0]
            start = int(parts[3])
            end = int(parts[4])
            
            attr = parts[8]
            name = "Unknown"
            if 'gene_name "' in attr:
                name = attr.split('gene_name "')[1].split('"')[0]
            elif 'gene_name=' in attr:
                name = attr.split('gene_name=')[1].split(';')[0]
                
            if chrom not in genes:
                genes[chrom] = []
            genes[chrom].append((start, end, name))

print("2. Associazione e pulizia del bug 'chrchr'...")
cluster_to_gene = {}
with open(effect_file, 'r') as f:
    header = f.readline()
    for line in f:
        intron = line.split('\\t')[0]
        parts = intron.split(':')
        if len(parts) < 4: continue
        
        # CORREZIONE DEL BUG: sostituiamo chrchr con chr
        clean_chrom = parts[0].replace("chrchr", "chr")
        start = int(parts[1])
        end = int(parts[2])
        
        clu_id_rest = ":".join(parts[3:-1]) if len(parts) >= 5 else parts[3]
        unified_clu_id = "{}:{}".format(clean_chrom, clu_id_rest) 
        
        if unified_clu_id not in cluster_to_gene:
            matched_genes = set()
            if clean_chrom in genes:
                for g_start, g_end, g_name in genes[clean_chrom]:
                    if start <= g_end and end >= g_start:
                        matched_genes.add(g_name)
            
            if matched_genes:
                cluster_to_gene[unified_clu_id] = ",".join(matched_genes)
            else:
                cluster_to_gene[unified_clu_id] = "Intergenic"

print("3. Scrittura del referto annotato...")
with open(sig_file, 'r') as f, open("Annotated_LeafCutter_Significance.tsv", 'w') as out:
    header = f.readline().strip()
    out.write(header + "\\tGene_Name\\n")
    for line in f:
        parts = line.strip().split('\\t')
        
        # Puliamo anche l'ID nel file finale per avere risultati eleganti
        clean_cluster = parts[0].replace("chrchr", "chr")
        parts[0] = clean_cluster
        
        g = cluster_to_gene.get(clean_cluster, "Unknown")
        out.write("\\t".join(parts) + "\\t" + g + "\\n")
EOF

    SIG_FILE=\$(find ${leafcutter_results} -name "*cluster_significance.txt" | head -n 1)
    EFF_FILE=\$(find ${leafcutter_results} -name "*effect_sizes.txt" | head -n 1)

    python annotate_clusters.py ${gtf} \$SIG_FILE \$EFF_FILE
    """
}
