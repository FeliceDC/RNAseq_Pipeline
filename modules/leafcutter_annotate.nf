process LEAFCUTTER_ANNOTATE {
    tag "Annotating LeafCutter"
    label 'process_low'
    
    // Usiamo lo stesso container così non devi scaricare nulla di nuovo
    container 'lifebitai/leafcutter:latest'

    input:
    path leafcutter_results
    path gtf

    output:
    path "Annotated_LeafCutter_Significance.tsv", emit: annotated_results

    script:
    """
    # Usiamo cat con 'EOF' per blindare lo script dagli errori di sintassi
    cat << 'EOF' > annotate_clusters.py
import sys

gtf_file = sys.argv[1]
sig_file = sys.argv[2]
effect_file = sys.argv[3]

print("1. Parsing GTF per estrarre le coordinate dei geni...")
genes = {}
with open(gtf_file, 'r') as f:
    for line in f:
        if line.startswith('#'): continue
        parts = line.split('\\t')
        # Prendiamo solo le righe che definiscono i confini di un gene intero
        if len(parts) >= 9 and parts[2] == 'gene':
            chrom = parts[0]
            start = int(parts[3])
            end = int(parts[4])
            
            # Estrazione del nome del gene
            attr = parts[8]
            name = "Unknown"
            if 'gene_name "' in attr:
                name = attr.split('gene_name "')[1].split('"')[0]
            elif 'gene_name=' in attr:
                name = attr.split('gene_name=')[1].split(';')[0]
                
            if chrom not in genes:
                genes[chrom] = []
            genes[chrom].append((start, end, name))

print("2. Associazione delle coordinate di LeafCutter ai Geni...")
cluster_to_gene = {}
with open(effect_file, 'r') as f:
    header = f.readline()
    for line in f:
        intron = line.split('\\t')[0] # Formato: chr1:1000:2000:clu_1_+
        parts = intron.split(':')
        if len(parts) < 4: continue
        
        chrom = parts[0]
        start = int(parts[1])
        end = int(parts[2])
        
        # Ricostruiamo l'ID del cluster per incrociarlo con il file significance
        clu_id_rest = ":".join(parts[3:])
        clu_id = "{}:{}".format(chrom, clu_id_rest) 
        
        if clu_id not in cluster_to_gene:
            matched_genes = set()
            if chrom in genes:
                for g_start, g_end, g_name in genes[chrom]:
                    # La logica matematica dell'Overlap genomico
                    if start <= g_end and end >= g_start:
                        matched_genes.add(g_name)
            
            if matched_genes:
                cluster_to_gene[clu_id] = ",".join(matched_genes)
            else:
                cluster_to_gene[clu_id] = "Intergenic" # Se non cade in nessun gene noto

print("3. Scrittura del referto annotato...")
with open(sig_file, 'r') as f, open("Annotated_LeafCutter_Significance.tsv", 'w') as out:
    header = f.readline().strip()
    out.write(header + "\\tGene_Name\\n")
    for line in f:
        parts = line.strip().split('\\t')
        clu_id = parts[0]
        g = cluster_to_gene.get(clu_id, "Unknown")
        out.write(line.strip() + "\\t" + g + "\\n")
EOF

    # Troviamo automaticamente i file di input generati dal processo precedente
    SIG_FILE=\$(find ${leafcutter_results} -name "*cluster_significance.txt" | head -n 1)
    EFF_FILE=\$(find ${leafcutter_results} -name "*effect_sizes.txt" | head -n 1)

    # Lanciamo l'annotazione
    python annotate_clusters.py ${gtf} \$SIG_FILE \$EFF_FILE
    """
}
