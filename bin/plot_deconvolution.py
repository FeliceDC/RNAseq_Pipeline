#!/usr/bin/env python3
import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

input_excel = sys.argv[1]

print(f"Lettura dei risultati da {input_excel}...")
df = pd.read_excel(input_excel, index_col=0)
df_t = df.T

#PLOT 1
print("Generazione dello Stacked Barplot con colori personalizzati...")
cell_colors = {
    # LINFOCITI B (Tonalità di Blu)
    'Bnaive': '#08306b', 'Breg': '#08519c', 'FOB': '#2171b5', 'GC_B': '#4292c6',
    'MZB': '#6baed6', 'memoryB': '#9ecae1', 'plasma': '#c6dbef', 'plasmablast': '#deebf7', 'exhaustedB': '#17408B',

    # LINFOCITI T CD4+ (Tonalità di Verde)
    'CD4Tcm': '#00441b', 'CD4Tem': '#006d2c', 'CD4Temra': '#238b45', 'CD4Tnaive': '#41ab5d',
    'CD4Trm': '#74c476', 'Tfh': '#a1d99b', 'Th1': '#c7e9c0', 'Th1/Th17': '#e5f5e0',
    'Th17': '#f7fcf5', 'Th2': '#3F704D', 'Tr1': '#8F9779', 'Treg': '#4F7942',

    # LINFOCITI T CD8+ E NK (Tonalità di Viola/Lilla)
    'CD8Tcm': '#3f007d', 'CD8Tem': '#54278f', 'CD8Temra': '#6a51a3', 'CD8Tnaive': '#807dba',
    'CD8Trm': '#9e9ac8', 'Tc': '#bcbddc', 'exhausted_T': '#dadaeb', 'cytotoxicNK': '#e0ecf4', 'regulatoryNK': '#bfd3e6',

    # MACROFAGI E MONOCITI (Tonalità di Arancione/Marrone)
    'CMonocyte': '#7f2704', 'IMonocyte': '#a63603', 'NMonocyte': '#d94801', 'TAM': '#f16913',
    'M0': '#fd8d3c', 'M1': '#fdae6b', 'M2': '#fdd0a2', 'MDSCs': '#feedde', 'Langerhans': '#8B4513',

    # CELLULE DENDRITICHE E GRANULOCITI (Tonalità di Rosso/Giallo)
    'cDC1': '#67000d', 'cDC2': '#a50f15', 'monoDC': '#cb181d', 'pDC': '#ef3b2c',
    'basophils': '#fb6a4a', 'eosinophils': '#fc9272', 'neutrophils': '#fcbba1', 'mast_cell': '#fee0d2',

    # ILC E ALTRI (Tonalità Grigio/Neutro)
    'ILC1': '#525252', 'ILC2': '#737373', 'ILC3': '#969696', 'MAIT': '#bdbdbd', 'NKT': '#d9d9d9', 'gdT': '#f0f0f0'
}
ordered_colors = [cell_colors.get(cell, '#cccccc') for cell in df_t.columns]
fig, ax = plt.subplots(figsize=(14, 8))
df_t.plot(kind='bar', stacked=True, ax=ax, color=ordered_colors, width=0.85)
plt.legend(title='Cell Types (ImmuCellAI v2)', bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=7, ncol=3)
plt.title("Tumor Microenvironment (TME) Cellular Composition", fontsize=16, fontweight='bold', pad=20)
plt.ylabel("Relative Abundance (Fraction)", fontsize=12, fontweight='bold')
plt.xlabel("Samples", fontsize=12, fontweight='bold')
plt.xticks(rotation=45, ha='right')
plt.tight_layout()
plt.savefig("TME_Stacked_Barplot.pdf", dpi=300, bbox_inches='tight')
plt.savefig("Stacked_Barplot_mqc.png", dpi=300, bbox_inches='tight')
plt.close()

#PLOT 2
sns.set(font_scale=0.7)
g = sns.clustermap(df, cmap="YlGnBu", figsize=(12, 15), standard_scale=1, method='ward')
g.ax_heatmap.set_title("Clustered TME Profile (Standardized)", fontsize=15, fontweight='bold', pad=40)
plt.savefig("TME_Clustered_Heatmap.pdf", dpi=300, bbox_inches='tight')
plt.close()

#PLOT3
print("Generazione Correlation Matrix...")
# Calcola la correlazione tra i diversi tipi cellulari
corr = df_t.corr()
plt.figure(figsize=(12, 10))
# Plot con cerchi/quadrati per vedere quali cellule correlano tra loro
sns.heatmap(corr, cmap="RdBu_r", center=0, square=True, linewidths=.5, cbar_kws={"shrink": .5})
plt.title("Cell-Type Correlation Matrix", fontsize=15, fontweight='bold')
plt.tight_layout()
plt.savefig("TME_Cell_Correlation.pdf", dpi=300, bbox_inches='tight')
plt.close()

#PLOT4
print("Generazione Boxplot...")
# Trasformiamo i dati per Seaborn
df_melted = df.reset_index().melt(id_vars='index')
df_melted.columns = ['CellType', 'Sample', 'Fraction']
plt.figure(figsize=(15, 8))
sns.boxplot(data=df_melted, x='CellType', y='Fraction', palette='viridis')
plt.xticks(rotation=90, fontsize=8)
plt.title("Cell Type Distribution across all samples", fontsize=15, fontweight='bold')
plt.tight_layout()
plt.savefig("TME_Abundance_Boxplot.pdf", dpi=300, bbox_inches='tight')
plt.close()
