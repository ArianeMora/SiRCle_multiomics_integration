# Imports
import os
from scircm import *  # Note if you have a mac M1 use from sircle import * and you won't be able to do 7,8
import seaborn as sns
from sciutil import *
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams

# Setup file locations and label of the cancer
u = SciUtil()

#  ----------------------------------------------------------------------------------------------
#                  File and naming setup
#  ----------------------------------------------------------------------------------------------
cancer = 'ClearCellRenalCellCarcinoma'

output_dir = 'Output_Data'
fig_dir = 'Output_Figures'
supp_dir = 'Required_Refs'

# Files
DNA_methylation_file = os.path.join('Input_Methylation', f'{cancer}_filtered_DCpG.csv')
RNA_file = os.path.join('Input_RNAseq', f'{cancer}_filtered_DE_RNA.csv')
protein_file = os.path.join('Input_Protein', f'{cancer}_filtered_DA_Protein.csv')

# Names of columns in each of the above files
logFC_RNA_column = "logFC_rna"
padj_RNA_column = "padj_rna"
betadiff_DNA_methylation_column = "beta_diff"
padj_DNA_methylation_column = 'adj.P.Val'
logFC_protein_column = "logFC_protein"
padj_protein_column = "padj_protein"

# NOTE: all of the above files MUST have this column i.e. they must all have the samely named gene ID column
gene_id_column = "gene_name"

# Set cutoffs for RCM
rna_padj_cutoff = 0.05
prot_padj_cutoff = 0.05
meth_padj_cutoff = 0.05

rna_logfc_cutoff = 1.0
prot_logfc_cutoff = 0.5
meth_diff_cutoff = 0.1
bg_type = 'P|(M&R)'

# Give your run a logical name
label = f'PorMandR_{cancer}'
sircle_file_name = f'sircle_{label}.csv'


#  ----------------------------------------------------------------------------------------------
#                   Run RCM
#  ----------------------------------------------------------------------------------------------
# Join the DCpG file with the epic manifest and then filter the file
# Annotate the gene names to entrez gene IDs using annotation file from HG38
annot = pd.read_csv(os.path.join(supp_dir, 'hsapiens_gene_ensembl-GRCh38.p13.csv'))
annot = annot.dropna(subset=['external_gene_name', 'entrezgene_id'])
annot = annot.drop_duplicates(subset='external_gene_name')
name_to_entrez = dict(zip(annot.external_gene_name, annot.entrezgene_id))

# RUN
rcm = SciRCM(meth_file=DNA_methylation_file,
             rna_file=RNA_file,
             proteomics_file=protein_file,
             rna_logfc=logFC_RNA_column,
             rna_padj=padj_RNA_column,
             meth_diff=betadiff_DNA_methylation_column,
             meth_padj=padj_DNA_methylation_column,
             prot_logfc=logFC_protein_column,
             prot_padj=padj_protein_column,
             gene_id=gene_id_column,
             sep=',',
             rna_padj_cutoff=rna_padj_cutoff,
             prot_padj_cutoff=prot_padj_cutoff,
             meth_padj_cutoff=meth_padj_cutoff,
             rna_logfc_cutoff=rna_logfc_cutoff,
             prot_logfc_cutoff=prot_logfc_cutoff,
             meth_diff_cutoff=meth_diff_cutoff,
             output_dir=output_dir,
             non_coding_genes=['None'],
             output_filename=sircle_file_name,
             bg_type=bg_type
             )
rcm.run()
df = rcm.get_df()

# Fix missing ID issues
df['entrezgene_id'] = [name_to_entrez.get(g) for g in df['gene_name'].values]
df.to_csv(os.path.join(output_dir, f'{sircle_file_name}'), index=False)

#  ----------------------------------------------------------------------------------------------
#                   Make figures
#  ----------------------------------------------------------------------------------------------
# Make some plots
reg_label = 'RG2_Changes'
# figure size in inches
rcParams['figure.figsize'] = 3, 2
sns.set(rc={'figure.figsize': (3, 2)}, style='ticks')

colour_map = {'MDS': '#d8419b', 'MDS_TMDE': '#e585c0', 'MDS_ncRNA': '#d880b4',
              'MDE': '#6aaf44', 'MDE_TMDS': '#0e8e6d', 'MDE_ncRNA': '#9edb77',
              'TMDE': '#fe2323', 'TMDS': '#2952ff', 'TPDE': '#e68e25', 'TPDE_TMDS': '#844c0f',
              'TPDS': '#462d76', 'TPDS_TMDE': '#9b29b7'}

sns.set_palette("Greys_r")
rcm_labels = ["TMDE", "TMDS", "TPDE_TMDS", "TPDE", "TPDS_TMDE", "TPDS", "MDS_TMDE", "MDE",  "MDE_TMDS", "MDS"]
colours = [colour_map[c] for c in rcm_labels]

sns.catplot(data=df, x=reg_label + '_filtered', kind="count", order=rcm_labels, palette=sns.color_palette(colours),
            height=4)
plt.xticks(rotation=45, ha='right')
plt.title(f'{label} {reg_label}')
plt.savefig(os.path.join(fig_dir, f'Figure4A_barplot_{label}_TvN.svg'))
plt.show()

reg_label = 'RG3_Translation'
sns.catplot(data=df, x=reg_label + '_filtered', kind="count",
            order=[c for c in set(df[reg_label].values) if c != 'None'], height=4)
plt.xticks(rotation=45, ha='right')
plt.title(f'{label} {reg_label}')
plt.savefig(os.path.join(fig_dir, f'Figure4A_barplot_{label}_{reg_label}_TvN.svg'))
plt.show()

