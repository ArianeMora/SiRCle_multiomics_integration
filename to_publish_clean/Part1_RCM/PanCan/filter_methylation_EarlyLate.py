# Imports
import pandas as pd
import numpy as np
from sciutil import *
import os

# Setup file locations and label of the cancer
u = SciUtil()

#  ----------------------------------------------------------------------------------------------
#                  File and naming setup
#  ----------------------------------------------------------------------------------------------
cancer = 'PanCan'
input_dir = 'Input_Methylation'
output_dir = 'Output_Data'
fig_dir = 'Output_Figures'
supp_dir = 'Required_Refs'

#  ----------------------------------------------------------------------------------------------
#                   Filter DNA methylation data if you haven't done this already
#  ----------------------------------------------------------------------------------------------

# Read in data
epic_manifest = pd.read_csv(os.path.join(supp_dir, f'infinium-methylationepic-v-1-0-b5-manifest-file.csv'), comment='#')
epic_manifest.set_index('IlmnID', inplace=True)

# Join the DCpG file with the epic manifest and then filter the file
# Annotate the gene names to entrez gene IDs using annotation file from HG38
annot = pd.read_csv(os.path.join(supp_dir, f'hsapiens_gene_ensembl-GRCh38.p13.csv'))
annot = annot.dropna(subset=['external_gene_name', 'entrezgene_id'])
annot = annot.drop_duplicates(subset='external_gene_name')
name_to_entrez = dict(zip(annot.external_gene_name, annot.entrezgene_id))

# Map the CpGs to their genes and then filter
dcpg = os.path.join(input_dir, f'{cancer}_filtered_DMC_CpG_EarlyLate.csv')
dcpg = pd.read_csv(dcpg, index_col=0)

cpg_data_df = pd.read_csv(os.path.join(input_dir, f'{cancer}_CpG.csv'))
cpg_data_df.set_index('id', inplace=True)
dcpg = dcpg.join(cpg_data_df, how='left')

# Read in the sample DF and then compute the Beta difference
cpg_samples = pd.read_csv(os.path.join(input_dir, f'{cancer}_samples_CpG.csv'))
# First make sure it is only Tumour samples
cpg_samples = cpg_samples[cpg_samples['CondID'] == 1]
tumor_samples = list(cpg_samples[cpg_samples['Stage'] == 'Late']['Sample'].values)
normal_samples = list(cpg_samples[cpg_samples['Stage'] == 'Early']['Sample'].values)
dcpg['beta_diff'] = np.mean(cpg_data_df[tumor_samples].values, axis=1) - np.mean(cpg_data_df[normal_samples].values, axis=1)

# Filter DNA methylation to genes
def filter_methylation_data_by_genes(cpg_df, gene_id, p_val, logfc):
    cpg_df_grped = cpg_df.groupby(gene_id)
    rows = []
    num_cpgs = []
    for cpg in cpg_df_grped:
        cpg = cpg[1]
        cpg = cpg[cpg[p_val] < 0.05]
        num_cpgs.append(len(cpg))
        if len(cpg) > 0:
            if len(cpg) < 3:
                add_row = True
            else:
                pos_cpg = cpg[cpg[logfc] > 0]
                neg_cpg = cpg[cpg[logfc] < 0]
                num_pos = len(pos_cpg)
                num_neg = len(neg_cpg)
                add_row = False
                if num_pos and num_pos/len(cpg) > 0.6:
                    cpg = pos_cpg
                    add_row = True
                elif num_neg and num_neg/len(cpg) > 0.6:
                    cpg = neg_cpg
                    add_row = True
            if add_row:
                max_cpg_idx = None
                max_t_value = 0  # absolute
                idxs = cpg.index
                for xi, t in enumerate(cpg[logfc].values):
                    if abs(t) > abs(max_t_value):
                        max_t_value = t
                        max_cpg_idx = xi
                rows.append(cpg[cpg.index == idxs[max_cpg_idx]].values[0])
    new_cpg_df = pd.DataFrame(rows, columns=cpg_df.columns)
    u.dp(['Originally had: ', len(cpg_df_grped), 'genes.\n', 'Filtered DF now has: ', len(new_cpg_df), ' genes.'])
    return new_cpg_df


dcpg = dcpg.join(epic_manifest, how='left')
dcpg['Name'] = dcpg.index
# For each of the overlapping genes
rows = []
values = dcpg.values
for i, g in enumerate(dcpg['UCSC_RefGene_Name'].values):
    if isinstance(g, str):
        genes = g.split(';')
        for gene_name in genes:
            vals = list(values[i])
            vals.append(gene_name)
            rows.append(vals)

columns = list(dcpg.columns)
columns.append('gene_name')
meth_df = pd.DataFrame(data=rows, columns=columns)
meth_df.rename(columns={'logFC': 'M_diff'}, inplace=True)

filtered_dcpg = filter_methylation_data_by_genes(meth_df, 'gene_name', 'adj.P.Val', 'beta_diff')

filtered_dcpg.to_csv(os.path.join(input_dir, f'{cancer}_filtered_DCpG_EarlyLate.csv'), index=False)
