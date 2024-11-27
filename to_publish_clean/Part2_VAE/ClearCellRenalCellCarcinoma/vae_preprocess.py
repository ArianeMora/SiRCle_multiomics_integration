# Imports
import pandas as pd
import numpy as np
from sciutil import SciUtil
import os

# Setup file locations and label of the cancer
u = SciUtil()

input_dir = 'Input_RCM'
output_dir = 'Output_Data'
supp_dir = 'Required_Refs'
fig_dir = 'Output_Figures'
cancer = 'ClearCellRenalCellCarcinoma'
gene_id = 'gene_name'

# Set regulatory clusters you're interested in investigating
rcm_labels = ["MDS", "MDS_TMDE", "MDE", "MDE_TMDS", "TMDE", "TMDS", "TPDE", "TPDE_TMDS",  "TPDS", "TPDS_TMDE"]
reg_label = 'RG2_Changes_filtered'  # Use these clusters to build a VAE for each one could be something else!

# Load in the files that we used for DE analysis and remove any redundant columns
# the VAE expects just the ID as the index, and then the values #ClearCellRenalCellCarcinoma_filtered_samples_RNA
rna_sample_file = pd.read_csv(os.path.join(input_dir, f'{cancer}_filtered_samples_RNA.csv'))
prot_sample_file = pd.read_csv(os.path.join(input_dir, f'{cancer}_filtered_samples_Protein.csv'))
# Now we want to merge the clinical info with the cases from the sample df
meth_sample_file = pd.read_csv(os.path.join(input_dir, f'{cancer}_filtered_samples_CpG.csv'))

# Read in the RCM and filter to only include the regulatory labels of interest
r_df = pd.read_csv(os.path.join(input_dir, f'sircle_PorMandR_{cancer}.csv'))
rcm_df = r_df[r_df[reg_label].isin(rcm_labels)]

# The sample files need to be formatted consistently
meth_sample_file['FullLabel'] = [c for c in meth_sample_file['Sample'].values]
prot_sample_file['FullLabel'] = [c for c in prot_sample_file['Sample'].values]
rna_sample_file['FullLabel'] = [c for c in rna_sample_file['Sample'].values]

meth_sample_file.to_csv(os.path.join(input_dir, f'samples_CpG_{cancer}_VAE.csv'), index=False)
rna_sample_file.to_csv(os.path.join(input_dir, f'samples_RNA_{cancer}_VAE.csv'), index=False)
prot_sample_file.to_csv(os.path.join(input_dir, f'samples_protein_{cancer}_VAE.csv'), index=False)

# Create a patient info sample file
patient_info = pd.DataFrame()
cases = list(set(list(meth_sample_file.SafeCases.values) + list(prot_sample_file.SafeCases.values) + list(
    rna_sample_file.SafeCases.values)))
cols = ['_'.join(c.split('_')[:3]) for c in r_df.columns if len(c.split('_')) > 3]
print(cols)
normal_protein, normal_rna, normal_cpg, tumour_protein, tumour_rna, tumour_cpg, matching, counts = [], [], [], \
                                                                                                   [], [], [], [], []
for c in cases:
    count = 0
    if f'{c}_Normal_RNA' in cols:
        normal_rna.append(f'{c}_Normal_RNA')
        count += 1
    else:
        normal_rna.append(None)
    if f'{c}_Tumor_RNA' in cols:
        tumour_rna.append(f'{c}_Tumor_RNA')
        count += 1
    else:
        tumour_rna.append(None)

    if f'{c}_Normal_Protein' in cols:
        normal_protein.append(f'{c}_Normal_Protein')
        count += 1
    else:
        normal_protein.append(None)

    if f'{c}_Tumor_Protein' in cols:
        tumour_protein.append(f'{c}_Tumor_Protein')
        count += 1
    else:
        tumour_protein.append(None)

    if f'{c}_Normal_CpG' in cols:
        normal_cpg.append(f'{c}_Normal_CpG')
        count += 1
    else:
        normal_cpg.append(None)

    if f'{c}_Tumor_CpG' in cols:
        tumour_cpg.append(f'{c}_Tumor_CpG')
        count += 1
    else:
        tumour_cpg.append(None)

    if count == 6:
        matching.append(True)
    else:
        matching.append(False)
    counts.append(count)

patient_info['SafeCases'] = cases
patient_info['Normal_Protein'] = normal_protein
patient_info['Normal_RNA'] = normal_rna
patient_info['Normal_CpG'] = normal_cpg
patient_info['Tumour_Protein'] = tumour_protein
patient_info['Tumour_RNA'] = tumour_rna
patient_info['Tumour_CpG'] = tumour_cpg
patient_info['Matching'] = matching
patient_info['Sample Counts'] = counts

protein_info = prot_sample_file.copy()
# Only keep those that have at least protein tumour data
patient_info = patient_info.dropna(subset=['Tumour_Protein', 'Normal_Protein'])
protein_info.set_index('SafeCases', inplace=True)
patient_info.set_index('SafeCases', inplace=True)
patient_info = patient_info.join(protein_info, how='left')

# --------- Add in early and late stages for the last reviewer
stages = []
for stage in patient_info['TumorStage'].values:
    if stage == 'Stage I' or stage == 'Stage II':
        stages.append('Early')
    elif stage == 'Stage III' or stage == 'Stage IV':
        stages.append('Late')
    else:
        stages.append(None)
patient_info['Stage'] = stages
# Drop patients with no stage!
patient_info = patient_info[patient_info['Stage'] != None]

patient_info.to_csv(os.path.join(input_dir, f'patient_info_{cancer}.csv'))

# Add ID column in for all datasets
rcm_df['id'] = rcm_df[gene_id].values

meta_cols = ['id', gene_id, 'entrezgene_id', "logFC_rna", "padj_rna",
             "beta_diff", "adj.P.Val", "logFC_protein", "padj_protein",
             "RG2_Changes_filtered"]

# We need to standardise the protein values --> this may not be the case for your data but for ours they were not
# sufficiently close to normal
# for c in rcm_df.columns:
#     if 'Protein' in c and 'C3' in c:
#         mask = np.sign(rcm_df[c].values)
#         rcm_df[c] = mask * np.log2(abs(rcm_df[c].values) + 1)

# Save to input dir
save_input_data = True
rcm_df = rcm_df.drop_duplicates(subset=['gene_name'])
rcm_df = rcm_df[rcm_df[reg_label] != 'None']
rcm_df = rcm_df.fillna(0)
meta_rcm_df = rcm_df[meta_cols].copy()
meta_rcm_df.set_index('id', inplace=True)

# Save to the input data dir folder.
if save_input_data:
    rcm_df[['id'] + [c for c in rcm_df.columns if c in prot_sample_file['FullLabel'].values]].to_csv(
        os.path.join(input_dir, f'CPTAC_protein_{cancer}.csv'), index=False)
    rcm_df[['id'] + [c for c in rcm_df.columns if c in rna_sample_file['FullLabel'].values]].to_csv(
        os.path.join(input_dir, f'CPTAC_rna_{cancer}.csv'), index=False)
    rcm_df[['id'] + [c for c in rcm_df.columns if c in meth_sample_file['FullLabel'].values]].to_csv(
        os.path.join(input_dir, f'CPTAC_cpg_{cancer}.csv'), index=False)
    meta_rcm_df.to_csv(os.path.join(input_dir, f'RCM_{cancer}.csv'))
