# Imports
import pandas as pd
from tqdm import tqdm
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from sciviso import Scatterplot
# Imports
from sciutil import SciUtil
import os


# Setup file locations and label of the cancer
u = SciUtil()
from collections import defaultdict

input_dir = f'Input_RCM'
output_dir = f'Output_Data/'
supp_dir = f'Required_Refs/'
fig_dir = f'Output_Figures/'

train = True
# Set regulatory clusters you're interested in investigating
reg_label = 'RG2_Changes_filtered'  # Use these clusters to build a VAE for each one could be something else!
cancer = 'PanCan'
missing_method = 'mean'

column_id = 'FullLabel'
condition_column = 'CondID'
patient_id_column = 'SafeCases'
clinical_label = 'Stage'

feature_columns = ['RNA-LogFC',
                  'Protein-LogFC',
                  'CpG-LogFC',
                  'RNA-Tumor',
                  'RNA-Normal',
                  'Protein-Tumor',
                  'Protein-Normal']
meta_rcm_df = pd.read_csv(os.path.join(input_dir, f'RCM_{cancer}.csv'), index_col=0)
files = os.listdir(output_dir)
stats_files = [f for f in files if 'stats' in f]

rcm_labels = ["MDS", "MDS_TMDE", "MDE", "MDE_TMDS", "TMDE", "TMDS", "TPDE", "TPDE_TMDS", "TPDS", "TPDS_TMDE"]
comparisons = ['Stage IV-Stage I'] #'Late-Early']

vae_rcm_df = meta_rcm_df.copy()  # Make a copy that we'll add to
encoded_df = pd.read_csv(os.path.join(output_dir, f'encoded_df_{cancer}.csv'))
raw_input_df = pd.read_csv(os.path.join(output_dir, f'raw_input_df_{cancer}.csv'))
patient_info = pd.read_csv(os.path.join(input_dir, f'patient_info_{cancer}.csv'))

for comp in comparisons:
    cond0 = comp.split('-')[1]
    cond1 = comp.split('-')[0]
    for f in tqdm(stats_files):
        if cancer in f and cond0 in f and cond1 in f and 'stats' in f:
            df = pd.read_csv(os.path.join(output_dir, f'{f}'), index_col=0)
            df.set_index('id', inplace=True)
            df_merged = df.join(vae_rcm_df, how='left', lsuffix='_')  # Join to the
            df_merged = df_merged.drop_duplicates()
            # Map gene names and move folder
            df_merged = df_merged[[
                reg_label,
                f'Integrated padj ({cond1}-{cond0})',
                f'Integrated pval ({cond1}-{cond0})',
                f'Integrated diff ({cond1}-{cond0})',
                f'mannwhitneyu stat ({cond1}-{cond0})',
                f'Integrated mean ({cond0})',
                f'Integrated mean ({cond1})',
                f'Protein-LogFC mean ({cond1}-{cond0})',
                f'RNA-LogFC mean ({cond1}-{cond0})',
                f'CpG-LogFC mean ({cond1}-{cond0})',
                'entrezgene_id',
                'logFC_rna',
                'padj_rna',
                'beta_diff',
                'adj.P.Val',
                'logFC_protein',
                'padj_protein']]  # Reorder columns
            df_merged.to_csv(os.path.join(output_dir, f'{missing_method}_Integrated_comparison_{cond1}-{cond0}.csv'))
            print(comp)
            for c in df_merged.columns:
                if 'Integrated diff' in c or 'Integrated padj' in c or 'mean' in c or 'Integrated pval' in c:
                    vae_rcm_df.loc[df_merged.index.values, c] = df_merged[c].values

diff_df = vae_rcm_df[[c for c in vae_rcm_df if ' diff' in c or 'LogFC' in c]]
corr = diff_df.corr()
save_fig = True
sns.clustermap(corr, vmin=-1, vmax=1,
               xticklabels=corr.columns.values,
               yticklabels=corr.columns.values, cmap='RdBu_r', row_cluster=True, col_cluster=True)
if save_fig:
    plt.savefig(os.path.join(fig_dir, f'{missing_method}_VAE_data_corr.svg'), transparent=True)

p_cutoff = 0.05

for comparison in comparisons:
    fig, axs = plt.subplots(10, 3, figsize=(6, 21))
    vae_change = f'Integrated diff ({comparison})'
    vae_p = f'Integrated padj ({comparison})'

    protein_change = f'Protein-LogFC mean ({comparison})'
    rna_change = f'RNA-LogFC mean ({comparison})'
    cpg_change = f'CpG-LogFC mean ({comparison})'

    for fi, r in enumerate(rcm_labels):
        try:
            tmp_df = vae_rcm_df[vae_rcm_df[reg_label] == r]

            opts = {'title_font_weight': 'bold', "title_font_size": 9, "label_font_size": 7, "figsize": (6, 9),
                    "s": 20}
            corr_c = round(np.corrcoef(tmp_df[vae_change], tmp_df[cpg_change])[0][1], 2)
            scatter = Scatterplot(tmp_df, vae_change, cpg_change,
                                  colour='orange', config=opts)
            grp_labels = ['Integrated']
            vae_sig = []
            vae_int_change = abs(tmp_df[vae_change].values)
            vae_cutoff = 0.5
            for i, vae_padj in enumerate(tmp_df[vae_p].values):
                if vae_padj < p_cutoff and vae_int_change[i] > vae_cutoff:
                    vae_sig.append(i)
            grp_idxs = [vae_sig]
            scatter.plot_groups_2D(axs[fi][0], grp_labels, grp_idxs, ['orange'])

            corr_r = round(np.corrcoef(tmp_df[vae_change], tmp_df[rna_change])[0][1], 2)
            scatter = Scatterplot(tmp_df, vae_change, rna_change, colour='rebeccapurple',
                                  xlabel='Integrated Rank (Stage 4 - Stage 1)',
                                  config=opts)
            scatter.plot_groups_2D(axs[fi][1], grp_labels, grp_idxs, ['orange', 'rebeccapurple', 'darkred'])

            corr_p = round(np.corrcoef(tmp_df[vae_change], tmp_df[protein_change])[0][1], 2)
            scatter = Scatterplot(tmp_df, vae_change, protein_change, colour='grey',
                                  config=opts)
            scatter.plot_groups_2D(axs[fi][2], grp_labels, grp_idxs, ['orange', 'rebeccapurple', 'darkred'])
            plt.tight_layout()

            axs[fi][0].set_title(f'CpG ({corr_c})', fontweight="bold")
            axs[fi][1].set_title(f'{r} \n \n RNA ({corr_r})', fontweight="bold")
            axs[fi][2].set_title(f'Protein ({corr_p})', fontweight="bold")

            axs[fi][0].set(xlabel='Integrated diff', ylabel='CpG diff')
            axs[fi][1].set(xlabel='Integrated diff', ylabel='RNA LogFC')
            axs[fi][2].set(xlabel='Integrated diff', ylabel='Protein LogFC')
        except:
            print(fi)
    plt.savefig(os.path.join(fig_dir, f'VAE_data_corr_PLOTS_change_{cancer}.svg'), transparent=True)
    plt.show()

for comparison in comparisons:
    fig, axs = plt.subplots(10, 2, figsize=(4, 24))
    cond0 = comparison.split('-')[1]
    cond1 = comparison.split('-')[0]
    for fi, r in enumerate(rcm_labels):
        tmp_df = vae_rcm_df[vae_rcm_df[reg_label] == r]

        axs[fi][0].hist(tmp_df[f'Integrated diff ({comparison})'], color='orange')

        axs[fi][1].hist(tmp_df[f'Integrated mean ({cond0})'], color='blue', alpha=0.2)
        axs[fi][1].hist(tmp_df[f'Integrated mean ({cond1})'], color='red', alpha=0.2)
        axs[fi][1].set_title(r)
    if save_fig:
        plt.savefig(os.path.join(fig_dir, f'{missing_method}_VAE_hists_{comparison}_change_{cancer}.svg'), transparent=True)
    plt.show()

values = defaultdict(lambda: defaultdict(list))
for r in rcm_labels:
    encoding_df = encoded_df[encoded_df['RG2_Changes_filtered'] == r]
    input_df = raw_input_df[raw_input_df['RG2_Changes_filtered'] == r]  # Get the input DF
    # Make a df for each of the data types we want box plots for
    for column in input_df.columns:
        if 'C' in column:
            try:
                c = column.split('_')
                case_info = patient_info[patient_info['SafeCases'] == c[0]]
                # Stage	BAP1&PBRM1 gender	TumorStage	AgeGrouped
                values[c[1]]['Stage'] += [case_info['Stage'].values[0]] * len(input_df)
                values[c[1]]['Tumor Stage'] += [case_info['TumorStage'].values[0]] * len(input_df)

                values[c[1]]['values'] += [v for v in input_df[column].values]
                values[c[1]]['GeneId'] += [v for v in input_df.id.values]
                values[c[1]]['Cluster'] += [v for v in input_df['RG2_Changes_filtered'].values]
            except:
                print("ERROR", column)
    # Now we also want an encoding values
    for column in encoding_df.columns:
        if 'C' in column:
            try:
                case_info = patient_info[patient_info['SafeCases'] == column]

                values['Integrated']['Stage'] += [case_info['Stage'].values[0]] * len(encoding_df)
                values['Integrated']['Tumor Stage'] += [case_info['TumorStage'].values[0]] * len(encoding_df)

                values['Integrated']['values'] += [c for c in encoding_df[column].values]
                values['Integrated']['GeneId'] += [c for c in encoding_df['id'].values]
                values['Integrated']['Cluster'] += [v for v in encoding_df['RG2_Changes_filtered'].values]
            except:
                print("ERROR", column)
# Just formatting the data for a boxplot
dfs = {}
for v in values:
    print(v)
    box_df = pd.DataFrame()
    for c in values[v]:
        box_df[c] = values[v][c]
    dfs[v] = box_df

for d in dfs:
    dfs[d].to_csv(os.path.join(output_dir, f'Boxplot_data_{d}.csv'), index=False)

