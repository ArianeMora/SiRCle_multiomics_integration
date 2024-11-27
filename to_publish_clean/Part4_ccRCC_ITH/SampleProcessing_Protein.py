import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sciviso import *
from sciutil import SciUtil
import os
import seaborn as sns
# Have a look at clustering each of these
from sklearn.decomposition import PCA

OPP = ['#51AE66', '#AE5199', '#4325DA', '#BCDA25']

u = SciUtil()

disease = 'ClearCellRenalCellCarcinoma'

# https://cptac-data-portal.georgetown.edu/cptac/documents/CDAP_Results_Overview_rev_09152014.pdf
input_folder = '../../input_data/ClearCellRenalCellCarcinoma_PDC000127'
disease_space = 'Clear Cell Renal Cell Carcinoma'
sample_df = pd.read_csv(os.path.join(input_folder, f'PDC_study_biospecimen_07172023_114806.tsv'), sep='\t')
df = pd.read_csv(os.path.join(input_folder, f'CPTAC3_Clear_Cell_Renal_Cell_Carcinoma_Proteome.tmt10.tsv'), sep='\t')

# Make a log file
output_dir = f''
outlier_threshold = 3.0
save_fig = True
plot_fig = True

fig_dir = f'figs'

logfile = os.path.join(output_dir, f'{disease}_Protein_filterlog.tsv')
logfile = open(logfile, 'w+')

# Get the summary file and the biospecimen file
case_df = pd.read_csv(os.path.join(output_dir, f'{disease}_case.csv'))
sample_df['SafeCases'] = [c.replace('-', '.') for c in sample_df['Case Submitter ID'].values]
sample_df.set_index('SafeCases', inplace=True)
case_df.set_index('SafeCases', inplace=True)
sample_df = sample_df.join(case_df, how='left', lsuffix='_protein')

# ----------- Ensure all the disease is the same
sample_df = sample_df[sample_df['Disease Type'] == disease_space]

# ----------- Remove std row
gene_drop_columns = ['Mean', 'Median', 'StdDev']
df = df[~df['Gene'].isin(gene_drop_columns)]

plt.rcParams["figure.figsize"] = (20, 3)
cols = [c for c in df.columns if 'Log' in c and 'Unshared' not in c]
X = df[cols].values
df.boxplot(column=cols)
plt.xticks(rotation=90)
plt.show()

# Drop columns that dont have CPT in them
df = df[['Gene'] + [c for c in df.columns if 'CPT' in c and 'Unshared' not in c]]

# Drop outliers from each gene (this could affect the DA)
cols = [c for c in df.columns if 'Log' in c and 'Unshared' not in c]

df.boxplot(column=cols)
plt.xticks(rotation=90)
plt.show()

# Check the distribution
plt.rcParams["figure.figsize"] = (4, 3)
plt.hist(np.nanmedian(df[cols].values, axis=1), bins=20)
plt.show()

# Set NaNs for values that are 2SD or greater from the mean of the gene
median_protein = np.nanmedian(df[cols].values, axis=1)
std_protein = np.nanstd(df[cols].values, axis=1)
print(median_protein, std_protein)
upper_cutoff = median_protein + (outlier_threshold * std_protein)
lower_cutoff = median_protein - (outlier_threshold * std_protein)
for col in cols:
    df[col][df[col] > upper_cutoff] = None
    df[col][df[col] < lower_cutoff] = None

plt.rcParams["figure.figsize"] = (20, 3)

df.boxplot(column=cols)
plt.xticks(rotation=90)
plt.show()

plt.rcParams["figure.figsize"] = (20, 3)

df.boxplot(column=cols)
plt.xticks(rotation=90)
plt.show()

# Drop genes with > 25% missing values (including those with Nulls that we put in )
df = df[df.isnull().sum(axis=1) < len(df.values[0]) / 4]

# -------------------------------------------
#       Rename the columns so that they are in the same format as the other data types
# -------------------------------------------
# C3L.00094_Tumor_RNA_LungAdenocarcinoma_167b35b1.7858.4b23.8b87.17c660b5d658
col_map = {'Gene': 'gene_name'}
sample_types = sample_df['Sample Type'].values
sample_df['SafeCases'] = sample_df.index
submitter_ids = sample_df['SafeCases'].values
tumour_type = []
samples = []
cond_ids = []
for i, aliquot in enumerate(sample_df['Aliquot Submitter ID'].values):
    if sample_types[i] == 'Primary Tumor':
        sample_type = 'Tumor'
    elif sample_types[i] == 'Solid Tissue Normal':
        sample_type = 'Normal'
    else:
        print(sample_types[i])
        sample_type = 'NA'
    tumour_type.append(sample_type)
    col_map[f'{aliquot} Log Ratio'] = f'{submitter_ids[i]}_{sample_type}_Protein_{aliquot}'
    samples.append(f'{submitter_ids[i]}_{sample_type}_Protein_{aliquot}')
    cond_id = 0 if sample_type == 'Normal' else 1
    cond_ids.append(cond_id)
sample_df['Sample'] = samples
sample_df['CondID'] = cond_ids
sample_df['CondLabel'] = tumour_type

sample_df = sample_df[sample_df['CondLabel'] != 'NA']

df = df.rename(columns=col_map)

# ---------------------------------------------
#       Compute correlation for tumour samples
# ---------------------------------------------
cols = [c for c in df.columns if c != 'gene_name' and 'Tumor' in c]

corr = df[cols].corr()

# Print out the minimum correlation:
mean_cor = np.nanmedian(corr, axis=1)
corr['mean_corr'] = mean_cor
corr.sort_values(by=['mean_corr'])

m_corr = np.nanmedian(mean_cor)
# Plot out the mean correlation values so we can choose a good filter.
h = Histogram(corr, x='mean_corr', title=f'Mean corr. {m_corr}')
if save_fig:
    plt.savefig(os.path.join(fig_dir, f'{disease}_Protein_corrHist.svg'))
if plot_fig:
    plt.show()

u.dp(['Mean corr: ', m_corr, 'std corr', np.std(mean_cor)])
logfile.write(f'Mean Pearsons correlation for Tumour samples\t{m_corr}\n')
logfile.write(f'Std Pearsons correlation for Tumour samples\t{np.std(mean_cor)}\n')

# -------------------------------------------
#       Filter patients on correlation
# -------------------------------------------
corr_sorted = corr.sort_values(by=['mean_corr'])
cutoff = np.mean(corr_sorted.mean_corr) - (outlier_threshold * np.std(corr_sorted.mean_corr))
corr_sorted = corr_sorted[corr_sorted['mean_corr'] < cutoff]

u.dp(['Protein size after correlation filter: ', np.nanmedian(corr_sorted.mean_corr), cutoff, df.shape])

cols_to_omit = [c for c in corr_sorted.index]

logfile.write(f'Protein columns to omit\t{",".join(cols_to_omit)}\n')
u.dp(['Protein columns to omit: '])
print('\n'.join(cols_to_omit))

cols_to_keep = [c for c in df.columns if c not in cols_to_omit]
df = df[cols_to_keep]

u.dp(['Protein shape after dropping tumour columns:', df.shape])
logfile.write(f'Protein size after correlation filter\t{df.shape}\n')

# -------------------------------------------
#    Filter sample df to only include samples passing QC
# -------------------------------------------
sample_df = sample_df[sample_df['Sample'].isin(cols_to_keep)]

# ---------------------------------------------
#       Compute correlation for Normal samples
# ---------------------------------------------
cols = [c for c in df.columns if c != 'gene_name' and 'Normal' in c]

corr = df[cols].corr()

# Print out the minimum correlation:
mean_cor = np.nanmedian(corr, axis=1)
corr['mean_corr'] = mean_cor
corr.sort_values(by=['mean_corr'])

m_corr = np.nanmedian(mean_cor)
# Plot out the mean correlation values so we can choose a good filter.
h = Histogram(corr, x='mean_corr', title=f'Mean corr. {m_corr}')
if save_fig:
    plt.savefig(os.path.join(fig_dir, f'{disease}_Protein_corrHist.svg'))
if plot_fig:
    plt.show()

u.dp(['Mean corr: ', m_corr, 'std corr', np.std(mean_cor)])
logfile.write(f'Mean Pearsons correlation for Tumour samples\t{m_corr}\n')
logfile.write(f'Std Pearsons correlation for Tumour samples\t{np.std(mean_cor)}\n')

# -------------------------------------------
#       Filter patients on correlation
# -------------------------------------------
corr_sorted = corr.sort_values(by=['mean_corr'])
cutoff = np.nanmedian(corr_sorted.mean_corr) - (outlier_threshold * np.std(corr_sorted.mean_corr))
corr_sorted = corr_sorted[corr_sorted['mean_corr'] < cutoff]

u.dp(['Protein size after correlation filter: ', np.nanmedian(corr_sorted.mean_corr), cutoff, df.shape])

cols_to_omit = [c for c in corr_sorted.index]

logfile.write(f'Protein columns to omit\t{",".join(cols_to_omit)}\n')
u.dp(['Protein columns to omit: '])
print('\n'.join(cols_to_omit))

cols_to_keep = [c for c in df.columns if c not in cols_to_omit]
df = df[cols_to_keep]

u.dp(['Protein shape after dropping tumour columns:', df.shape])
logfile.write(f'Protein size after correlation filter\t{df.shape}\n')

# -------------------------------------------
#    Filter sample df to only include samples passing QC
# -------------------------------------------
sample_df = sample_df[sample_df['Sample'].isin(cols_to_keep)]

# -------------------------------------------
#    Visualise to check using PCA
# -------------------------------------------
cols = list(sample_df['Sample'].values)

vals = df[cols].values.T

# For vis replace Nan with mean - we'll impute later
vals = np.nan_to_num(vals, np.nanmedian(np.nanmedian(vals)))

pca = PCA(n_components=2)
pca_values = pca.fit_transform(vals)
var_ratio = pca.fit(vals).explained_variance_ratio_
plt.rcParams['figure.figsize'] = [4, 4]
vis_df = pd.DataFrame()
vis_df['PC_1'] = pca_values[:, 0]
vis_df['PC_2'] = pca_values[:, 1]
vis_df['Stage'] = sample_df['TumorStage'].values
vis_df['Disease'] = sample_df['Disease'].values
vis_df['CondID'] = sample_df['CondID'].values
vis_df['Colour'] = ['#C24B3D' if c == 1 else '#3DB4C2' for c in vis_df['CondID'].values]
stage_c_map = {'Stage I': '#3139ba', 'Stage II': '#7169E4', 'Stage III': '#2A9D8F',
               'Stage IV': '#264653'}
vis_df['StageColour'] = [stage_c_map.get(c) if stage_c_map.get(c) else '#808080' for c in vis_df['Stage'].values]
vis_df['Primary Site'] = sample_df['Primary Site'].values

sc = Scatterplot(vis_df, x='PC_1', y='PC_2', title=f'Sample Type', xlabel='PC 1', ylabel='PC 2',
                 add_legend=True, colour=vis_df['Colour'].values,
                 config={'s': 40, 'opacity': 0.6, 'figsize': (3, 3)})
sc.plot()
if save_fig:
    plt.savefig(os.path.join(fig_dir, f'{disease}_Protein_scatterPCASampleType.svg'))
if plot_fig:
    plt.show()

sc = Scatterplot(vis_df, x='PC_1', y='PC_2', title=f'Stage', xlabel='PC 1', ylabel='PC 2',
                 add_legend=True, colour=vis_df['StageColour'].values,
                 config={'s': 40, 'opacity': 0.6, 'figsize': (3, 3)})
sc.plot()
if save_fig:
    plt.savefig(os.path.join(fig_dir, f'{disease}_Protein_scatterPCAStage.svg'))
if plot_fig:
    plt.show()

sc = Scatterplot(vis_df, x='PC_1', y='PC_2', title=f'Primary Site', xlabel='PC 1', ylabel='PC 2',
                 add_legend=True, color_col='Primary Site',
                 config={'s': 40, 'opacity': 0.6, 'figsize': (3, 3)})
sc.plot()
if save_fig:
    plt.savefig(os.path.join(fig_dir, f'{disease}_Protein_scatterPCASite.svg'))
if plot_fig:
    plt.show()

# -------------------------------------------
# Drop duplicates based on the case ID and the condition type of samples
# -------------------------------------------
u.dp(['Before dropping duplicate samples for a patient', sample_df.shape])
sample_df_dedup = sample_df.drop_duplicates(subset=['SafeCases', 'CondID'])
u.dp(['After dropping duplicate samples for a patient', sample_df_dedup.shape])
logfile.write(f'After dropping duplicates from samples\t{sample_df_dedup.shape}\n')
logfile.write(
    f'Duplicate samples that were dropped\t{",".join([c for c in sample_df.Sample if c not in list(sample_df_dedup.Sample)])}\n')
u.dp(['Cases included in dataset', len(list(set(sample_df_dedup.SafeCases.values)))])

logfile.write(f'Cases included in dataset\t{",".join(list(set(sample_df_dedup.SafeCases.values)))}\n')

samples = list(sample_df_dedup['Sample'].values)

# Drop the outliers and shift the 0 to be real numbered so that we can have a proper logFC
df_min = min(df[samples].min())

for col in samples:
    df[col] = df[col].values + abs(df_min)

plt.rcParams["figure.figsize"] = (20, 3)

df.boxplot(column=cols)
plt.xticks(rotation=90)
plt.show()

sample_df_dedup.to_csv(os.path.join(output_dir, f'{disease}_filtered_samples_Protein.csv'), index=False)
df[['gene_name'] + samples].to_csv(os.path.join(output_dir, f'{disease}_filtered_Protein.csv'), index=False)

logfile.close()
