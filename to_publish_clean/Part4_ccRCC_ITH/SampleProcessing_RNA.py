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
outlier_threshold = 3.0
# Make a log file
output_dir = f'.'

save_fig = True
plot_fig = True

fig_dir = 'figs/'

logfile = os.path.join(output_dir, f'{disease}_RNA_filterlog.tsv')
logfile = open(logfile, 'w+')
sample_df = pd.read_csv(os.path.join(output_dir, f'{disease}_samples_RNA.csv'))
df = pd.read_csv(os.path.join(output_dir, f'{disease}_RNA.csv'))

# -------------------------------------------
#       Write the original sizes out
# -------------------------------------------

cols = list(sample_df['Sample'].values)
logfile.write(f'Original Samples\t{",".join(cols)}\n')
logfile.write(f'Original Size\t{df.shape}\n')
u.dp(['RNA size: ', df.shape])

# -------------------------------------------
#       Drop genes with no gene name
# -------------------------------------------
df = df.dropna(subset=['gene_name'])
u.dp(['RNA size after dropping missing gene names: ', df.shape])
logfile.write(f'Size after dropping genes with no name\t{df.shape}\n')

# -------------------------------------------
#      Drop genes with 50% missing values
# -------------------------------------------
# Drop RNA samples with > 50% missing values
df = df[df.isnull().sum(axis=1) < len(df.values[0]) / 2]
df = df[(df == 0).astype(int).sum(axis=1) < len(df.values[0]) / 2]

u.dp(['After dropping rows with 50% nulls:', df.shape])
logfile.write(f'Size after dropping genes with 50% missing values\t{df.shape}\n')

df = df.fillna(0)

# -------------------------------------------
#       Compute mean RNA
# -------------------------------------------
mean_rna = np.nanmean(df[cols].values, axis=1)
u.dp(['RNA size: ', df.shape, 'Mean counts:', np.mean(mean_rna)])
logfile.write(f'Mean RNA\t{np.mean(mean_rna)}\n')

df = df[mean_rna > 10]
u.dp(['RNA size after 10 count filter: ', df.shape])
logfile.write(f'RNA size after 10 count mean filter:\t{df.shape}\n')

# -------------------------------------------
#       Normalise RNA for correlation
# -------------------------------------------
for c in cols:
    df[c] = np.log2(df[c].values + 1)

corr = df[cols].corr()

# Print out the minimum correlation:
mean_cor = np.nanmedian(corr, axis=1)
corr['mean_corr'] = mean_cor
corr.sort_values(by=['mean_corr'])

m_corr = np.nanmedian(mean_cor)
# Plot out the mean correlation values so we can choose a good filter.
h = Histogram(corr, x='mean_corr', title=f'Mean corr. {m_corr}')
if save_fig:
    plt.savefig(os.path.join(fig_dir, f'{disease}_RNA_corrHist.svg'))
if plot_fig:
    plt.show()

u.dp(['Mean corr: ', m_corr, 'std corr', np.std(mean_cor)])
logfile.write(f'Mean Pearsons correlation\t{m_corr}\n')
logfile.write(f'Std Pearsons correlation\t{np.std(mean_cor)}\n')


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
    plt.savefig(os.path.join(fig_dir, f'{disease}_Tumor_RNA_corrHist.svg'))
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

u.dp(['RNA size after correlation filter: ', np.nanmean(corr_sorted.mean_corr), cutoff, df.shape])

cols_to_omit = [c for c in corr_sorted.index]

logfile.write(f'RNA columns to omit\t{",".join(cols_to_omit)}\n')
u.dp(['RNA columns to omit: '])
print('\n'.join(cols_to_omit))

cols_to_keep = [c for c in df.columns if c not in cols_to_omit]
df = df[cols_to_keep]

u.dp(['RNA shape after dropping tumour columns:', df.shape])
logfile.write(f'RNA size after correlation filter\t{df.shape}\n')

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
    plt.savefig(os.path.join(fig_dir, f'{disease}_Normal_RNA_corrHist.svg'))
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

u.dp(['RNA size after correlation filter: ', np.nanmedian(corr_sorted.mean_corr), cutoff, df.shape])

cols_to_omit = [c for c in corr_sorted.index]

logfile.write(f'RNA columns to omit\t{",".join(cols_to_omit)}\n')
u.dp(['RNA columns to omit: '])
print('\n'.join(cols_to_omit))

cols_to_keep = [c for c in df.columns if c not in cols_to_omit]
df = df[cols_to_keep]

u.dp(['RNA shape after dropping tumour columns:', df.shape])
logfile.write(f'RNA size after correlation filter\t{df.shape}\n')

# -------------------------------------------
#    Filter sample df to only include samples passing QC
# -------------------------------------------
sample_df = sample_df[sample_df['Sample'].isin(cols_to_keep)]

# -------------------------------------------
#    Visualise to check using PCA
# -------------------------------------------
cols = list(sample_df['Sample'].values)
vals = df[cols].values.T

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

sc = Scatterplot(vis_df, x='PC_1', y='PC_2', title=f'Sample Type', xlabel='PC 1', ylabel='PC 2',
                 add_legend=True, colour=vis_df['Colour'].values,
                 config={'s': 40, 'opacity': 0.6, 'figsize': (3, 3)})
sc.plot()
if save_fig:
    plt.savefig(os.path.join(fig_dir, f'{disease}_RNA_scatterPCASampleType.svg'))
if plot_fig:
    plt.show()

sc = Scatterplot(vis_df, x='PC_1', y='PC_2', title=f'Stage', xlabel='PC 1', ylabel='PC 2',
                 add_legend=True, colour=vis_df['StageColour'].values,
                 config={'s': 40, 'opacity': 0.6, 'figsize': (3, 3)})
sc.plot()
if save_fig:
    plt.savefig(os.path.join(fig_dir, f'{disease}_RNA_scatterPCAStage.svg'))
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

sample_df_dedup.to_csv(os.path.join(output_dir, f'{disease}_filtered_samples_RNA.csv'), index=False)

# -------------------------------------------
#      Read back in then save (i.e. so we don't get logged values)
# -------------------------------------------
# -------------------------------------------
#      Drop genes with 50% missing values
# -------------------------------------------
df = pd.read_csv(os.path.join(output_dir, f'{disease}_RNA.csv'))
df = df[df.isnull().sum(axis=1) < len(df.values[0]) / 2]
df = df[(df == 0).astype(int).sum(axis=1) < len(df.values[0]) / 2]
df = df.fillna(0)
# -------------------------------------------
#       Compute mean RNA
# -------------------------------------------
mean_rna = np.nanmean(df[cols].values, axis=1)
df = df[mean_rna > 10]
# -------------------------------------------
#      Save DF
# -------------------------------------------
df[['gene_id', 'gene_name'] + samples].to_csv(os.path.join(output_dir, f'{disease}_filtered_RNA.csv'), index=False)

logfile.close()
