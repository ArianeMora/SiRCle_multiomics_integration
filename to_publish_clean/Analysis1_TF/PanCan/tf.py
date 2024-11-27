# Imports
from sciutil import *
import os
from scimotf import *

# Setup file locations and label of the cancer
u = SciUtil()
cancer = 'PanCan'
input_dir = 'Input_RCM'
output_dir = 'Output_Data'
supp_dir = 'Required_Refs'
fig_dir = 'Output_Figures/'
regLabel = 'RG2_Changes_filtered'
clusters = ['TPDE_TMDS', 'TPDS_TMDE', 'MDE_TMDS', 'MDS_TMDE', 'TMDE', 'TMDS']

tf_file = os.path.join(supp_dir, 'DoRoTHea.csv')
sircleFileName = os.path.join(input_dir, f'sircle_PorMandR_{cancer}.csv')

mo = SciMotf_Doro(doro_file=tf_file, cluster_file=sircleFileName, cluster_id=regLabel,
                  cluster_gene_id='gene_name', # got to match motif
                  tf_in_dataset=False,
                   padj_protein='padj_protein', logfc_protein='logFC_protein', padj_rna='padj_rna',
                  logfc_rna='logFC_rna', output_dir=output_dir)

df = mo.run(doro_level=['A', 'B', 'C', 'D'], rcm_clusters=clusters)
df.to_csv(os.path.join(output_dir, f'scimotif_DORO_ABC_filtered.csv'))

# Plot the tfs
plot_cluster_tf(os.path.join(output_dir, f'scimotif_DORO_ABC_filtered.csv'), gene_ratio_min=1, padj_max=0.1, title='',
                fig_dir=fig_dir, label_font_size=9, figsize=(3, 3), axis_font_size=6,
                rcm_labels=clusters,
                save_fig=True)