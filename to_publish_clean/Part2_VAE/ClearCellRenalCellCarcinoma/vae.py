# Imports
from scircm import RCMStats
from sciutil import SciUtil
import os


# Setup file locations and label of the cancer
u = SciUtil()

input_dir = f'Input_RCM'
output_dir = f'Output_Data/'
supp_dir = f'Required_Refs/'
fig_dir = f'Output_Figures/'

train = False
# Set regulatory clusters you're interested in investigating
reg_label = 'RG2_Changes_filtered'  # Use these clusters to build a VAE for each one could be something else!
cancer = 'ClearCellRenalCellCarcinoma'
missing_method = 'mean'

column_id = 'FullLabel'
condition_column = 'CondID'
patient_id_column = 'SafeCases'
clinical_label = 'TumorStage'

feature_columns = ['RNA-LogFC',
                  'Protein-LogFC',
                  'CpG-LogFC',
                   # 'CpG-Tumor',
                   # 'CpG-Normal',
                  'RNA-Tumor',
                  'RNA-Normal',
                  'Protein-Tumor',
                  'Protein-Normal']
#  ----------------------------------------------------------------------------------------------
#                         VAE configuration
#  ----------------------------------------------------------------------------------------------

epochs = 100
batch_size = 16
num_nodes = 5
mmd_weight = 0.5
loss = {'loss_type': 'mse', 'distance_metric': 'mmd', 'mmd_weight': mmd_weight}

config = {"loss": loss,
          "encoding": {"layers": [{"num_nodes": num_nodes, "activation_fn": "selu"}]},
          "decoding": {"layers": [{"num_nodes": num_nodes, "activation_fn": "selu"}]},
          "latent": {"num_nodes": 1},
          "optimiser": {"params": {'learning_rate': 0.01}, "name": "adam"},
          "epochs": epochs,
          "batch_size": batch_size,
          "scale_data": False
          }

#  ----------------------------------------------------------------------------------------------
#                         Start the VAE/RCM stats
#  ----------------------------------------------------------------------------------------------

sv = RCMStats(rcm_file=os.path.join(input_dir, f'RCM_{cancer}.csv'),
              patient_sample_file=os.path.join(input_dir, f'patient_info_{cancer}.csv'),
              meth_file=os.path.join(input_dir, f'CPTAC_cpg_{cancer}.csv'),
              meth_sample_file=os.path.join(input_dir, f'samples_CpG_{cancer}_VAE.csv'),
              rna_file=os.path.join(input_dir, f'CPTAC_rna_{cancer}.csv'),
              rna_sample_file=os.path.join(input_dir, f'samples_RNA_{cancer}_VAE.csv'),
              protein_file=os.path.join(input_dir, f'CPTAC_protein_{cancer}.csv'),
              protein_sample_file=os.path.join(input_dir, f'samples_protein_{cancer}_VAE.csv'),
              output_folder=output_dir,
              column_id=column_id,
              condition_column=condition_column,
              patient_id_column=patient_id_column,
              clinical_label=clinical_label,
              regulatory_label=reg_label,
              run_name=cancer,
              normalise='rows',
              verbose=True,
              iid=False,
              missing_method=missing_method)

# Get the patient info that has been compiled from the provided sample files
patient_info = sv.patient_clinical_df

#  ----------------------------------------------------------------------------------------------
#                        Pick the cases you want to use for traininh
#  ----------------------------------------------------------------------------------------------

# Select the cases with 5 samples
matching_patient_info = patient_info[patient_info['Sample Counts'] == 6]
matching_cases = matching_patient_info[patient_id_column].values
print("total number of patients: ", len(patient_info), " vs number with matching data: ",
      len(matching_patient_info))

sv.feature_columns = feature_columns

if train:
    sv.train_vae(cases=list(set(matching_cases)), config=config)
    sv.save()  # Save the information we have generated.
else:
    sv.load_saved_vaes()
    sv.load_saved_encodings(os.path.join(output_dir, f'encoded_df_{cancer}.csv'))
    sv.load_saved_inputs(os.path.join(output_dir, f'vae_input_df_{cancer}.csv'))
    sv.load_saved_raws(os.path.join(output_dir, f'raw_input_df_{cancer}.csv'))

#  ----------------------------------------------------------------------------------------------
#                         Make statistics between comparisons
#  ----------------------------------------------------------------------------------------------

u.warn_p(['S4 vs S1'])
sv.run_vae_stats(cond_label='TumorStage', cond0='Stage I', cond1='Stage IV')

u.warn_p(['Late vs Early'])
sv.run_vae_stats(cond_label='Stage', cond0='Early', cond1='Late')