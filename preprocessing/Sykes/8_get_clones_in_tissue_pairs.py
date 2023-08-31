"""

"""

import pandas as pd
import numpy as np
from pathlib import Path, PureWindowsPath
import pickle
from numpy import *
import json

subjects_dict = {   
                2: 'Pt19_R',
                3: 'Pt21_R',
                4: 'Pt23_R',
                8: 'Pt14_R',
                9: 'Pt20_R',
                10: 'Pt17_R',
                17: 'Pt25_R',
                }
        
# =============================================================================
# This function returns two dicts:
# samples_per_tissue: dict[tissue] : [list of sample_ids in that tissue]    
# tissues_per_sample: dict[sample_id] : [tissue in that sample_id] ###this will be a list with a single element (the tissue)
# =============================================================================
def get_tissue_sample_ids(patient_ID, patient_name):
  
    tissues_per_sample = {}
    samples_per_tissue = {}
    
    if patient_name in subjects_dict.values():
        samples_table = pd.read_pickle(r"/path_to_sykes_children_data_folder/samples.pkl")
        metadata_table = pd.read_pickle(r"/path_to_sykes_children_data_folder/sample_metadata.pkl")
        
        subject_samples = samples_table['id'].loc[samples_table['subject_id'] == patient_ID].unique() 
        
        metadata_table = metadata_table.loc[metadata_table['sample_id'].isin(subject_samples)]
        metadata_table['value'] = metadata_table['value'].replace("Colon ", "Colon")
        metadata_table['value'] = metadata_table['value'].replace("MLN_allograft", "AxLN")
        metadata_table['value'] = metadata_table['value'].replace("POD37 44", "POD37")
        ###Change Ileum to Ileum_allograft for Pt21_R POD1145
        POD1145_sample_ids = metadata_table['sample_id'].loc[metadata_table['value'] == 'POD1145']
        metadata_table['value'].loc[metadata_table['sample_id'].isin(POD1145_sample_ids)] =  metadata_table['value'].loc[metadata_table['sample_id'].isin(POD1145_sample_ids)].replace("Ileum", "Ileum_allograft")
        
        sample_ids = metadata_table['sample_id'].unique()
        tissues = metadata_table['value'].loc[metadata_table['key'] == 'sample_origin'].unique()
                 
        for sample_id in sample_ids:
            tissues_per_sample[sample_id] = list(metadata_table['value'].loc[(metadata_table['key'] == 'sample_origin') & (metadata_table['sample_id'] == sample_id)])[0]
        for tissue in tissues:
            samples_per_tissue[tissue] = list(metadata_table['sample_id'].loc[metadata_table['value'] == tissue].unique())
            
    return tissues_per_sample, samples_per_tissue, tissues

# =============================================================================
# Make a dictionary containing every combination of tissues
# =============================================================================
def get_clones_in_tissue_pairs():
 
    for patient_ID, patient_name in subjects_dict.items():
    
        print("Collecting data for: " + patient_name)
    
        seq_table = pd.read_pickle(r"/path_to_sykes_children_data_folder/sequences_table/{}_sequences_table.pkl".format(patient_ID))
        seq_table.dropna(subset=["clone_id"], inplace=True)       
        seq_table['clone_id'] = seq_table['clone_id'].astype(int)
        
        tissues_per_sample, samples_per_tissue, tissues = get_tissue_sample_ids(patient_ID, patient_name)
        
        results_1 = {}
        for tissue1 in range(0, len(tissues) - 1):
            
            for tissue2 in range(tissue1 + 1, len(tissues)):
            
                results_1[tissues[tissue1] + "+" + tissues[tissue2]] = []
                   
        ###Make list of clones_ids to be analyzed
        sample_id_per_clone_df = seq_table.groupby(['clone_id'])['sample_id'].unique().reset_index()
        
        for clone_id, clone_samples in zip(sample_id_per_clone_df['clone_id'], sample_id_per_clone_df['sample_id']):
            
            tissues_in_clone = []
            for sample_id in clone_samples:
                tissue = tissues_per_sample[sample_id]
                if tissue not in tissues_in_clone:
                    tissues_in_clone.append(tissue)
                
            for tissue1 in range(0, len(tissues) - 1):                
                for tissue2 in range(tissue1 + 1, len(tissues)):          
                    if tissues[tissue1] in tissues_in_clone and tissues[tissue2] in tissues_in_clone:
                        results_1[tissues[tissue1] + "+" + tissues[tissue2]].append(clone_id)
                                    
        ###covert results lists to sets
        results_2 = {}
        for tissue_pair, clone_list in results_1.items():
            temp = set(clone_list)
            results_2[tissue_pair] = temp
            
        # =============================================================================
        # Save the result as a pickle    
        # =============================================================================
        path = Path(r"/path_to_sykes_children_data_folder/{}_clones_per_tissue_pairs_new_1.pkl".format(subjects_dict[patient_ID]))
        pickle_out = open(path,"wb")
        pickle.dump(results_2, pickle_out)
        pickle_out.close() 
                 
    print("Finished!")
                      
get_clones_in_tissue_pairs()

