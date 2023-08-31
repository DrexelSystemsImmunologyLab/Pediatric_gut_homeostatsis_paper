"""

"""

import json
import pickle
import pandas as pd
from pathlib import Path, PureWindowsPath
import os

lp16_subjects_dict = {
                2: 'D145',
                3: 'D149',
                5: 'D168',
                6: 'D181',
                7: 'D182',
                8: 'D207', 
                    }
          
# =============================================================================
# Make a dictionary containing every combination of tissues
# =============================================================================
def get_clones_in_tissue_pairs():
 
    # =============================================================================
    # This function returns two dicts:
    # samples_per_tissue: dict[tissue] : [list of sample_ids in that tissue]    
    # tissues_per_sample: dict[sample_id] : [tissue in that sample_id] ###this will be a list with a single element (the tissue)
    # =============================================================================
    def get_tissue_sample_ids(patient_name):
      
        tissues_per_sample = {}
        samples_per_tissue = {}
        
        if patient_name in lp16_subjects_dict.values():
            samples_table = pd.read_pickle(r"/path_to_lp16_data_folder/samples.pkl")
            metadata_table = pd.read_pickle(r"/path_to_lp16_data_folder/sample_metadata.pkl")
            subject_samples = samples_table['id'].loc[samples_table['subject_id'] == patient_ID].unique() 
            metadata_table = metadata_table.loc[metadata_table['sample_id'].isin(subject_samples)]
            sample_ids = metadata_table['sample_id'].unique()
            tissues = metadata_table['value'].loc[metadata_table['key'] == 'tissue'].unique()
                     
            for sample_id in sample_ids:
                tissues_per_sample[sample_id] = list(metadata_table['value'].loc[(metadata_table['key'] == 'tissue') & (metadata_table['sample_id'] == sample_id)])[0]
            for tissue in tissues:
                samples_per_tissue[tissue] = list(metadata_table['sample_id'].loc[metadata_table['value'] == tissue].unique())
                
        return tissues_per_sample, samples_per_tissue, tissues
    
    for patient_ID, patient_name in lp16_subjects_dict.items():
    
        print("Making list of clones for: " + patient_name)
    
        seq_table = pd.read_pickle(r"/path_to_lp16_data_folder/sequences_table/{}_sequences_table.pkl".format(patient_ID))
        seq_table.dropna(subset=["clone_id"], inplace=True)       
        seq_table['clone_id'] = seq_table['clone_id'].astype(int)
        tissues_per_sample, samples_per_tissue, tissues = get_tissue_sample_ids(patient_name)
        
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
        path = Path(r"/path_to_lp16_data_folder/{}_clones_per_tissue_pairs_new_1.pkl".format(patient_name))
        pickle_out = open(path,"wb")
        pickle.dump(results_2, pickle_out)
        pickle_out.close() 
       
    print("Finished!")
                      
get_clones_in_tissue_pairs()

