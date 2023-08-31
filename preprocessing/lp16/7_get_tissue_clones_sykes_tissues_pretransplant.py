# -*- coding: utf-8 -*-
"""
Created on Thu Dec  5 14:20:09 2019

@author: thsia


"""

import pandas as pd
import pickle
import numpy as np
from pathlib import Path, PureWindowsPath
from numpy import *

lp16_subjects_dict = {
                2: 'D145',
                3: 'D149',
                5: 'D168',
                6: 'D181',
                7: 'D182',
                8: 'D207', 
                    }

def get_tissue_clones():
    
    # =============================================================================
    # This function returns two dicts:
    # samples_per_tissue: dict[tissue] : [list of sample_ids in that tissue]    
    # tissues_per_sample: dict[sample_id] : [tissue in that sample_id] ###this will be a list with a single element (the tissue)
    # =============================================================================
    def get_tissue_sample_ids(patient_name, patient_ID):
      
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
                
        return tissues_per_sample, samples_per_tissue
                    
    for patient_ID, patient_name in lp16_subjects_dict.items():
        print("Loading data for: " + lp16_subjects_dict[patient_ID])
        df = pd.read_pickle(r"/path_to_lp16_data_folder/sequences_table/{}_sequences_table.pkl".format(patient_ID))
        df.dropna(subset=["clone_id"], inplace=True)
        df = df.loc[df['functional'] == 1]
        df['clone_id'] = df['clone_id'].astype(int)

        print("Getting tissue clones for: " + lp16_subjects_dict[patient_ID])
        def get_tissue_clones(df, patient_ID, patient_name):
            result = {}
            ###get sample_ids for each tissue
            tissues_per_sample, samples_per_tissue = get_tissue_sample_ids(patient_name, patient_ID)
            
            tissues = list(samples_per_tissue.keys())
            for tissue in tissues:
                
                tissue_samples = samples_per_tissue[tissue]
                clones_in_tissue = set(df['clone_id'].loc[df['sample_id'].isin(tissue_samples)].tolist())
                result[tissue] = clones_in_tissue

            return result
        
        tissue_clones = get_tissue_clones(df, patient_ID, patient_name)
        
        for tissue, value in tissue_clones.items():
            print("There are {} clones in {}.".format(len(value), tissue))
        
        print("Finished with: " + patient_name)
        # =============================================================================
        # Save result
        # =============================================================================
        path = Path(r"/path_to_lp16_data_folder/{}_tissue_clones_pretransplant.pkl".format(patient_name))
         
        pickle_out = open(path,"wb")
        pickle.dump(tissue_clones, pickle_out)
        pickle_out.close()

    print("Finished!")
        
get_tissue_clones()










