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

def get_tissue_clones():

    # =============================================================================
    # This function returns two dicts:
    # samples_per_tissue: dict[tissue] : [list of sample_ids in that tissue]    
    # tissues_per_sample: dict[sample_id] : [tissue in that sample_id] ###this will be a list with a single element (the tissue)
    # =============================================================================
    def get_tissue_sample_ids(patient_name, patient_ID):
      
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
                        
            # =============================================================================
            # Add "_pretransplant" to tissues at time zero        
            # =============================================================================
            pretransplant_samples = metadata_table['sample_id'].loc[metadata_table['value'] == 'POD0'].tolist()
            pretrans_tissues = metadata_table['value'].loc[(metadata_table['key'] == 'sample_origin') & (metadata_table['sample_id'].isin(pretransplant_samples))]
            for tissue in pretrans_tissues:         
                metadata_table['value'].loc[(metadata_table['value'] == tissue) & (metadata_table['sample_id'].isin(pretransplant_samples))] = tissue + "_pretransplant"
            
            sample_ids = metadata_table['sample_id'].unique()
            tissues = metadata_table['value'].loc[metadata_table['key'] == 'sample_origin'].unique()
                     
            for sample_id in sample_ids:
                tissues_per_sample[sample_id] = list(metadata_table['value'].loc[(metadata_table['key'] == 'sample_origin') & (metadata_table['sample_id'] == sample_id)])[0]
            for tissue in tissues:
                samples_per_tissue[tissue] = list(metadata_table['sample_id'].loc[metadata_table['value'] == tissue].unique())

            return tissues_per_sample, samples_per_tissue
                    
    for patient_ID, patient_name in subjects_dict.items():
        print("Loading data for: " + subjects_dict[patient_ID])
        df = pd.read_pickle(r"/path_to_sykes_children_data_folder/sequences_table/{}_sequences_table.pkl".format(patient_ID))
        df.dropna(subset=["clone_id"], inplace=True)
        df = df.loc[df['functional'] == 1]
        df['clone_id'] = df['clone_id'].astype(int)

        print("Getting tissue clones for: " + subjects_dict[patient_ID])
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
        path = Path(r"/path_to_sykes_children_data_folder/{}_tissue_clones_pretransplant.pkl".format(patient_name))
         
        pickle_out = open(path,"wb")
        pickle.dump(tissue_clones, pickle_out)
        pickle_out.close()

    print("Finished!")
        
get_tissue_clones()
