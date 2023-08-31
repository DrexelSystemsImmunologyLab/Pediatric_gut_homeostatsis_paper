# -*- coding: utf-8 -*-
"""
Created on Thu Dec  5 14:20:09 2019

@author: thsia
"""

import pickle
import pandas as pd
import numpy as np
from numpy import *
import json
from pathlib import Path, PureWindowsPath
import subprocess

lp16_subjects_dict = {
                2: 'D145',
                3: 'D149',
                5: 'D168',
                6: 'D181',
                7: 'D182',
                8: 'D207', 
                    }
    
# =============================================================================
# Load tables          
# =============================================================================
def load_data(patient_name, patient_ID):

    if patient_name in lp16_subjects_dict.values():
        samples_table = pd.read_pickle(r"/path_to_lp16_data_folder/samples.pkl")
        metadata_table = pd.read_pickle(r"/path_to_lp16_data_folder/sample_metadata.pkl")           
        seq_table = pd.read_pickle(r"/path_to_lp16_data_folder/sequences_table/{}_sequences_table.pkl".format(patient_ID))
        seq_table.dropna(subset=["clone_id"], inplace=True)       
        seq_table = seq_table.loc[seq_table['functional'] == 1]
        seq_table['clone_id'] = seq_table['clone_id'].astype(int)
        ###Get functional clones using the clone_stats table
        clone_stats = pd.read_pickle(r"/path_to_lp16_data_folder/clone_stats/{}_clone_stats.pkl".format(patient_ID))
        functional_clones = clone_stats['clone_id'].loc[clone_stats['functional'] == 1].unique().tolist()
        seq_table = seq_table.loc[seq_table['clone_id'].isin(functional_clones)]                      
        subject_samples = samples_table['id'].loc[samples_table['subject_id'] == patient_ID].unique()   
        subject_metadata_table = metadata_table.loc[metadata_table['sample_id'].isin(subject_samples)]        
        tissues = subject_metadata_table['value'].loc[subject_metadata_table['key'] == 'tissue'].unique()   

    return seq_table, tissues

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


def get_data():
 
    pd.options.mode.chained_assignment = None  # Removing warning for: A value is trying to be set on a copy of a slice from a DataFrame.

    for patient_ID, patient_name in lp16_subjects_dict.items():
        result = {}          

        print("Collecting data for: " + patient_name)
        
        seq_table, tissues = load_data(patient_name, patient_ID)
        tissues_per_sample, samples_per_tissue = get_tissue_sample_ids(patient_name, patient_ID)  
        patient_tissues = list(samples_per_tissue.keys())
        patient_tissues = ['overall'] + patient_tissues
        
        for tissue in patient_tissues:
            
            print(tissue)
            if tissue == 'overall':
                df = seq_table[['clone_id', 'v_mutation_fraction']]    
            else:
                tissue_samples = samples_per_tissue[tissue]
                df = seq_table.loc[seq_table['sample_id'].isin(tissue_samples)]
                df = df[['clone_id', 'v_mutation_fraction']]    
                       
            df['percent_mutated'] = df['v_mutation_fraction'].apply(lambda x: x * 100)
            df_mutation_freq = df.groupby(['clone_id'])['percent_mutated'].mean().reset_index()
            df_mutation_freq = df_mutation_freq.dropna()
            
            clones = df['clone_id'].tolist()
            mutation_freq = df['percent_mutated'].tolist()
            
            tissue_result = dict(zip(clones, mutation_freq))
            result[tissue] = tissue_result

        print("Saving result...")

        # =============================================================================
        # Save the result as a pickle    
        # =============================================================================
        path = Path(r"/path_to_lp16_data_folder/{}_mutation_freq_per_clone_per_tissue_pretransplant_version.pkl".format(patient_name))
        pickle_out = open(path,"wb")
        pickle.dump(result, pickle_out)
        pickle_out.close() 
        
    print("Finished!")
           
get_data()    


