# -*- coding: utf-8 -*-
"""
Created on Thu Dec  5 14:20:09 2019

@author: thsia
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
# Load tables          
# =============================================================================
def load_data_pretransplant_version(patient_name, patient_ID):
    if patient_name in subjects_dict.values():
        samples_table = pd.read_pickle(r"/path_to_sykes_children_data_folder/samples.pkl")
        metadata_table = pd.read_pickle(r"/path_to_sykes_children_data_folder/sample_metadata.pkl")    
        metadata_table['value'] = metadata_table['value'].replace("Colon ", "Colon")
        metadata_table['value'] = metadata_table['value'].replace("MLN_allograft", "AxLN")
        metadata_table['value'] = metadata_table['value'].replace("POD37 44", "POD37")
        ###Change Ileum to Ileum_allograft for Pt21_R POD1145
        POD1145_sample_ids = metadata_table['sample_id'].loc[metadata_table['value'] == 'POD1145']
        metadata_table['value'].loc[metadata_table['sample_id'].isin(POD1145_sample_ids)] =  metadata_table['value'].loc[metadata_table['sample_id'].isin(POD1145_sample_ids)].replace("Ileum", "Ileum_allograft")
        seq_table = pd.read_pickle(r"/path_to_sykes_children_data_folder/sequences_table/{}_sequences_table.pkl".format(patient_ID))
        seq_table.dropna(subset=["clone_id"], inplace=True)       
        seq_table = seq_table.loc[seq_table['functional'] == 1]
        seq_table['clone_id'] = seq_table['clone_id'].astype(int)
        
        # =============================================================================
        # Load rejection samples and remove these samples from data
        # =============================================================================
        path = Path(r"/path_to_sykes_children_data_folder/{}_samples_with_rejection.pkl".format(patient_name))
        pickle_in = open(path,"rb")
        rejection_samples = pickle.load(pickle_in)
        pickle_in.close()
        
        seq_table = seq_table.loc[~seq_table['sample_id'].isin(rejection_samples)]
        metadata_table = metadata_table.loc[~metadata_table['sample_id'].isin(rejection_samples)]
        samples_table = samples_table.loc[~samples_table['id'].isin(rejection_samples)]
        clone_stats = pd.read_pickle(r"/path_to_sykes_children_data_folder/clone_stats/{}_clone_stats.pkl".format(patient_ID))
        functional_clones = clone_stats['clone_id'].loc[clone_stats['functional'] == 1].unique().tolist()
        seq_table = seq_table.loc[seq_table['clone_id'].isin(functional_clones)]
        subject_samples = samples_table['id'].loc[samples_table['subject_id'] == patient_ID].unique()   
        subject_metadata_table = metadata_table.loc[metadata_table['sample_id'].isin(subject_samples)]     
        
        # =============================================================================
        # Add "_pretransplant" to tissues at time zero        
        # =============================================================================
        pretransplant_samples = subject_metadata_table['sample_id'].loc[subject_metadata_table['value'] == 'POD0'].tolist()
        pretrans_tissues = subject_metadata_table['value'].loc[(subject_metadata_table['key'] == 'sample_origin') & (subject_metadata_table['sample_id'].isin(pretransplant_samples))]
        for tissue in pretrans_tissues:         
            subject_metadata_table['value'].loc[(subject_metadata_table['value'] == tissue) & (subject_metadata_table['sample_id'].isin(pretransplant_samples))] = tissue + "_pretransplant"
                             
        tissues = subject_metadata_table['value'].loc[subject_metadata_table['key'] == 'sample_origin'].unique()

    return tissues, seq_table
         
# =============================================================================
# This function returns two dicts:
# samples_per_tissue: dict[tissue] : [list of sample_ids in that tissue]    
# tissues_per_sample: dict[sample_id] : [tissue in that sample_id] ###this will be a list with a single element (the tissue)
# =============================================================================
def get_tissue_sample_ids_pretransplant_version(patient_name, patient_ID):
  
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
        # Load rejection samples and remove these samples from data
        # =============================================================================
        path = Path(r"/path_to_sykes_children_data_folder/{}_samples_with_rejection.pkl".format(patient_name))
        pickle_in = open(path,"rb")
        rejection_samples = pickle.load(pickle_in)
        pickle_in.close()
        metadata_table = metadata_table.loc[~metadata_table['sample_id'].isin(rejection_samples)]
        samples_table = samples_table.loc[~samples_table['id'].isin(rejection_samples)]
                    
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
                 
def get_timepoint_sample_ids_pretransplant_version(patient_name, patient_ID):
  
    timepoints_per_sample = {}
    samples_per_timepoint = {}
    
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
        # Load rejection samples and remove these samples from data
        # =============================================================================
        path = Path(r"/path_to_sykes_children_data_folder/{}_samples_with_rejection.pkl".format(patient_name))
        pickle_in = open(path,"rb")
        rejection_samples = pickle.load(pickle_in)
        pickle_in.close()
        metadata_table = metadata_table.loc[~metadata_table['sample_id'].isin(rejection_samples)]
        samples_table = samples_table.loc[~samples_table['id'].isin(rejection_samples)]

        sample_ids = metadata_table['sample_id'].unique()
        timepoints = metadata_table['value'].loc[metadata_table['key'] == 'pod'].unique().tolist()
        timepoints_original = metadata_table['value'].loc[metadata_table['key'] == 'pod'].unique().tolist()
        timepoints = [int(''.join(i for i in x if i.isdigit())) for x in timepoints]
        timepoints, timepoints_original = (list(t) for t in zip(*sorted(zip(timepoints, timepoints_original))))
        
        for sample_id in sample_ids:
            result = list(metadata_table['value'].loc[(metadata_table['key'] == 'pod') & (metadata_table['sample_id'] == sample_id)])[0]
            result = int(''.join(i for i in result if i.isdigit()))
            timepoints_per_sample[sample_id] = result                
        for timepoint, timepoint_original in zip(timepoints, timepoints_original):
            samples_per_timepoint[str(timepoint)] = list(metadata_table['sample_id'].loc[metadata_table['value'] == timepoint_original].unique())
            
    return timepoints_per_sample, samples_per_timepoint
              
def get_mutation_freq_per_timepoint():

    pod_list = []
    tissue_list = []
    mutation_result = []   
    patient_list = [] 
    total_clones_list = []
    for patient_ID, patient_name in subjects_dict.items():
        print("Collecting data for: " + patient_name)
    
        tissues, seq_table = load_data_pretransplant_version(patient_name, patient_ID)
        tissues_per_sample, samples_per_tissue = get_tissue_sample_ids_pretransplant_version(patient_name, patient_ID)        
        pod_per_sample, samples_per_timepoint = get_timepoint_sample_ids_pretransplant_version(patient_name, patient_ID) 
                     
        df = seq_table[['clone_id', 'sample_id', 'v_mutation_fraction']]
        
        patient_pod_list = [int(x) for x in samples_per_timepoint.keys()]
        patient_pod_list = sorted(patient_pod_list)
        patient_pod_list_strings = [str(x) for x in patient_pod_list]      
        
        for pod_int, pod_string in zip(patient_pod_list, patient_pod_list_strings):
            pod_samples = samples_per_timepoint[pod_string]
            pod_df = df.loc[df['sample_id'].isin(pod_samples)]
            for tissue in tissues:
                tissue_samples = samples_per_tissue[tissue]
                tissue_df = pod_df.loc[pod_df['sample_id'].isin(tissue_samples)]
                tissue_clones = set(tissue_df['clone_id'].unique().tolist())
                total_clones = len(tissue_clones)
    
                if total_clones > 0:  
                    
                    tissue_df['percent_mutated'] = tissue_df['v_mutation_fraction'].apply(lambda x: x * 100)
                    df_mutation_freq = tissue_df.groupby(['clone_id'])['percent_mutated'].mean().reset_index()
                    df_mutation_freq = df_mutation_freq.dropna()
                    df_mutation_freq = df_mutation_freq.loc[df_mutation_freq['percent_mutated'] > 2]
                    num_clones = df_mutation_freq['clone_id'].nunique()
                    
                    fraction_greater_than_2 = num_clones / total_clones
                    
                    mutation_result.append(fraction_greater_than_2)                
                    pod_list.append(pod_int)
                    tissue_list.append(tissue)
                    patient_list.append(patient_name)
                    total_clones_list.append(total_clones)

    result_df = pd.DataFrame()
    result_df['pod'] = pod_list
    result_df['tissue'] = tissue_list
    result_df['percent_mutation_freq_greater_than_2'] = mutation_result
    result_df['patient'] = patient_list
    result_df['total_clones'] = total_clones_list
    
    # =============================================================================
    # Save result
    # =============================================================================
    path = Path(r"/path_to_sykes_children_data_folder/Percent_mutation_freq_greater_than_2_by_timepoint_and_node_min.pkl")
    pickle_out = open(path,"wb")
    pickle.dump(result_df, pickle_out)
    pickle_out.close()
      
    print("Finished!")
        
get_mutation_freq_per_timepoint()
    

 
