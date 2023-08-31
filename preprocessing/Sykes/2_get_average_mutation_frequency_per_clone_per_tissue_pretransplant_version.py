# -*- coding: utf-8 -*-
"""
Created on Thu Dec  5 14:20:09 2019

@author: thsia
"""

import pickle
import pandas as pd
from pathlib import Path, PureWindowsPath

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
def load_data(patient_name, patient_ID):
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
        
        ###Get functional clones using the clone_stats table
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

    return seq_table, tissues

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
        tissues = metadata_table['value'].loc[metadata_table['key'] == 'sample_origin'].unique()
                 
        for sample_id in sample_ids:
            tissues_per_sample[sample_id] = list(metadata_table['value'].loc[(metadata_table['key'] == 'sample_origin') & (metadata_table['sample_id'] == sample_id)])[0]
        for tissue in tissues:
            samples_per_tissue[tissue] = list(metadata_table['sample_id'].loc[metadata_table['value'] == tissue].unique())
      
    return tissues_per_sample, samples_per_tissue

def get_data():
 
    pd.options.mode.chained_assignment = None  # Removing warning for: A value is trying to be set on a copy of a slice from a DataFrame.

    for patient_ID, patient_name in subjects_dict.items():

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
        path = Path(r"/path_to_sykes_children_data_folder/{}_mutation_freq_per_clone_per_tissue_pretransplant_version.pkl".format(patient_name))
        pickle_out = open(path,"wb")
        pickle.dump(result, pickle_out)
        pickle_out.close() 
        
    print("Finished!")
           
get_data()    


