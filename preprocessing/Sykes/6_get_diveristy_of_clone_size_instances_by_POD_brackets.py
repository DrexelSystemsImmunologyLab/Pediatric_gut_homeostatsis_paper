# -*- coding: utf-8 -*-
"""
Created on Tue Sep 22 23:20:01 2020

@author: thsia

Note: this script requires the diversity program, which is a Haskell program, to be installed on the machine.
For more information see: https://github.com/GregorySchwartz/diversity

IMPORTANT:
    
    This script requires a temp folder to be made where files will be created and deleted. Make sure this folder is empty to start with because all files will be deleted in it.
    Change the following path throughout the script (it appears multiple times):
        
        /path_to_sykes_children_data_folder/temp/
        
    to whereever your temp folder is located.
    
"""


import json
import pickle
import pandas as pd
import numpy as np
from pathlib import Path, PureWindowsPath
import os
import re
import subprocess
import copy   
import sys

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
# Converts a list of sequences into fasta format.    
# =============================================================================
def convert_to_fasta(diversity_input):
    n = 0
    f = diversity_input
    output = r'/path_to_sykes_children_data_folder/temp/diversity/input_for_diversity.fasta'
    with open(output, 'w') as out:
        for seq in f:
            n += 1
            out.write('\n>item ' + str(n) + "" + '\n' + seq)        
 
    
# =============================================================================
# This function accepts a list of sequences and outputs the diversity in a separate list, one for order, window, position, weight, diversity.
# e.g. [seq1, seq2, seq3, seq4]
# Note: the position of the seq starts at 1. For example, if your first nucleotide only appears at position 3, your seq needs two dashes to fill 1 and 2: --nucleotide
# =============================================================================
def get_diversity(diversity_input, diversity_order):
    
    ###change directory to the directly where results will be saved
    os.chdir(r"/path_to_sykes_children_data_folder/temp/diversity") 
    subprocess.run(r'find . -type f -delete', shell=True, check=True, stdout=subprocess.PIPE, universal_newlines=True)

    convert_to_fasta(diversity_input)
    
    path, dirs, files = next(os.walk('/path_to_sykes_children_data_folder/temp/diversity'))
    
    # file_count = len(files)    
    # print("Calculating diversity. There are {} files.".format(file_count))
    # =============================================================================
    # This command runs the diversity haskell program from the terminal on the file produced above        
    # =============================================================================
    file_name_linux_syntax = r'${file%%.*}'
    str_command = r'for file in *.fasta; do cat {}.fasta | diversity -i {}.fasta -o {}.csv -r {} -L -a; done'.format(file_name_linux_syntax, file_name_linux_syntax, file_name_linux_syntax, diversity_order)    
    subprocess.run(str_command, shell=True, check=True, stdout=subprocess.PIPE, universal_newlines=True)

    order = []
    window = []
    position = []
    weight = []
    diversity = []
        
    
    # =============================================================================
    # Results are in the following order: label,order,window,position,weight,diversity    
    # =============================================================================
    for entry in os.scandir(r"/path_to_sykes_children_data_folder/temp/diversity"):    
        if entry.path.endswith(".csv") and entry.is_file():
            with open(entry) as reader:    
                ###skip headers using [1:]
                for line in reader.readlines()[1:]:  
                    # print(line)
                    ###remove the \n
                    line = line.rsplit()[0]
                    line = line.rsplit(",")
                    
                    order.append(line[1])
                    window.append(line[2])
                    position.append(line[3])
                    weight.append(line[4])
                    diversity.append(line[5])
                    
    return order, window, position, weight, diversity
    
# =============================================================================
# Load tables          
# =============================================================================
def load_data_convert_ileum_to_ileum_allograft(patient_ID, patient_name):
    if patient_name in subjects_dict.values():
        samples_table = pd.read_pickle(r"/path_to_sykes_children_data_folder/samples.pkl")
        metadata_table = pd.read_pickle(r"/path_to_sykes_children_data_folder/sample_metadata.pkl")    
        metadata_table['value'] = metadata_table['value'].replace("Colon ", "Colon")
        metadata_table['value'] = metadata_table['value'].replace("Ileum", "Ileum_allograft") #Convert to Ileum_allograft so it's considered a single tissue
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
        tissues = subject_metadata_table['value'].loc[subject_metadata_table['key'] == 'sample_origin'].unique()

    return seq_table, tissues, rejection_samples
        

# =============================================================================
# POD grouping function to group PODs together that are within 2 days of each other
# Input is name of patient, such as "Pt19_R"
# Output is list of lists, with each sublist containing the grouped PODs (or single POD if there are no PODs that are grouped with it)
# And a second output in the same format, except each sublist is a list of sample_ids that belong to the corresponding POD
# Example: [[0], [20, 22], [100], [110, 112, 114]], [[sample_ids that belong to 0], [sample_ids that belong to 20, 22], [etc], [etc]]
# =============================================================================
def group_PODs(patient_ID, patient_name):
    samples_table = pd.read_pickle(r"/path_to_sykes_children_data_folder/samples.pkl")
    metadata_table = pd.read_pickle(r"/path_to_sykes_children_data_folder/sample_metadata.pkl")
    metadata_table['value'] = metadata_table['value'].replace("POD37 44", 'POD37')
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
    subject_samples = samples_table['id'].loc[samples_table['subject_id'] == patient_ID].unique()
    metadata_table = metadata_table.loc[metadata_table['sample_id'].isin(subject_samples)]
    pod_list_strings = metadata_table['value'].loc[metadata_table['key'] == 'pod'].unique()
    pod_list = [int(''.join(i for i in x if i.isdigit())) for x in pod_list_strings]     
    ###sort both using the int version as the guide
    pod_list, pod_list_strings = (list(t) for t in zip(*sorted(zip(pod_list, pod_list_strings))))
    
    result = []
    result_samples = []
    current_pod = pod_list[0]
    current_samples = metadata_table['sample_id'].loc[metadata_table['value'] == pod_list_strings[pod_list.index(current_pod)]].unique().tolist()
    temp = []
    temp_samples = []
    for next_pod in pod_list[1:]: 
        next_samples = metadata_table['sample_id'].loc[metadata_table['value'] == pod_list_strings[pod_list.index(next_pod)]].unique().tolist()
        if len(temp) == 0:
            temp.append(current_pod)
            for sample in current_samples:
                temp_samples.append(sample)
        if next_pod - current_pod <= 2:       
            temp.append(next_pod)
            for sample in next_samples:
                temp_samples.append(sample)                
        else:
            result.append(temp)
            result_samples.append(temp_samples)
            temp_samples = []
            temp = []            
        current_pod = next_pod
        current_samples = next_samples
    if len(pod_list) > 1:
        if pod_list[-1] not in result[-1]:
            result.append([pod_list[-1]])
            result_samples.append(metadata_table['sample_id'].loc[metadata_table['value'] == pod_list_strings[pod_list.index(pod_list[-1])]].unique().tolist())

    return result, result_samples

# =============================================================================
# This function returns two dicts:
# samples_per_tissue: dict[tissue] : [list of sample_ids in that tissue]    
# tissues_per_sample: dict[sample_id] : [tissue in that sample_id] ###this will be a list with a single element (the tissue)
# =============================================================================
def get_tissue_sample_ids_convert_ileum_to_ileum_allograft(patient_ID, patient_name):
  
    tissues_per_sample = {}
    samples_per_tissue = {}
    
    if patient_name in subjects_dict.values():
        samples_table = pd.read_pickle(r"/path_to_sykes_children_data_folder/samples.pkl")
        metadata_table = pd.read_pickle(r"/path_to_sykes_children_data_folder/sample_metadata.pkl")
        
        subject_samples = samples_table['id'].loc[samples_table['subject_id'] == patient_ID].unique() 
        
        metadata_table = metadata_table.loc[metadata_table['sample_id'].isin(subject_samples)]
        metadata_table['value'] = metadata_table['value'].replace("Colon ", "Colon")
        metadata_table['value'] = metadata_table['value'].replace("Ileum", "Ileum_allograft") #Convert to Ileum_allograft so it's considered a single tissue
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
        tissues = metadata_table['value'].loc[metadata_table['key'] == 'sample_origin'].unique()
                 
        for sample_id in sample_ids:
            tissues_per_sample[sample_id] = list(metadata_table['value'].loc[(metadata_table['key'] == 'sample_origin') & (metadata_table['sample_id'] == sample_id)])[0]
        for tissue in tissues:
            samples_per_tissue[tissue] = list(metadata_table['sample_id'].loc[metadata_table['value'] == tissue].unique())
            
    return tissues_per_sample, samples_per_tissue

def get_timepoint_sample_ids_convert_ileum_to_ileum_allograft(patient_ID, patient_name):
  
    timepoints_per_sample = {}
    samples_per_timepoint = {}
    
    if patient_name in subjects_dict.values():
        samples_table = pd.read_pickle(r"/path_to_sykes_children_data_folder/samples.pkl")
        metadata_table = pd.read_pickle(r"/path_to_sykes_children_data_folder/sample_metadata.pkl")     
        subject_samples = samples_table['id'].loc[samples_table['subject_id'] == patient_ID].unique()        
        metadata_table = metadata_table.loc[metadata_table['sample_id'].isin(subject_samples)]
        metadata_table['value'] = metadata_table['value'].replace("Colon ", "Colon")
        metadata_table['value'] = metadata_table['value'].replace("Ileum", "Ileum_allograft") #Convert to Ileum_allograft so it's considered a single tissue
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

                
def collect_diversity_per_tissue():
    
    pod_list = []
    tissue_list = []
    diversity_result_list = []   
    patient_list = [] 
    diversity_order_list = [] 
    total_clones_list = []    
    for patient_ID, patient_name in subjects_dict.items():
        print("Collecting data for: " + patient_name)
    
        seq_table, tissues, rejection_samples = load_data_convert_ileum_to_ileum_allograft(patient_ID, patient_name)
        timepoints_per_sample, samples_per_timepoint = get_timepoint_sample_ids_convert_ileum_to_ileum_allograft(patient_ID, patient_name)
        tissues_per_sample, samples_per_tissue = get_tissue_sample_ids_convert_ileum_to_ileum_allograft(patient_ID, patient_name)
    
        df = seq_table[['clone_id', 'sample_id', 'ai']]
        
        patient_pod_list = [int(x) for x in samples_per_timepoint.keys()]
        patient_pod_list = sorted(patient_pod_list)
        patient_pod_list_strings = [str(x) for x in patient_pod_list]   
        
        for pod_bracket in ['0', '1-90', '91-365', '>365']:
            if "-" in pod_bracket:
                pod_min = int(pod_bracket.split("-")[0])
                pod_max = int(pod_bracket.split("-")[1])
            elif pod_bracket == '0':
                pod_min = 0
                pod_max = 0
            else:
                pod_min = 366
                pod_max = 9999999
            samples_in_bracket = []
            for pod_int, pod_string in zip(patient_pod_list, patient_pod_list_strings):
                if pod_int >= pod_min and pod_int <= pod_max:
                    pod_samples = samples_per_timepoint[pod_string]
                    samples_in_bracket += pod_samples
            pod_df = df.loc[df['sample_id'].isin(samples_in_bracket)]
                        
            print("On pod bracke: " + pod_bracket)
            for tissue in tissues:
                tissue_samples = samples_per_tissue[tissue]
                tissue_df = pod_df.loc[pod_df['sample_id'].isin(tissue_samples)]
                total_clones = tissue_df['clone_id'].nunique()
    
                if total_clones > 0:  
                        
                    grouped_df = tissue_df.groupby('clone_id')['ai'].nunique().reset_index()
                    # =============================================================================
                    # Build a list of clones that will be used to calculate the diversity.
                    # The number of times a clone appears is equal to the number of instances it has (i.e., the ['ai'] count)
                    # =============================================================================
                    diversity_input = []
                    for clone_id, num_instances in zip(grouped_df['clone_id'], grouped_df['ai']):
                        for x in range(num_instances):
                            diversity_input.append(str(clone_id)) #Names must be in string format for the diversity program
                            
                    for diversity_order in [0, 1, 2, 3, 4, 5, 6]:
                        order, window, position, weight, diversity = get_diversity(diversity_input, diversity_order)
                        order = int(float(order[0]))
                        diversity = float(diversity[0])
                        
                        pod_list.append(pod_bracket)
                        tissue_list.append(tissue)
                        diversity_result_list.append(diversity)
                        patient_list.append(patient_name)
                        diversity_order_list.append(order)   
                        total_clones_list.append(total_clones)
                            
    result_df = pd.DataFrame()
    result_df['patient'] = patient_list
    result_df['pod'] = pod_list
    result_df['tissue'] = tissue_list
    result_df['diversity'] = diversity_result_list
    result_df['diversity_order'] = diversity_order_list   
    result_df['total_clones'] = total_clones_list    
    
    # =============================================================================
    # Save result
    # =============================================================================
    path = Path(r"/path_to_sykes_children_data_folder/Repertoire_diversity_by_POD_BRACKETS_and_tissue_and_diversity_order.pkl")
    pickle_out = open(path,"wb")
    pickle.dump(result_df, pickle_out)
    pickle_out.close() 

collect_diversity_per_tissue()






