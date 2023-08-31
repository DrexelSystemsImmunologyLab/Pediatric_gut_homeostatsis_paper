# -*- coding: utf-8 -*-
"""
This script is performing the following simple calculation for calculating a "trunk":
    Find the top most-frequent mutations using the clone_stats table.
    for trunks with 5 mutations, simply take the 5th mutation (the least frequent of the top 5) and if it's in 85%+ of the nodes, then count it as a trunk with 5.
    The percent of unique sequences can be easily calculated from the last row in the clone stats table. The number of unique sequences for each mutation is indicated by "unique" in the mutation's dictionary.
    The total number of unique sequences in the clone is shown in the last row in the clone_stats table under "unique_cnt", so if the mutation's unqiue/unique_count >= 85%, then the clone has a trunk.'
    If the 5th mutation is not in 85%+, check the 4th, 3rd, etc.
    If a mutation is found in 85%+, then the remaining most-frequent mutations will automatically also be in 85%+.


@author: thsia
"""

import json
import pickle
import pandas as pd
import numpy as np
from pathlib import Path, PureWindowsPath
import os
import re
import statistics
import subprocess
import copy

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
        df_of_trees = pd.read_pickle(r"/path_to_sykes_children_data_folder/clones_table/{}_clones_table.pkl".format(patient_ID))
        seq_table = pd.read_pickle(r"/path_to_sykes_children_data_folder/sequences_table/{}_sequences_table.pkl".format(patient_ID))
        seq_table.dropna(subset=["clone_id"], inplace=True)       
        seq_table['clone_id'] = seq_table['clone_id'].astype(int)
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
        tissues = subject_metadata_table['value'].loc[subject_metadata_table['key'] == 'sample_origin'].unique().tolist()          
        
    return seq_table, tissues, clone_stats, df_of_trees
         
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
        # Add "_pretransplant" to tissues at time zero        
        # =============================================================================
        pretransplant_samples = metadata_table['sample_id'].loc[metadata_table['value'] == 'POD0'].tolist()
        pretrans_tissues = metadata_table['value'].loc[(metadata_table['key'] == 'sample_origin') & (metadata_table['sample_id'].isin(pretransplant_samples))]
        for tissue in pretrans_tissues:         
            metadata_table['value'].loc[(metadata_table['value'] == tissue) & (metadata_table['sample_id'].isin(pretransplant_samples))] = tissue + "_pretransplant"
        
        sample_ids = metadata_table['sample_id'].unique()
        tissues = metadata_table['value'].loc[metadata_table['key'] == 'sample_origin'].unique().tolist()
                 
        for sample_id in sample_ids:
            tissues_per_sample[sample_id] = list(metadata_table['value'].loc[(metadata_table['key'] == 'sample_origin') & (metadata_table['sample_id'] == sample_id)])[0]
        for tissue in tissues:
            samples_per_tissue[tissue] = list(metadata_table['sample_id'].loc[metadata_table['value'] == tissue].unique())
            
    return tissues_per_sample, samples_per_tissue
        
def check_if_unmutated(json_tree):
    clone = json.loads(json_tree)
    nodes_by_mutation = {}
    node_counter = {'count': 0}
    all_unique_mutations = []
    def traverse_tree(tree, node_counter, nodes_by_mutation, ancestor_mutations):
        cur_node_num = node_counter['count']
        cur_mutations = copy.deepcopy(ancestor_mutations)
        for mutation in tree['data']['mutations']:
            mutation_name = mutation['from'] + ' ' + str(mutation['pos']) + ' ' + mutation['to']
            if cur_node_num not in nodes_by_mutation.keys():
                nodes_by_mutation[cur_node_num] = []
            if len(tree['data']['seq_ids']) > 0:
                nodes_by_mutation[cur_node_num].append(mutation_name)
            else:
                ###remove nodes with no sequences
                del nodes_by_mutation[cur_node_num]
            cur_mutations.append(mutation_name)
            if mutation_name not in all_unique_mutations:
                all_unique_mutations.append(mutation_name)
        for mutation_name in ancestor_mutations:
            if len(tree['data']['seq_ids']) > 0:
                nodes_by_mutation[cur_node_num].append(mutation_name) 
        for child in tree['children']:
            node_counter['count'] += 1
            traverse_tree(child, node_counter, nodes_by_mutation, cur_mutations)       
    traverse_tree(clone['tree'], node_counter, nodes_by_mutation, [])                   

    itemsets = []
    for node, mutations in nodes_by_mutation.items():
        itemsets.append(mutations)
    num_nodes = len(itemsets)
    itemsets_with_at_least_5 = 0                   
    for itemset in itemsets:
        temp = len(itemset)
        if temp >= 5:
            itemsets_with_at_least_5 += 1
            
    if itemsets_with_at_least_5 < num_nodes * 0.85: 
        result = True
    else:
        result = False
    return result
        

def calculate_trunk_length():

    # =============================================================================
    # Body of function starts here
    # =============================================================================
    master_result = {} 
    for patient_ID, patient_name in subjects_dict.items():
        
        print("Currently analyzing " + patient_name)
        
        seq_table, tissues, clone_stats, clones_table = load_data_pretransplant_version(patient_name, patient_ID)
        tissues_per_sample, samples_per_tissue = get_tissue_sample_ids_pretransplant_version(patient_name, patient_ID)
        
        clone_stats = clone_stats.loc[np.isnan(clone_stats['sample_id'])]     
        clone_stats = clone_stats.merge(clones_table[['id', 'tree']], left_on="clone_id", right_on="id", how="left")

        results = {}
        results['unmutated'] = {}
        results['pretrunk'] = {}
        results['trunk'] = {}
                    
        clone_counter = 0
        for clone_id, clone_mutations, unique_cnt, tree in zip(clone_stats['clone_id'], clone_stats['mutations'], clone_stats['unique_cnt'], clone_stats['tree']):
            clone_counter += 1
     
            find_freq_mutations = {}
            for mutations in [clone_mutations]:
                mutations = json.loads(mutations)
                if len(mutations['regions']) > 0:
                    ALL_mutations = mutations['regions']['ALL']
                    for mutation_type, mutation_values in ALL_mutations.items():
                        for individual_mutations in mutation_values:
                            ###position gets +1 because the database adds +1 to the position so that it starts at 1 and not at 0
                            mutation_name = individual_mutations['from_nt'] + " " + str(individual_mutations['pos'] + 1) + " " + individual_mutations['to_nt']
                            num_unique_seqs_with_mutation = individual_mutations['unique']

                            if mutation_name not in find_freq_mutations.keys():
                                find_freq_mutations[mutation_name] = num_unique_seqs_with_mutation
                            else:
                                find_freq_mutations[mutation_name] += num_unique_seqs_with_mutation
            clone_num_nodes = unique_cnt
                  
            most_freq_mutations = {k: v for k, v in sorted(find_freq_mutations.items(), key=lambda item: item[1], reverse=True)}  

            num_nodes_for_85percent_threshold = clone_num_nodes * .85
            
            mutations_result = []
            fraction_of_nodes_result = []
            node_count_result = []
            for mutation_name, num_nodes in most_freq_mutations.items():
                if num_nodes >= num_nodes_for_85percent_threshold:
                    mutations_result.append(mutation_name)
                    fraction = round(num_nodes / clone_num_nodes, 2)
                    fraction_of_nodes_result.append(fraction)
                    node_count_result.append(num_nodes)
            num_mutations_in_trunk = len(mutations_result)
            
            temp_clone_result = {
                                'mutations' : mutations_result, 
                                'node_cnt' : node_count_result,
                                'fraction_of_nodes' : fraction_of_nodes_result,
                                'trunk_length' : num_mutations_in_trunk,                                 
                                }        

            is_clone_unmutated = check_if_unmutated(tree)
            if is_clone_unmutated == True:
                results['unmutated'][clone_id] = temp_clone_result
                
            elif num_mutations_in_trunk >= 5:
                results['trunk'][clone_id] = temp_clone_result
                                  
            else:
                results['pretrunk'][clone_id] = temp_clone_result
                             
        print('unmutated:', len(results['unmutated']))
        print('pretrunk:', len(results['pretrunk']))
        print('trunk:', len(results['trunk']))
            
        master_result[patient_name] = results
                
    # =============================================================================
    # Save result
    # =============================================================================
    path = Path(r"/path_to_sykes_children_data_folder/trunk_length_simple_method_all_clones_with_tissues_combined.pkl")
     
    pickle_out = open(path,"wb")
    pickle.dump(master_result, pickle_out)
    pickle_out.close() 
               
calculate_trunk_length()   
