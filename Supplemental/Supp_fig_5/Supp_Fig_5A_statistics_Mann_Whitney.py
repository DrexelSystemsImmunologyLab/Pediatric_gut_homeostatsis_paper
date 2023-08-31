# -*- coding: utf-8 -*-
"""
Created on Thu Jul  2 13:24:30 2020

@author: thsia
"""

import pandas as pd
import numpy as np
from math import hypot
from itertools import combinations
import pickle
from numpy import *
from pathlib import Path, PureWindowsPath
import json
import plotly.io as pio
import plotly.graph_objects as go
import subprocess
import random
import statistics
from scipy.stats import mannwhitneyu

tissue_color_dict = {
               'PBMC_pretransplant' : '#67001f',
               'PBMC' : '#67001f',
               'PBL' : '#67001f',               
               'Ileum' : '#4393c3',   
               'Ileum_allograft' : '#4393c3',
               'Colon' : '#2166ac',
               'Colon_allograft' : '#2166ac',
                }

subjects_dict = {   
                2: 'Pt19_R',
                3: 'Pt21_R',
                4: 'Pt23_R',
                8: 'Pt14_R',
                9: 'Pt20_R',
                10: 'Pt17_R',
                17: 'Pt25_R',
                }

subjects_dict_shapes = { 
            'Pt25_R' : 'hourglass', 
            'Pt14_R': 'diamond',
            'Pt17_R': 'cross',
            'Pt19_R': 'circle',
            'Pt20_R': 'square',
            'Pt21_R': 'hexagram',
            'Pt23_R': 'star-triangle-up',
            
            'D168': 'x-thin-open',
            'D181': 'hash-open',
            'D145': 'y-up-open',
            'D182': 'circle-open',
            'D149': 'cross-thin-open',
            'D207': 'asterisk-open',                
            }
        
lp16_subjects_dict = {
                2: 'D145',
                3: 'D149',
                5: 'D168',
                6: 'D181',
                7: 'D182',
                8: 'D207', 
                }

subjects_dict_children = {   
                'Pt20_R' : 1.7,
                'Pt23_R' : 2.2,
                'Pt14_R' : 2.3,
                'Pt21_R' : 2.6,
                'Pt19_R' : 3.3,
                'Pt17_R' : 5.4,     
                'Pt25_R': 9.32,
                }

lp16_tissue_name_convert = {
               'PBL' : 'PBMC',
               'Ileum' : 'Ileum',
               'Colon' : 'Colon',
                }

lp16_tissue_name_revert = {
                'PBMC' : 'PBMC',
                'Ileum_allograft' : 'Ileum',
                'Colon' : 'Colon',
                'Colon_allograft' : 'Colon',
                }

lp16_tissue_name_revert_original_names = {
                'PBMC' : 'PBL',
                'Ileum_allograft' : 'Ileum',
                'Colon' : 'Colon',
                'Colon_allograft' : 'Colon',
                }

lp16_colors = {   
                'D207' : 'black',                    
                'D181' : 'black',
                'D182' : 'black',
                'D149' : 'black',
                'D168' : 'black',
                'D145' : 'black',
            }    

def load_data(patient_name):
    
    # =============================================================================
    # Load data.
    # =============================================================================
    if patient_name in subjects_dict.values():     
        # =============================================================================
        # tissue clones
        # =============================================================================
        path = Path(r"/path_to_sykes_children_data_folder/{}_tissue_clones_pretransplant.pkl".format(patient_name))
        pickle_in = open(path,"rb")
        tissue_clones = pickle.load(pickle_in)
        pickle_in.close()
        
        node_min_clone_dict = pd.read_pickle(r"/path_to_sykes_children_data_folder/{}_clones_per_node_min_using_tree.pkl".format(patient_name))

        # =============================================================================
        # Mutation freq >2% clones
        # =============================================================================
        path = Path(r"/path_to_sykes_children_data_folder/{}_clones_with_mutation_freq_greater_than_2.pkl".format(patient_name))
        pickle_in = open(path,"rb")
        mutation_freq_greater_than_2_clones = pickle.load(pickle_in)
        pickle_in.close()
        
        # =============================================================================
        # Load clumpiness data.
        # =============================================================================
        path = Path(r"/path_to_sykes_children_data_folder/{}_clumpiness_all_tissue_pair_clones_updated_node_min_script.JSON".format(patient_name))
        with open(path) as json_file:
            clumpiness_data = json.load(json_file)   
            
        # =============================================================================
        # Load clones with rejection
        # =============================================================================
        path = Path(r"/path_to_sykes_children_data_folder/{}_clones_containing_rejection.pkl".format(patient_name))
        pickle_in = open(path,"rb")
        rejection_clones = pickle.load(pickle_in)
        pickle_in.close()
                    
            
    elif patient_name in lp16_subjects_dict.values():  
        # =============================================================================
        # tissue clones
        # =============================================================================
        path = Path(r"/path_to_lp16_data_folder/{}_tissue_clones_pretransplant.pkl".format(patient_name))
        pickle_in = open(path,"rb")
        tissue_clones = pickle.load(pickle_in)
        pickle_in.close()
        
        node_min_clone_dict = pd.read_pickle(r"/path_to_lp16_data_folder/{}_clones_per_node_min_using_tree.pkl".format(patient_name))
        
        # =============================================================================
        # Mutation freq >2% clones
        # =============================================================================
        path = Path(r"/path_to_lp16_data_folder/{}_clones_with_mutation_freq_greater_than_2.pkl".format(patient_name))
        pickle_in = open(path,"rb")
        mutation_freq_greater_than_2_clones = pickle.load(pickle_in)
        pickle_in.close()
                
        # =============================================================================
        # Load clumpiness data.
        # =============================================================================
        path = Path(r"/path_to_lp16_data_folder/{}_clumpiness_all_tissue_pair_clones_updated_node_min_script.JSON".format(patient_name))
        with open(path) as json_file:
            clumpiness_data = json.load(json_file)   
            
        # =============================================================================
        # Load clones with rejection
        # =============================================================================
        rejection_clones = set()
                                           
    return tissue_clones, node_min_clone_dict, clumpiness_data, rejection_clones, mutation_freq_greater_than_2_clones

# =============================================================================
# Mann-Whitney U test  
# =============================================================================
def stats_test_Mann_Whitney(dist1, dist2):
    if len(dist1) > 0 and len(dist2) > 0:
        stat, p = mannwhitneyu(dist1, dist2)
        print('Statistics=%.3f, p=%.22f' % (stat, p))
        # interpret
        alpha = 0.05
        if p > alpha:
        	 print('Same distribution (fail to reject H0)')
        else:
        	 print('Different distribution (reject H0)')        
        return f'{p:.21f}'

    else:
        return np.nan
        print("Distribution has no values.")

def get_statistics():
    
    # =============================================================================
    # Set the minimum number of nodes per clone which will be used to filter clones.    
    # =============================================================================
    node_min = 3
        
    # =============================================================================
    # Make multi-categorical x-axis    
    # =============================================================================
    label_order = {
                    'PBMC+Ileum_allograft child' : 1, 
                    'PBMC+Ileum adult' : 2, 
                    'Ileum_allograft+Colon_allograft child' : 3, 
                    'Ileum+Colon adult' : 4, 
                    }
    
    children = list(subjects_dict_children.keys())
    adults = list(lp16_subjects_dict.values())
    
    statistics_values = {}
    
    for x_axis_name, x_axis_value in label_order.items():
        
        statistics_values[x_axis_name] = []
        
        current_tissue = x_axis_name.split(" ")[0]
        first_tissue = current_tissue.split("+")[0]
        second_tissue = current_tissue.split("+")[1]   
        
        if "child" in x_axis_name:
            
            for patient in children:
                
                tissue_clones, node_min_clone_dict, clumpiness_data, rejection_clones, mutation_freq_greater_than_2_clones = load_data(patient)
                
                if first_tissue + "+" + second_tissue not in clumpiness_data.keys() and second_tissue + "+" + first_tissue not in clumpiness_data.keys():
                    continue
                                
                #Find the order of the tissue names
                if first_tissue + "+" + second_tissue in clumpiness_data.keys():
                    current_tissue = first_tissue + "+" + second_tissue
                elif second_tissue + "+" + first_tissue in clumpiness_data.keys():
                    current_tissue = second_tissue + "+" + first_tissue
                                     
                clumpiness_clones = set([int(x) for x in clumpiness_data[current_tissue].keys()])
                node_min_clones = node_min_clone_dict[node_min]
                total_clones_set = clumpiness_clones.intersection(node_min_clones, mutation_freq_greater_than_2_clones)
                
                print("Num clones before removing rejection clones: ", len(total_clones_set))
                # =============================================================================
                # Remove clones with rejection               
                # =============================================================================
                total_clones_set = total_clones_set.difference(rejection_clones)
                print("Num clones after removing rejection clones: ", len(total_clones_set))
                                    
                total_clones = len(total_clones_set)
                if total_clones > 5:
                    first_tissue_key = current_tissue.split("+")[0]
                    second_tissue_key = current_tissue.split("+")[1]   
                    tissue_key = first_tissue_key + "," + second_tissue_key
                    temp = []
                    for clone in total_clones_set:
                        clone_value = clumpiness_data[current_tissue][str(clone)]
                        tissue_value = clone_value[tissue_key]
                        temp.append(tissue_value)
                    median = statistics.median(temp)    
                        
                    statistics_values[x_axis_name].append(median)
                                    
        elif "adult" in x_axis_name:
            
            if first_tissue == 'PBMC':
                first_tissue = 'PBL'
                current_tissue = first_tissue + '+' + second_tissue
            
            for patient in adults:
                             
                tissue_clones, node_min_clone_dict, clumpiness_data, rejection_clones, mutation_freq_greater_than_2_clones = load_data(patient)
                
                if first_tissue + "+" + second_tissue not in clumpiness_data.keys() and second_tissue + "+" + first_tissue not in clumpiness_data.keys():
                    continue
                                
                #Find the order of the tissue names
                if first_tissue + "+" + second_tissue in clumpiness_data.keys():
                    current_tissue = first_tissue + "+" + second_tissue
                elif second_tissue + "+" + first_tissue in clumpiness_data.keys():
                    current_tissue = second_tissue + "+" + first_tissue
                                     
                clumpiness_clones = set([int(x) for x in clumpiness_data[current_tissue].keys()])
                node_min_clones = node_min_clone_dict[node_min]                
                total_clones_set = clumpiness_clones.intersection(node_min_clones, mutation_freq_greater_than_2_clones)
                total_clones = len(total_clones_set)
                
                if total_clones > 5:
                    first_tissue_key = current_tissue.split("+")[0]
                    if first_tissue_key == 'PBMC': ###The outer key uses PBMC but the inner keys use PBL for the PPG adults
                        first_tissue_key = 'PBL'                        
                    second_tissue_key = current_tissue.split("+")[1] 
                    if second_tissue_key == 'PBMC':
                        second_tissue_key = 'PBL'                         
                    tissue_key = first_tissue_key + "," + second_tissue_key
                    temp = []
                    for clone in total_clones_set:
                        clone_value = clumpiness_data[current_tissue][str(clone)]
                        tissue_value = clone_value[tissue_key]
                        temp.append(tissue_value)
                    median = statistics.median(temp)    
                        
                    statistics_values[x_axis_name].append(median)

    categories = list(statistics_values.keys())
    
    tissue_pair_result = []
    p_result = []  
    star_result = []    
    index_tracker = 1
    for first_cat in categories:
        for second_cat in categories[index_tracker:]:
            
            cat_pair = first_cat + " + " + second_cat
            
            #Remove comparisons that are between the same cohort. We will use the Wilcoxon test to test between the same cohort.
            if cat_pair.count('child') == 2 or cat_pair.count('adult') == 2:
                continue   
            
            first_distribution = statistics_values[first_cat]
            second_distribution = statistics_values[second_cat]
            
            p_value = stats_test_Mann_Whitney(first_distribution, second_distribution)
            
            p_result.append(p_value)
            tissue_pair_result.append(cat_pair)
            
            if float(p_value) > .05:
                star_result.append("ns")        
            elif float(p_value) < .05 and float(p_value) > 0.01:
                star_result.append("*")
            elif float(p_value) < .01 and float(p_value) > 0.001:
                star_result.append("**")
            elif float(p_value) < .001:
                star_result.append("***")
            else:
                star_result.append("N/A")
        index_tracker += 1
    
    table_result = pd.DataFrame()
    table_result['Tissue Pair'] = tissue_pair_result 
    table_result['star notation'] = star_result    
    table_result['p value'] = p_result       
    table_result['p value'] = table_result['p value'].apply(lambda x: round(float(x), 5)) 

    # =============================================================================
    # Make table
    # =============================================================================
    df = table_result
    fig = go.Figure(data=[go.Table(
        columnorder = [1,2,3,4],
        columnwidth = [120,50,50,350],    
        header=dict(values=list(df.columns),
                    fill_color='gray',
                    font=dict(color='white', size=20),
                    height=40,                                            
                    align='left'),
        cells=dict(values=[df['Tissue Pair'], df['star notation'], df['p value']],
                   fill_color='whitesmoke',
                    font=dict(size=20),
                    height=40,                                            
                   align='left'))
    ])
    
    pio.write_html(fig, file=r'/path_to_output_save_folder/Supp_Fig_5A_statistics_Mann_Whitney.html')
    print("Finished!")

get_statistics()
