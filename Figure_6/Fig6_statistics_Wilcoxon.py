# -*- coding: utf-8 -*-
"""
Created on Thu Jul  2 13:24:30 2020

@author: thsia
"""

import json
import pandas as pd
import numpy as np
from pathlib import Path, PureWindowsPath
import os
import re
import plotly.io as pio
import plotly.graph_objects as go
import statistics
import pickle
from scipy.stats import wilcoxon

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
        # Load clumpiness data.
        # =============================================================================
        path = Path(r"/path_to_lp16_data_folder/{}_clumpiness_all_tissue_pair_clones_updated_node_min_script.JSON".format(patient_name))
        with open(path) as json_file:
            clumpiness_data = json.load(json_file)   
            
        # =============================================================================
        # Load clones with rejection
        # =============================================================================
        rejection_clones = set()

    return tissue_clones, node_min_clone_dict, clumpiness_data, rejection_clones

def get_statistics():
    
    # =============================================================================
    # Set the minimum number of nodes per clone which will be used to filter clones.    
    # =============================================================================
    node_min = 3
    
    # =============================================================================
    # Make axis labels  
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
    statistics_values_patients = {}
    
    for x_axis_name, x_axis_value in label_order.items():
        
        statistics_values[x_axis_name] = []
        statistics_values_patients[x_axis_name] = []
                   
        current_tissue = x_axis_name.split(" ")[0]
        first_tissue = current_tissue.split("+")[0]
        second_tissue = current_tissue.split("+")[1]   
        
        if "child" in x_axis_name:
            
            for patient in children:
                
                tissue_clones, node_min_clone_dict, clumpiness_data, rejection_clones = load_data(patient)
                
                if first_tissue + "+" + second_tissue not in clumpiness_data.keys() and second_tissue + "+" + first_tissue not in clumpiness_data.keys():
                    continue
                                
                #Find the order of the tissue names
                if first_tissue + "+" + second_tissue in clumpiness_data.keys():
                    current_tissue = first_tissue + "+" + second_tissue
                elif second_tissue + "+" + first_tissue in clumpiness_data.keys():
                    current_tissue = second_tissue + "+" + first_tissue
                                     
                clumpiness_clones = set([int(x) for x in clumpiness_data[current_tissue].keys()])
                node_min_clones = node_min_clone_dict[node_min]
                total_clones_set = clumpiness_clones.intersection(node_min_clones)
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
                    statistics_values_patients[x_axis_name].append(patient)
                    
        elif "adult" in x_axis_name:
            
            if first_tissue == 'PBMC':
                first_tissue = 'PBL'
                current_tissue = first_tissue + '+' + second_tissue
            
            for patient in adults:
                             
                tissue_clones, node_min_clone_dict, clumpiness_data, rejection_clones = load_data(patient)
                
                if first_tissue + "+" + second_tissue not in clumpiness_data.keys() and second_tissue + "+" + first_tissue not in clumpiness_data.keys():
                    continue
                                
                #Find the order of the tissue names
                if first_tissue + "+" + second_tissue in clumpiness_data.keys():
                    current_tissue = first_tissue + "+" + second_tissue
                elif second_tissue + "+" + first_tissue in clumpiness_data.keys():
                    current_tissue = second_tissue + "+" + first_tissue
                                     
                clumpiness_clones = set([int(x) for x in clumpiness_data[current_tissue].keys()])
                node_min_clones = node_min_clone_dict[node_min]
                total_clones_set = clumpiness_clones.intersection(node_min_clones)
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
                    statistics_values_patients[x_axis_name].append(patient)

    categories = list(statistics_values.keys())
    
    print("Now calculating Wilcoxon...")
    
    comparison_list = []
    patients_compared_list = []
    p_values_dict = {'two-sided' : [], 'greater' : [], 'less' : []}
    star_notation_dict = {'two-sided' : [], 'greater' : [], 'less' : []}
    index_tracker = 1
    for first_cat in categories:
        first_cat_patients = statistics_values_patients[first_cat]
        for second_cat in categories[index_tracker:]:
            second_cat_patients = statistics_values_patients[second_cat]
            
            overlapping_patients = set(first_cat_patients).intersection(set(second_cat_patients))
            if len(overlapping_patients) < 2:
                continue
                   
            cat_pair = first_cat + " + " + second_cat
            
            first_distribution = [x for x, y in zip(statistics_values[first_cat], first_cat_patients) if y in overlapping_patients]
            second_distribution = [x for x, y in zip(statistics_values[second_cat], second_cat_patients) if y in overlapping_patients]
  
            for alternative in ['two-sided', 'greater', 'less']:
                
                w, p_value = wilcoxon(x=first_distribution, y=second_distribution, alternative=alternative)
                
                if float(p_value) > .05:
                    star_notation_dict[alternative].append("ns")        
                elif float(p_value) < .05 and float(p_value) > 0.01:
                    star_notation_dict[alternative].append("*")
                elif float(p_value) < .01 and float(p_value) > 0.001:
                    star_notation_dict[alternative].append("**")
                elif float(p_value) < .001:
                    star_notation_dict[alternative].append("***")
                else:
                    star_notation_dict[alternative].append("N/A")      
                p_values_dict[alternative].append(p_value)            
    
            patients_compared_list.append(list(overlapping_patients))
            comparison_list.append(cat_pair)

        index_tracker += 1
    
    table_result = pd.DataFrame()
    table_result['Comparison'] = comparison_list 
    table_result['two-sided star notation'] = star_notation_dict['two-sided']    
    table_result['two-sided p value'] = p_values_dict['two-sided']       
    table_result['two-sided p value'] = table_result['two-sided p value'].apply(lambda x: round(float(x), 5)) 
    table_result['greater star notation'] = star_notation_dict['greater']    
    table_result['greater p value'] = p_values_dict['greater']           
    table_result['greater p value'] = table_result['greater p value'].apply(lambda x: round(float(x), 5)) 
    table_result['less star notation'] = star_notation_dict['less']        
    table_result['less p value'] = p_values_dict['less']           
    table_result['less p value'] = table_result['less p value'].apply(lambda x: round(float(x), 5))         
    table_result['patients_compared'] = patients_compared_list 
       
    # =============================================================================
    # Make table
    # =============================================================================
    df = table_result
    fig = go.Figure(data=[go.Table(
        columnwidth = [30 if x not in ['Comparison', 'patients_compared'] else 200 for x in df.columns],    
        header=dict(values=list(df.columns),
                    fill_color='gray',
                    font=dict(color='white', size=20),
                    height=40,                                            
                    align='left'),
        cells=dict(values=[df[x] for x in df.columns],
                   fill_color='whitesmoke',
                    font=dict(size=20),
                    height=40,                                            
                   align='left'))
    ])
    
    pio.write_html(fig, file=r'/path_to_output_save_folder/Fig_6_statistics_Wilcoxon.html') #Creates an .html file containing results of statistics
    print("Finished!")
    
get_statistics()    
    