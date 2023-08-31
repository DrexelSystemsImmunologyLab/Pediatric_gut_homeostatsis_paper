# -*- coding: utf-8 -*-
"""
Created on Thu Dec  5 14:20:09 2019

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
import random
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
        # Mutation freq data
        # =============================================================================
        path = Path(r"/path_to_sykes_children_data_folder/{}_mutation_freq_per_clone_per_tissue_pretransplant_version.pkl".format(patient_name))
        pickle_in = open(path,"rb")
        data = pickle.load(pickle_in)
        pickle_in.close()
                    
    elif patient_name in lp16_subjects_dict.values():  
        # =============================================================================
        # Mutation freq data
        # =============================================================================
        path = Path(r"/path_to_lp16_data_folder/{}_mutation_freq_per_clone_per_tissue_pretransplant_version.pkl".format(patient_name))
        pickle_in = open(path,"rb")
        data = pickle.load(pickle_in)
        pickle_in.close()
              
    return data

# =============================================================================
# Mann-Whitney
# =============================================================================
def stats_test_Mann_Whitney(dist1, dist2):
    if len(dist1) > 0 and len(dist2) > 0:
        stat, p = mannwhitneyu(dist1, dist2)
        print('Statistics=%.3f, p=%.22f' % (stat, p))
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
    # Define x_axis positions and categories.    
    # =============================================================================
    label_order = {
                    "PBMC_pretransplant child" : 1,
                    'PBMC child' : 2, 
                    'PBMC adult' : 3, 
                    'Ileum_allograft child' : 4, 
                    'Ileum adult' : 5, 
                    'Colon child' : 6, 
                    'Colon_allograft child' : 7,
                    'Colon adult' : 8,
                   }
    
    children = list(subjects_dict_children.keys())
    adults = list(lp16_subjects_dict.values())
        
    statistics_values = {}
    
    for x_axis_name, x_axis_value in label_order.items():
        
        statistics_values[x_axis_name] = []
        
        x_values = []
        y_values = []           
        current_tissue = x_axis_name.split(" ")[0]
        
        if "child" in x_axis_name:
            
            for patient in children:
                
                mutation_freq_data = load_data(patient)
                
                if current_tissue not in mutation_freq_data.keys():
                    continue
                                
                mutation_freq_data = mutation_freq_data[current_tissue]
                current_clones = set(mutation_freq_data.keys())
                total_clones = len(current_clones)
                    
                if total_clones > 5:
                    mutation_freq_greater_than_2 = [x for x, y in mutation_freq_data.items() if y > 2 and x in current_clones]
                    percent_greater_than_2 = len(mutation_freq_greater_than_2) / total_clones
    
                    x_values.append(x_axis_value)
                    y_values.append(percent_greater_than_2)
                    statistics_values[x_axis_name].append(percent_greater_than_2)
    
        elif "adult" in x_axis_name:
            
            if current_tissue == 'PBMC':
                current_tissue = 'PBL'
            if current_tissue == 'Spleen':
                current_tissue = 'SPL'
                                
            for patient in adults:
                
                mutation_freq_data = load_data(patient)
                
                if current_tissue not in mutation_freq_data.keys():
                    continue
                                
                mutation_freq_data = mutation_freq_data[current_tissue]
                current_clones = set(mutation_freq_data.keys())
                total_clones = len(current_clones)
                    
                if total_clones > 5:
                    mutation_freq_greater_than_2 = [x for x, y in mutation_freq_data.items() if y > 2 and x in current_clones]
                    percent_greater_than_2 = len(mutation_freq_greater_than_2) / total_clones
    
                    x_values.append(x_axis_value)
                    y_values.append(percent_greater_than_2)
                    statistics_values[x_axis_name].append(percent_greater_than_2)
    
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
        
        pio.write_html(fig, file=r'/path_to_output_save_folder/Fig_4A_Statistics_table_Mann_Whitney.html') #Creates an .html file containing results of statistics
        print("Finished!")
        
get_statistics()
