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
import plotly.io as pio
import plotly.graph_objects as go
import json
import statistics
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
            }
        
subjects_dict_by_age = {   
                'Pt20_R' : 1.7,
                'Pt23_R' : 2.2,
                'Pt14_R' : 2.3,
                'Pt21_R' : 2.6,
                'Pt19_R' : 3.3,
                'Pt17_R' : 5.4,     
                'Pt25_R': 9.32,
            }


def load_data():
    
    df = pd.read_pickle(r"/path_to_sykes_children_data_folder/Repertoire_diversity_by_POD_BRACKETS_and_tissue_and_diversity_order.pkl")
   
    return df

df_original = load_data()
  
def get_statistics(df_original):    
    
    # =============================================================================
    # Prepare data
    # =============================================================================        
    df = df_original.loc[df_original['node_min'] == 1]
    df = df.loc[df['total_clones'] > 5]
    current_patients = list(subjects_dict_by_age.keys())
    #Remove patients who do not have at least 1 data point in Early, Mid, and Late
    df = df.loc[df['patient'].isin(current_patients)]
      
    # =============================================================================
    # Make graphs
    # =============================================================================
    comparison_list = []
    diversity_order_result = []
    patients_compared_list = []
    p_values_dict = {'two-sided' : [], 'greater' : [], 'less' : []}
    star_notation_dict = {'two-sided' : [], 'greater' : [], 'less' : []} 
    pod_brackets = ['0', '1-90', '91-365', '>365']
    for tissue in ['PBMC', 'Ileum_allograft']:
        
        if tissue == 'PBMC':
            str_contains = 'PBMC'
        elif tissue == 'Ileum_allograft':
            str_contains = 'Ileum'
        tissue_df = df.loc[df['tissue'].str.contains(str_contains)]
        
        for pod_bracket in pod_brackets:
            
            current_index = pod_brackets.index(pod_bracket)
            pod_df = tissue_df.loc[tissue_df['pod'] == pod_bracket]
            pod_df = pod_df.loc[pod_df['total_clones'] > 5]
            pod_df_patients = pod_df['patient'].unique().tolist()           
            
            for remaining_brackets in pod_brackets[current_index + 1:]:
              
                remaining_brackets_pod_df = tissue_df.loc[tissue_df['pod'] == remaining_brackets]
                remaining_brackets_pod_df = remaining_brackets_pod_df.loc[remaining_brackets_pod_df['total_clones'] > 5]
                remaining_brackets_pod_df_patients = remaining_brackets_pod_df['patient'].unique().tolist()     
                
                patients_in_both = [x for x in remaining_brackets_pod_df_patients if x in pod_df_patients]
                name_of_comparison = "{}: {} | {}".format(tissue, pod_bracket, remaining_brackets)
                print("{}; {} patients in this comparison.".format(name_of_comparison, len(patients_in_both)))
                
                #Filter for current patients
                upper_loop_df = pod_df.loc[pod_df['patient'].isin(patients_in_both)]
                lower_loop_df = remaining_brackets_pod_df.loc[remaining_brackets_pod_df['patient'].isin(patients_in_both)]
                  
                # =============================================================================
                # UPPER loop y_values
                # Get richness                 
                # =============================================================================
                richness_df_upper = upper_loop_df.loc[upper_loop_df['diversity_order'] == 0]
                richness_upper = richness_df_upper['diversity'].tolist()
                
                order_df_upper = upper_loop_df.loc[upper_loop_df['diversity_order'] == 1]
                y_values_upper_loop = order_df_upper['diversity'].tolist()

                # =============================================================================
                # Normalize by richess
                # =============================================================================
                y_values_upper_loop = [y/richness_val for y, richness_val in zip(y_values_upper_loop, richness_upper)]
                                        
                # =============================================================================
                # LOWER loop y_values
                # Get richness                 
                # =============================================================================
                richness_df_lower = lower_loop_df.loc[lower_loop_df['diversity_order'] == 0]
                richness_lower = richness_df_lower['diversity'].tolist()
                
                order_df_lower = lower_loop_df.loc[lower_loop_df['diversity_order'] == 1]
                y_values_lower_loop = order_df_lower['diversity'].tolist()

                # =============================================================================
                # Normalize by richess
                # =============================================================================
                y_values_lower_loop = [y/richness_val for y, richness_val in zip(y_values_lower_loop, richness_lower)]

                for alternative in ['two-sided', 'greater', 'less']:
                    
                    w, p_value = wilcoxon(x=y_values_upper_loop, y=y_values_lower_loop, alternative=alternative)
                    
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
                
                comparison_list.append(name_of_comparison)
                patients_compared_list.append(patients_in_both)    
                diversity_order_result.append(1)

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
    table_result['diversity_order'] = diversity_order_result    
             
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
    
    pio.write_html(fig, file=r'/path_to_output_save_folder/Fig_5A_statistics_Wilcoxon.html') #Creates an .html file containing results of statistics
    print("Finished!")

get_statistics(df_original)
