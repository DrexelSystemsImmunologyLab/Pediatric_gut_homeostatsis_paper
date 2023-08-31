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


def make_fig():
  
    # =============================================================================
    # Manually add jitter
    # =============================================================================
    def add_jitter(x_values, y_values):
        
        def distance(p1,p2):
            """Euclidean distance between two points."""
            x1,y1 = p1
            x2,y2 = p2
            return hypot(x2 - x1, y2 - y1)
 
        ###Try to minimize overlap
        result = {}
        num_tries = 0
        current_best = 0
        while num_tries < 100000:
            temp = [x + random.uniform(-0.15, 0.15 ) for x in x_values]
            list_of_coords = [(x, y) for x, y in zip(temp, y_values)]
            min_distance = min([distance(*combo) for combo in combinations(list_of_coords,2)])
            num_tries += 1
            if min_distance > current_best:
                current_best = min_distance
                result = temp
        return result
    
    fig = go.Figure()
    
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
    
    for x_axis_name, x_axis_value in label_order.items():
        
        x_values = []
        y_values = []
        marker_symbol = []
        marker_color = []
        marker_line_color = []
        marker_line_width = []
        hovertext = []            
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
                    marker_symbol.append(subjects_dict_shapes[patient])
                    marker_color.append(tissue_color_dict[current_tissue])
                    if 'allograft' in current_tissue:
                        marker_line_color.append('orchid') 
                        marker_line_width.append(2)                      
                    else:
                        marker_line_color.append('black') 
                        marker_line_width.append(1)                      
                    hovertext.append("Patient: {}<br>Clones: {}<br>Total clones: {}<br>Tissue: {}".format(patient, len(mutation_freq_greater_than_2), total_clones, current_tissue))
        
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
                    marker_symbol.append(subjects_dict_shapes[patient])
                    marker_color.append(tissue_color_dict[current_tissue])                   
                    marker_line_color.append(tissue_color_dict[current_tissue]) 
                    marker_line_width.append(2)                      
                    hovertext.append("Patient: {}<br>Clones: {}<br>Total clones: {}<br>Tissue: {}".format(patient, len(mutation_freq_greater_than_2), total_clones, current_tissue))
        
        if len(x_values) == 1:    
            x_values = x_values
        elif len(x_values) > 1:
            x_values = add_jitter(x_values, y_values)
        else:
            continue              
      
        fig.add_trace(go.Scatter(
            x=x_values, 
            y=y_values,
            mode='markers',
            marker_size=16,        
            marker_color=marker_color,
            marker_line_color=marker_line_color,
            marker_line_width=marker_line_width,                
            marker_symbol=marker_symbol,
            hovertext=hovertext,
            opacity=0.85,
            showlegend=False,
        ))

    fig.update_yaxes(range=[-0.04, 1.05], tickfont_size=28)
    fig.update_xaxes(tickangle=25,
                     tickvals=list(label_order.values()),
                     ticktext=list(label_order.keys()),
                     tickfont_size=28,                     
                     )          
     
    fig.update_layout(title_text="Fraction of clones with a mutation frequency >2%.",
                      yaxis_title="Fraction clones",
                      titlefont_size=12,
                      font=dict(
                            family="Calibri",
                            size=28,
                            color="black"
                           ),
                       autosize=False,
                       width=800,
                       height=500,
                      plot_bgcolor='whitesmoke',
    
    )
    
    pio.write_html(fig, file=r'/path_to_output_save_folder/Fig_4A.html') #produces an interactive .html graph. You may want to comment out lines "autosize=False, width=some_number, height=some_number" in fig.update_layout to allow the size to be automatically made for this file.
    fig.write_image('/path_to_output_save_folder/Fig_4A.png', scale=8) #produces a high res .png image
    print("Finished!")
        
make_fig()  

