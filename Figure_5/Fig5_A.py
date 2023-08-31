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
import statistics


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
        temp = [x + random.uniform(-0.05, 0.05 ) for x in x_values]
        list_of_coords = [(x, y) for x, y in zip(temp, y_values)]
        min_distance = min([distance(*combo) for combo in combinations(list_of_coords,2)])
        num_tries += 1
        if min_distance > current_best:
            current_best = min_distance
            result = temp
    return result
       
def make_graph(df_original):     
   
    # =============================================================================
    # Filter data
    # =============================================================================        
    df = df_original.loc[df_original['node_min'] == 1]
    df = df.loc[df['total_clones'] > 5]
    
    #Remove patients who do not have at least 1 data point in Early, Mid, and Late
    current_patients = list(subjects_dict_by_age.keys())
    df = df.loc[df['patient'].isin(current_patients)]
    
    # =============================================================================
    # Define graph layout
    # =============================================================================        
    fig = go.Figure()
    pod_brackets = ['0', '1-90', '91-365', '>365']
    pod_brackets_names = ['Pre', 'Early', 'Mid', 'Late']
    x_axis_positions = [0, 1, 2, 3]

    # =============================================================================
    # Make graphs
    # =============================================================================
    for tissue in ['PBMC', 'Ileum_allograft']:
        
        if tissue == 'PBMC':
            str_contains = 'PBMC'
        elif tissue == 'Ileum_allograft':
            str_contains = 'Ileum'
        tissue_df = df.loc[df['tissue'].str.contains(str_contains)]
        
        for pod_bracket, x_value, pod_brackets_name in zip(pod_brackets, x_axis_positions, pod_brackets_names):

            pod_df = tissue_df.loc[tissue_df['pod'] == pod_bracket]
                            
            # =============================================================================
            # Filter for richness and diversity order 1                
            # =============================================================================
            richness_df = pod_df.loc[pod_df['diversity_order'] == 0]
            richness = richness_df['diversity'].tolist()
            
            order_df = pod_df.loc[pod_df['diversity_order'] == 1] #This graph uses diversity order 1
            y_values = order_df['diversity'].tolist()

            # =============================================================================
            # Normalize by richess
            # =============================================================================
            y_values = [y/richness_val for y, richness_val in zip(y_values, richness)]
            
            x_values = [x_value for x in y_values] 
            marker_color = [tissue_color_dict[tissue] for x in y_values]                
            marker_symbol = [subjects_dict_shapes[current_patient] for current_patient in order_df['patient']]     
            hovertext = ['Patient: ' + patient_name + '<br>Tissue: ' + tissue + '<br>Total clones: ' + str(total_clones) for patient_name, tissue, total_clones in zip(order_df['patient'], order_df['tissue'], order_df['total_clones'])]
                
            if pod_bracket != '0' and tissue == 'Ileum_allograft':
                marker_line_color = ['orchid' for x in y_values]
                marker_line_width = [1 for x in y_values]
            else:
                marker_line_color = ['black' for x in y_values]
                marker_line_width = [0.5 for x in y_values]
                
            #Shift PBMC markers slightly to the left and ileum markers to the right
            shift_distance = 0.2
            if tissue == 'PBMC':
                x_values = [x - shift_distance for x in x_values]
                median_x = x_value - shift_distance
            elif tissue == 'Ileum_allograft':
                x_values = [x + shift_distance for x in x_values]
                median_x = x_value + shift_distance

            x_values = add_jitter(x_values, y_values)
            fig.add_trace(go.Scatter(
                x=x_values,
                y=y_values,
                mode='markers',
                marker_color=marker_color,
                marker_symbol=marker_symbol,
                marker_line_color=marker_line_color,
                marker_line_width=marker_line_width,
                marker_size=16,
                hovertext=hovertext,
                opacity=0.8,
                showlegend=False
            ))
            
            #Add median
            median = statistics.median(y_values)               
            fig.add_trace(go.Scatter(
                x=[median_x - 0.17, median_x + 0.17],
                y=[median, median],
                mode='lines',
                line_color='black',
                line_width=2.5,
                showlegend=False
            ))                
                         
    fig.update_yaxes(range=[-0.05, 1.05 ], tickfont_size=28, title='Evenness', title_font_size=28) 
    fig.update_xaxes(tickfont_size=28, ticktext=pod_brackets_names, tickvals=x_axis_positions, title='POD bracket', title_font_size=28) 
                                       
    fig.update_layout(title_text="Repertoire diversity by clone size (instances).<br>Diveristy order 1/richness",
                      titlefont_size=12,
                      font=dict(
                            family="Calibri",
                            color="black"
                           ),
                       width=700,
                       height=375,
                      plot_bgcolor='whitesmoke',
                      # margin_t=150,
                      )
    
    pio.write_html(fig, file=r'/path_to_output_save_folder/Graphs/Fig_5A.html') #produces an interactive .html graph. You may want to comment out lines "autosize=False, width=some_number, height=some_number" in fig.update_layout to allow the size to be automatically made for this file.
    fig.write_image('/path_to_output_save_folder/Fig_5A.png', scale=8) #produces a high res .png image
    print("Finished!")
    
make_graph(df_original)
