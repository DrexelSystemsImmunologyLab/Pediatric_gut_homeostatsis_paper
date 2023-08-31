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
from itertools import combinations
import random
from math import hypot
import statistics

tissue_color_dict = {
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
                }

subjects_dict_shapes = { 
            'Pt14_R': 'diamond',
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
            }


subjects_dict_colors = {   
                'Pt20_R' : '#DA16FF',
                'Pt23_R' : '#B68100',
                'Pt14_R' : '#750D86',
                'Pt21_R' : '#EB663B',
                'Pt19_R' : '#511CFB',
            }

subjects_dict_children = {   
                'Pt20_R' : 1.7,
                'Pt23_R' : 2.2,
                'Pt14_R' : 2.3,
                'Pt21_R' : 2.6,
                'Pt19_R' : 3.3,
            }

def load_data():
    
    # =============================================================================
    # Load mutation frequency >2% data made from the preprocessing script
    # =============================================================================
    path = Path(r"/path_to_sykes_children_data_folder/Percent_mutation_freq_greater_than_2_by_timepoint_and_node_min.pkl")
    pickle_in = open(path,"rb")
    df = pickle.load(pickle_in)
    pickle_in.close()
    
    df = df.loc[df['patient'].isin(subjects_dict_by_age.keys())]

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
    # Prepare data. We used clones with at least 1 node (i.e., all clones). The data has data for clones of various values of node_min, so we must filter those out.
    # =============================================================================        
    df = df_original.loc[df_original['node_min'] == 1] 
    
    # =============================================================================
    # Define axis labels
    # =============================================================================        
    fig = go.Figure()
    pod_brackets = ['0', '1-90', '91-365', '>365']
    pod_brackets_text = ['Pre', 'Early', 'Mid', 'Late']
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
        
        for pod_bracket, x_value in zip(pod_brackets, x_axis_positions):
            if "-" in pod_bracket:
                pod_min = int(pod_bracket.split("-")[0])
                pod_max = int(pod_bracket.split("-")[1])
            elif pod_bracket == '0':
                pod_min = pod_max = 0
            else:
                pod_min = 366
                pod_max = 9999999

            pod_df = tissue_df.loc[(tissue_df['pod'] >= pod_min) & (tissue_df['pod'] <= pod_max)]
            pod_df = pod_df.loc[pod_df['total_clones'] > 5]
            
            #Convert to median
            pod_df = pod_df.groupby(['patient'])['percent_mutation_freq_greater_than_2'].apply(lambda x: statistics.median(x)).reset_index()
            
            y_values = pod_df['percent_mutation_freq_greater_than_2'].tolist()            
            x_values = [x_value for x in y_values] 
            marker_color = [tissue_color_dict[tissue] for x in y_values]                
            marker_symbol = [subjects_dict_shapes[current_patient] for current_patient in pod_df['patient']]     
            hovertext = ['Patient: ' + patient_name + '<br>Tissue: ' + tissue for patient_name in pod_df['patient']]
            
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
                marker_size=16,
                marker_line_color=marker_line_color,
                marker_line_width=marker_line_width,
                hovertext=hovertext,
                opacity=0.85,
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
                         
    fig.update_yaxes(range=[-0.02, 1.06 ], tickfont_size=26) 
    fig.update_xaxes(tickfont_size=26, ticktext=pod_brackets_text, tickvals=x_axis_positions) 
        
    fig.update_layout(title_text="Fraction clones with a mutation frequency >2% per time point.",
                      yaxis_title="Fraction clones",
                      titlefont_size=12,
                      font=dict(
                            family="Calibri",
                            size=26,
                            color="black"
                           ),
                      autosize=False,
                      width=700,
                      height=375,
                      plot_bgcolor='whitesmoke',    
                      )
    
    pio.write_html(fig, file=r'/path_to_output_save_folder/Fig_4B.html') #produces an interactive .html graph. You may want to comment out lines "autosize=False, width=some_number, height=some_number" in fig.update_layout to allow the size to be automatically made for this file.
    fig.write_image('/path_to_output_save_folder/Fig_4B.png', scale=8) #produces a high res .png image
    print("Finished!")
    
make_graph(df_original)
