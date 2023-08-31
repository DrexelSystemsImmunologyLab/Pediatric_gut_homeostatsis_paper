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
import json
import plotly.io as pio
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import statistics
from scipy.stats import binomtest

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
                9: 'Pt20_R',
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
                'Pt21_R' : 2.6,
                'Pt19_R' : 3.3,
            }

subjects_dict_children = {   
                'Pt20_R' : 1.7,
                'Pt23_R' : 2.2,
                'Pt21_R' : 2.6,
                'Pt19_R' : 3.3,
            }

def load_data():
    
    df = pd.read_pickle(r"/path_to_sykes_children_data_folder/Clone_size_by_POD_bracket_INSTANCES_clones_spanning_two_consecutive_brackets.pkl")
   
    return df

df_original = load_data()
  
def make_graph(df_original):     

    
    # =============================================================================
    # Prepare data
    # =============================================================================        
    df = df_original
    patients = df['patient'].unique()
    
    #Change Ileum to Ileum_allograft so it is included in the analysis
    df['tissue'] = df['tissue'].replace('Ileum', 'Ileum_allograft') 

    # =============================================================================
    # Define graph layout
    # ============================================================================= 
    POD_brackets = ['0 -> 1-90', '1-90 -> 91-365', '91-365 -> >365']
    pod_bracket_sets = [('0', '1-90'), ('1-90', '91-365'), ('91-365', '>365')]
    POD_brackets_text_names = ['P -> E', 'E -> M', 'M -> L']
    titles = ['Blood', 'Ileum allograft']
    num_cols = len(titles)
    current_patients = [x for x in subjects_dict_by_age.keys() if x in patients]
    num_rows = len(current_patients)
    subplot_titles = []
    for title in titles:
        subplot_titles.append("<b>" + title + "</b>")
    fig = make_subplots(rows=num_rows, cols=num_cols, subplot_titles=subplot_titles, 
                        vertical_spacing=0.035,
                        horizontal_spacing=0.015,
                        shared_yaxes=True, 
                        shared_xaxes=True
                        )
    
    ###Add remaining empty titles
    Num_remaining = num_cols * num_rows - num_cols
    for remaining_titles in range(Num_remaining):
        subplot_titles.append("")
    
    for i in fig['layout']['annotations']:
        i['font'] = dict(size=18, family='Calibri', color='black')
        
    col_number = 1
    row_number = 1   
    # =============================================================================
    # Make graphs
    # =============================================================================
    categories = POD_brackets
    x_numerical_vals = [1, 2, 3]
    value_types = ["+", "-"]
    colors = ['green', 'red']        
    for patient_name in current_patients:
     
        patient_df = df.loc[df['patient'] == patient_name]
        
        #get max y value in patient to set range of graph
        max_patient_y_val = 0
        for tissue in ['PBMC', 'Ileum_allograft']:
            
            if tissue == 'PBMC':
                str_contains = 'PBMC'
            elif tissue == 'Ileum_allograft':
                str_contains = 'Ileum'
            tissue_df = patient_df.loc[patient_df['tissue'].str.contains(str_contains)]
                         
            for category, bracket_set in zip(categories, pod_bracket_sets):
                bracket_df = tissue_df.loc[tissue_df['pod_bracket_sets'] == bracket_set]
                max_val = tissue_df['size_change_sign'].value_counts()
                if len(max_val) == 0:
                    continue
                max_val = max(max_val["+"], max_val["-"])
                if max_val > max_patient_y_val:
                    max_patient_y_val = max_val


        for tissue in ['PBMC', 'Ileum_allograft']:
            
            if tissue == 'PBMC':
                str_contains = 'PBMC'
            elif tissue == 'Ileum_allograft':
                str_contains = 'Ileum'
            tissue_df = patient_df.loc[patient_df['tissue'].str.contains(str_contains)]
                                     
            max_y = 0
            for value, color in zip(value_types, colors):
                y_values = []
                x_values = []
                hovertext = []                    
                for category, bracket_set, x_val in zip(categories, pod_bracket_sets, x_numerical_vals):
                    bracket_df = tissue_df.loc[tissue_df['pod_bracket_sets'] == bracket_set]
                    counts = bracket_df['size_change_sign'].value_counts()
                    if len(counts) == 0 or "+" not in counts.keys() or "-" not in counts.keys():
                        continue
                    if value in counts.keys():
                        count = counts[value]
                        if count > max_y:
                            max_y = count
                    else:
                        count = 0
                    y_values.append(count)
                    x_values.append(x_val)
                    
                    # =============================================================================
                    # Add p-value using the binomial test               
                    # =============================================================================
                    if value != "tie":
                        num_pos = counts["+"]
                        num_neg = counts["-"]
                        total = num_pos + num_neg
                        binomtest_result = binomtest(num_pos, n=total, p=0.5)
                        p_value = binomtest_result.pvalue
                        
                    hovertext.append("Tissue: {}<br>Category: {}<br>Count: {}<br>p-value: {}".format(tissue, value, count, p_value))
                    
                fig.add_trace(go.Bar(
                    x=x_values,
                    y=y_values,
                    name=value,
                    marker_color=color,
                    hovertext=hovertext,
                    opacity=0.8,
                    showlegend=False
                ),row=row_number, col=col_number)
            
            # =============================================================================
            # Add p-value to graph                     
            # =============================================================================
            for category, bracket_set, x_val in zip(categories, pod_bracket_sets, x_numerical_vals):
                bracket_df = tissue_df.loc[tissue_df['pod_bracket_sets'] == bracket_set]
                counts = bracket_df['size_change_sign'].value_counts()
                if len(counts) == 0 or "+" not in counts.keys() or "-" not in counts.keys():
                    continue
                num_pos = counts["+"]
                num_neg = counts["-"]
                total = num_pos + num_neg
                binomtest_result = binomtest(num_pos, n=total, p=0.5)
                p_value = binomtest_result.pvalue
                        
                if p_value >= 0.05:
                    continue
                elif p_value < 0.05:
                    star_value = '*'
                    
                fig.add_trace(go.Scatter(
                    x=[x_val-.25, x_val, x_val+0.25],
                    y=[max_patient_y_val + max_patient_y_val*0.15, max_patient_y_val + max_patient_y_val*0.15, max_patient_y_val + max_patient_y_val*0.15],
                    mode="lines+text",
                    name="p-val for " + category,
                    text=["", star_value, ""],
                    textfont_size=16,
                    line_color='black',
                    textposition="top center",
                    showlegend=False,
                ),row=row_number, col=col_number)
                
            fig.update_yaxes(range=[0, max_patient_y_val + max_patient_y_val*0.5], row=row_number, col=col_number)
                                        
            if col_number == 1:
                fig.update_yaxes(title=patient_name, row=row_number, col=col_number, titlefont_size=18)
            col_number += 1

        col_number = 1
        row_number += 1    
    # # =============================================================================
    # # Add legend               
    # # =============================================================================
    fig.add_trace(go.Scatter(
        x=[None],
        y=[None],
        name='+',
        marker_color='green',
        mode='markers',
        ), 
        row=1, col=1)   
    
    fig.add_trace(go.Scatter(
        x=[None],
        y=[None],
        name='-',
        marker_color='red',
        mode='markers',
        ), 
        row=1, col=1)    
    
    star_values = ['*: p < 0.05']
    for star_val in star_values:
        fig.add_trace(go.Scatter(
            x=[None],
            y=[None],
            name=star_val,
            marker_color='black',
            mode='markers',
            ), 
            row=1, col=1) 
               
    fig.update_yaxes(tickfont_size=15, nticks=4, tickfont_family='Calibri', tickfont_color='black')
    fig.update_xaxes(tickfont_size=15, tickfont_family='Calibri', tickfont_color='black', tickvals=x_numerical_vals, ticktext=POD_brackets_text_names)

    fig.update_layout(title_text="Change in clone size (instances) by POD bracket.<br>Filtered for clones spanning two consecutive POD brackets.",
                      titlefont_size=12,
                      font=dict(
                            family="Calibri",
                            color="black"
                           ),
                       width=700,
                       height=500,                          
                      plot_bgcolor='whitesmoke',
                      barmode='group',
                      # margin_t=150,
    
    )
    
    pio.write_html(fig, file=r'/path_to_output_save_folder/Fig_5B.html') #produces an interactive .html graph. You may want to comment out lines "autosize=False, width=some_number, height=some_number" in fig.update_layout to allow the size to be automatically made for this file.
    fig.write_image('/path_to_output_save_folder/Fig_5B.png', scale=8) #produces a high res .png image
    print("Finished!")
    
make_graph(df_original)
