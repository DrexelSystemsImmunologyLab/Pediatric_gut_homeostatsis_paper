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
from sklearn.linear_model import LinearRegression
import subprocess
import statistics
from scipy.optimize import curve_fit

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
    # Sykes data
    # =============================================================================
    path = Path(r"/path_to_sykes_children_data_folder/Percent_mutation_freq_greater_than_2_by_timepoint_and_node_min.pkl")
    pickle_in = open(path,"rb")
    df = pickle.load(pickle_in)
    pickle_in.close()
    
    df = df.loc[df['patient'].isin(subjects_dict_by_age.keys())]
    
    # =============================================================================
    # Adult data
    # =============================================================================
    path = Path(r"/path_to_lp16_data_folder/Percent_mutation_freq_greater_than_2_by_timepoint_and_node_min_POD_liz_brackets.pkl")
    pickle_in = open(path,"rb")
    adult_df_original = pickle.load(pickle_in)   
    pickle_in.close() 

    return df, adult_df_original

df_original, adult_df_original = load_data()
  
def make_graph(df_original, adult_df_original):     
    
    # =============================================================================
    # Prepare data
    # =============================================================================  
    node_min = 1      
    df = df_original.loc[df_original['node_min'] == node_min]
    adult_df = adult_df_original.loc[adult_df_original['node_min'] == node_min]
    
    #Remove patients who do not have both PBMC and Ileum_allograft
    patients_to_be_removed = []
    for patient_name in subjects_dict_by_age.keys():
        
        patient_df = df.loc[df['patient'] == patient_name]
        patient_tissue_list = patient_df['tissue'].unique().tolist()

        skip_patients_that_do_not_have_PBMC_and_ileum_allograft = set(patient_tissue_list).intersection(set(['PBMC', 'Ileum_allograft']))
        if len(skip_patients_that_do_not_have_PBMC_and_ileum_allograft) < 2:
            patients_to_be_removed.append(patient_name)
    df = df.loc[~df['patient'].isin(patients_to_be_removed)]
        
    # =============================================================================
    # Define graph layout
    # =============================================================================        
    titles = ['Blood', 'Ileum allograft']
    num_cols = len(titles)
    num_rows = 1
    subplot_titles = []
    for title in titles:
        subplot_titles.append("<b>" + title + "</b>")
    fig = make_subplots(rows=num_rows, cols=num_cols, subplot_titles=subplot_titles, 
                        y_title='Fraction clones',
                        x_title="POD at sample",
                        vertical_spacing=0.03,
                        horizontal_spacing=0.01,
                        shared_yaxes=True, shared_xaxes=True
                        )
    for i in fig['layout']['annotations']:
        i['font'] = dict(size=28, family='Calibri', color='black')
        
    col_number = 1
    row_number = 1   
    # =============================================================================
    # Make graphs
    # =============================================================================
    all_x_values = []
    for tissue in ['PBMC', 'Ileum_allograft']:
        
        if tissue == 'PBMC':
            str_contains = 'PBMC'
        elif tissue == 'Ileum_allograft':
            str_contains = 'Ileum'
        tissue_df = df.loc[df['tissue'].str.contains(str_contains)]
           
        linear_regression_x_values = []
        linear_regression_y_values = []
        for patient_name in subjects_dict_by_age.keys():

            patient_df = tissue_df.loc[tissue_df['patient'] == patient_name]
            patient_tissue_list = patient_df['tissue'].unique().tolist()
            patient_pod_list = patient_df['pod'].unique().tolist()
                         
            y_values = patient_df['percent_mutation_freq_greater_than_2'].tolist()
            x_values = patient_df['pod'].tolist()
            tissues = patient_df['tissue'].tolist()
            linear_regression_x_values += x_values
            linear_regression_y_values += y_values
            all_x_values += x_values
            
            marker_line_color = ['black' for x in tissues]
            marker_color = [tissue_color_dict[tissue] if 'pretransplant' not in x else 'green' for x in tissues]
            marker_line_width = [1 for x in tissues]
            marker_symbol = [subjects_dict_shapes[patient_name] for x in tissues]
            hovertext = ['Patient: ' + patient_name + '<br>Tissue: ' + y for y in tissues]

            fig.add_trace(go.Scatter(
                x=x_values,
                y=y_values,
                # name=tissue,
                marker_color=marker_color,
                marker_symbol=marker_symbol,
                marker_size=12,
                marker_line_width=marker_line_width,
                marker_line_color=marker_line_color,
                hovertext=hovertext,
                mode='markers',
                opacity=0.8,
                showlegend=False
            ),row=row_number, col=col_number)
            
            """
        # =============================================================================
        # calculate linear regression https://realpython.com/linear-regression-in-python/
        # =============================================================================
        y = np.array(linear_regression_y_values)
        x = np.array(linear_regression_x_values)
        
        if len(y) > 0:
            
            # fit a fifth degree polynomial to the data
            # define the true objective function
            def objective(x, a, b, c, d, e, f):
                	return (a * x) + (b * x**2) + (c * x**3) + (d * x**4) + (e * x**5) + f
             
            # curve fit
            popt, _ = curve_fit(objective, x, y)
            # summarize the parameter values
            a, b, c, d, e, f = popt
            # define a sequence of inputs between the smallest and largest known inputs
            x_line = arange(min(x), max(x), 1)
            # calculate the output for the range
            y_line = objective(x_line, a, b, c, d, e, f)

            fig.add_trace(go.Scatter(
                x=x_line,
                y=y_line,
                name='Best Fit: ' + tissue,
                legendgroup='Best fit line',
                mode='lines',
                # marker_color=tissue_color_dict[tissue],
                line=dict(color='black', width=1),
                # marker_symbol=subjects_dict_shapes[patient],
                # marker_size=8,
                showlegend=False,
                # marker=dict(symbol='square')
            ),row=row_number, col=col_number)
            
            """
            
        # =============================================================================
        # Add median for Ileum and BLood
        # =============================================================================
        if tissue == 'PBMC':
            adult_tissue = 'PBL'
        elif tissue == 'Ileum_allograft':
            adult_tissue = 'Ileum'
            
        adult_tissue_df = adult_df.loc[adult_df['tissue'] == adult_tissue]
        adult_median = statistics.median(adult_tissue_df['percent_mutation_freq_greater_than_2'].tolist())
            
        fig.add_trace(go.Scatter(
            x=linear_regression_x_values,
            y=[adult_median for x in linear_regression_x_values],
            name='Adult Median',
            legendgroup='Adult Median',
            mode='lines+text',
            line=dict(color='black', width=1, dash='dot'),
            showlegend=False,
            opacity=0.8
        ),row=row_number, col=col_number)
               
        col_number += 1    

    # =============================================================================
    # Add legend               
    # =============================================================================
    fig.add_trace(go.Scatter(
        x=[None],
        y=[None],
        name='Pretransplant',
        marker_color='green',
        mode='markers',
        ), 
        row=1, col=1)     

    fig.add_trace(go.Scatter(
        x=[None],
        y=[None],
        name='Adult Median',
        legendgroup='Adult Median',
        mode='lines+text',
        line=dict(color='black', width=1, dash='dot'),
        textfont_size=8,
        showlegend=True,
        opacity=0.8
    ),row=1, col=1)
    
    fig.add_trace(go.Scatter(
        x=[None],
        y=[None],
        name='Best fit line',
        legendgroup='Best fit line',
        mode='lines',
        line=dict(color='black', width=1,),
        showlegend=True,
        opacity=0.8
    ),row=1, col=1)          
                                     
    fig.update_yaxes(range=[-0.02, 1.02], tickfont_size=28) 
    fig.update_xaxes(range=[min(all_x_values) - 300, max(all_x_values) + 300], nticks=5, tickfont_size=28) 
 
    fig.update_layout(title_text="Fraction clones with a mutation frequency >2% per time point.",
                      titlefont_size=12,
                      font=dict(
                            family="Calibri",
                            color="black"
                           ),
                      plot_bgcolor='whitesmoke',    
    )

    pio.write_html(fig, file=r'/path_to_output_save_folder/Supp_Fig_4.html') #produces an interactive .html graph. You may want to comment out lines "autosize=False, width=some_number, height=some_number" in fig.update_layout to allow the size to be automatically made for this file.
    fig.write_image('/path_to_output_save_folder/Supp_Fig_4.png', scale=8) #produces a high res .png image
    print("Finished!")
    
make_graph(df_original, adult_df_original)
