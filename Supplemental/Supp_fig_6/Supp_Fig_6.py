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
import subprocess
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
                'Pt21_R' : 2.6,
                'Pt19_R' : 3.3,
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
        path = Path(r"/path_to_sykes_children_data_folder/{}_clumpiness_POD_Specific.pkl".format(patient_name))
        pickle_in = open(path,"rb")
        clumpiness_data = pickle.load(pickle_in)
        pickle_in.close()
  
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

def make_fig():
  
    # =============================================================================
    # Set the minimum number of nodes per clone which will be used to filter clones.  
    # Note: You cannot set node_min to < 3 because the preprocessing script only gets clumpiness for 3+ node clones. You will have to change the preprocessing script 
    # if you want clones with fewer than 3 nodes. However, it is pointless to run clumpiness on clones with fewer than 3 nodes due to overly simple lineages skewing the results.
    # =============================================================================
    node_min = 3
          
    fig = go.Figure()
    
    patients = list(subjects_dict_children.keys())
    first_tissue = 'PBMC'
    second_tissue = 'Ileum_allograft'
    third_tissue = 'Ileum'
           
    linear_regression_x_values = []
    linear_regression_y_values = []
    for patient in patients:
        
        x_values = []
        y_values = []
        marker_symbol = []
        marker_color = []
        marker_line_color = []
        marker_line_width = []
        hovertext = []            
 
        tissue_clones, node_min_clone_dict, clumpiness_data, rejection_clones = load_data(patient)   
        node_min_clones = node_min_clone_dict[node_min]
        
        patient_PODS = list(clumpiness_data.keys())
        patient_PODS_ints = []
        for POD in patient_PODS:
            POD_list = eval(POD.split("POD")[1])
            patient_PODS_ints.append(POD_list[0])
        
        for POD_str, POD_int in zip(patient_PODS, patient_PODS_ints):
            POD_data = clumpiness_data[POD_str]
            tissue_pairs = set(POD_data.keys())
            if 'Ileum_allograft+PBMC' not in tissue_pairs  \
            and 'PBMC+Ileum_allograft' not in tissue_pairs \
            and 'PBMC+Ileum' not in tissue_pairs \
            and 'Ileum+PBMC' not in tissue_pairs:
                continue
    
            #Find the order of the tissue names
            if first_tissue + "+" + second_tissue in tissue_pairs:
                current_tissue = first_tissue + "+" + second_tissue
            elif second_tissue + "+" + first_tissue in tissue_pairs:
                current_tissue = second_tissue + "+" + first_tissue
            elif third_tissue + "+" + first_tissue in tissue_pairs:
                current_tissue = third_tissue + "+" + first_tissue
            elif first_tissue + "+" + third_tissue in tissue_pairs:
                current_tissue = first_tissue + "+" + third_tissue
                                 
            POD_values = []
            for clone_id, clumpiness_values in POD_data[current_tissue].items():
                if int(clone_id) in node_min_clones and int(clone_id) not in rejection_clones:
                    value = clumpiness_values["".join([x if x != "+" else "," for x in current_tissue])]
                    POD_values.append(value)
                
            total_clones = len(POD_values)
            if total_clones > 5:    
                median = statistics.median(POD_values)
                y_values.append(median)
                x_values.append(POD_int)
                marker_symbol.append(subjects_dict_shapes[patient])
                marker_color.append(tissue_color_dict[first_tissue])
                marker_line_color.append(tissue_color_dict[second_tissue])  
                marker_line_width.append(2)                      
                hovertext.append("Patient: {}<br>Clones: {}<br>Tissue: {}".format(patient, total_clones, current_tissue))
                
                linear_regression_x_values.append(POD_int)
                linear_regression_y_values.append(median)                  
                                
        fig.add_trace(go.Scatter(
            x=x_values, 
            y=y_values,
            mode='markers+lines',
            marker_size=16,        
            marker_color=marker_color,
            marker_line_color=marker_line_color,
            marker_line_width=marker_line_width,                
            marker_symbol=marker_symbol,
            hovertext=hovertext,
            opacity=0.85,
            showlegend=False,
        ))
        
    # =============================================================================
    # Add adult PPG median line       
    # =============================================================================
    adult_median_results = []
    adult_name_list = []
    for patient in lp16_subjects_dict.values():
        
        current_tissue = 'Ileum+PBL'
        tissue_clones, node_min_clone_dict, clumpiness_data, rejection_clones = load_data(patient)
                
        node_min_clones = node_min_clone_dict[node_min]

        if current_tissue not in clumpiness_data.keys():
            continue
        
        current_values = []
        for clone_id, clumpiness_values in clumpiness_data[current_tissue].items():
            if int(clone_id) in node_min_clones:
                value = clumpiness_values["".join([x if x != "+" else "," for x in current_tissue])]
                current_values.append(value)
            
        total_clones = len(current_values)
        if total_clones > 5:
            median = statistics.median(current_values)    
            adult_median_results.append(median)
            adult_name_list.append(patient) #just in case one wants to see which patient the median comes from
         
    median_of_adult_medians = statistics.median(adult_median_results)
    y_values = [median_of_adult_medians, median_of_adult_medians]
    x_values = [-.25, max(linear_regression_x_values) + 0.25]
    fig.add_trace(go.Scatter(
        x=x_values, 
        y=y_values,
        mode='lines',
        line_color='black',
        line_dash='dash',
        opacity=0.85,
        name='Median of adult data',
        showlegend=False,
    ))
               
    fig.update_yaxes(range=[0, 1.1], tickfont_size=28, nticks=4)

    fig.update_layout(title_text="Median clumpiness by POD. Node min: " + str(node_min),
                      yaxis_title="Median clumpiness",
                      xaxis_title="POD",
                      font=dict(
                            family="Calibri",
                            size=28,
                            color="black"
                           ),
                      plot_bgcolor='whitesmoke',
    
    )
    
    pio.write_html(fig, file=r'/path_to_output_save_folder/Supp_Fig_6.html')
    fig.write_image('/path_to_output_save_folder/Supp_Fig_6.png', scale=8)
    print("Finished!")
        
make_fig()  
