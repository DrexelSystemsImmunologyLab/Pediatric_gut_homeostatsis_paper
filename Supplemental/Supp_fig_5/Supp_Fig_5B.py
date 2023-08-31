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

def load_clone_type_data(patient_name):
    
    if patient_name in subjects_dict.values():
        # =============================================================================
        # Load result. 
        # =============================================================================
        path = Path(r"/path_to_sykes_children_data_folder/{}_unmutated_pretrunk_and_trunk_clones.pkl".format(patient_name))
        pickle_in = open(path,"rb")
        clone_type_data = pickle.load(pickle_in)
        pickle_in.close()
        
    elif patient_name in lp16_subjects_dict.values():
        # =============================================================================
        # Load result. 
        # =============================================================================
        path = Path(r"/path_to_lp16_data_folder/{}_unmutated_pretrunk_and_trunk_clones.pkl".format(patient_name))
        pickle_in = open(path,"rb")
        clone_type_data = pickle.load(pickle_in)
        pickle_in.close()

    return clone_type_data

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
        temp = [x + random.uniform(-0.10, 0.10) for x in x_values]
        list_of_coords = [(x, y) for x, y in zip(temp, y_values)]
        min_distance = min([distance(*combo) for combo in combinations(list_of_coords,2)])
        num_tries += 1
        if min_distance > current_best:
            current_best = min_distance
            result = temp
    return result

def make_fig():

    # =============================================================================
    # Set the minimum number of nodes per clone which will be used to filter clones.    
    # =============================================================================
    node_min = 3
      
    fig = go.Figure() 
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
    
    for x_axis_name, x_axis_value in label_order.items():
        
        x_values = []
        y_values = []
        marker_symbol = []
        marker_color = []
        marker_line_color = []
        marker_line_width = []
        hovertext = []            
        current_tissue = x_axis_name.split(" ")[0]
        first_tissue = current_tissue.split("+")[0]
        second_tissue = current_tissue.split("+")[1]   
        
        if "child" in x_axis_name:
            
            for patient in children:
                
                tissue_clones, node_min_clone_dict, clumpiness_data, rejection_clones = load_data(patient)
                clone_type_data = load_clone_type_data(patient)
                trunk_clones = clone_type_data['trunk']
                
                if first_tissue + "+" + second_tissue not in clumpiness_data.keys() and second_tissue + "+" + first_tissue not in clumpiness_data.keys():
                    continue
                                
                #Find the order of the tissue names
                if first_tissue + "+" + second_tissue in clumpiness_data.keys():
                    current_tissue = first_tissue + "+" + second_tissue
                elif second_tissue + "+" + first_tissue in clumpiness_data.keys():
                    current_tissue = second_tissue + "+" + first_tissue
                                     
                clumpiness_clones = set([int(x) for x in clumpiness_data[current_tissue].keys()])
                node_min_clones = node_min_clone_dict[node_min]
                total_clones_set = clumpiness_clones.intersection(node_min_clones, trunk_clones)
                
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
                        
                    x_values.append(x_axis_value)
                    y_values.append(median)
                    marker_symbol.append(subjects_dict_shapes[patient])
                    marker_color.append(tissue_color_dict[first_tissue]) 
                    marker_line_color.append(tissue_color_dict[second_tissue])  
                    marker_line_width.append(2)                      
                    hovertext.append("Patient: {}<br>Clones: {}<br>Total clones: {}<br>Tissue: {}".format(patient, len(temp), total_clones, current_tissue))
        
        elif "adult" in x_axis_name:
            
            if first_tissue == 'PBMC':
                first_tissue = 'PBL'
                current_tissue = first_tissue + '+' + second_tissue
            
            for patient in adults:
                             
                tissue_clones, node_min_clone_dict, clumpiness_data, rejection_clones = load_data(patient)
                clone_type_data = load_clone_type_data(patient)
                trunk_clones = clone_type_data['trunk']
                
                if first_tissue + "+" + second_tissue not in clumpiness_data.keys() and second_tissue + "+" + first_tissue not in clumpiness_data.keys():
                    print("Did not find tissue pair in {}".format(patient))
                    continue
                                
                #Find the order of the tissue names
                if first_tissue + "+" + second_tissue in clumpiness_data.keys():
                    current_tissue = first_tissue + "+" + second_tissue
                elif second_tissue + "+" + first_tissue in clumpiness_data.keys():
                    current_tissue = second_tissue + "+" + first_tissue
                                     
                clumpiness_clones = set([int(x) for x in clumpiness_data[current_tissue].keys()])
                node_min_clones = node_min_clone_dict[node_min]
                total_clones_set = clumpiness_clones.intersection(node_min_clones, trunk_clones)
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
                        
                    x_values.append(x_axis_value)
                    y_values.append(median)
                    marker_symbol.append(subjects_dict_shapes[patient])
                    marker_color.append('black')
                    marker_line_color.append('black')  
                    marker_line_width.append(2)                      
                    hovertext.append("Patient: {}<br>Clones: {}<br>Total clones: {}<br>Tissue: {}".format(patient, len(temp), total_clones, current_tissue))

        if len(x_values) == 0:
           continue
         
        if len(x_values) > 1:
            x_values = add_jitter(x_values, y_values)
        
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
            
    fig.update_yaxes(range=[0, 1.1], tickfont_size=28, nticks=5)
    fig.update_xaxes(tickangle=25,
                      tickvals=list(label_order.values()),
                      ticktext=list(label_order.keys()),
                      tickfont_size=28,                     
                      )          
     
    fig.update_layout(title_text="Median clumpiness. Node min: " + str(node_min),
                      yaxis_title="Median clumpiness",
                      font=dict(
                            family="Calibri",
                            size=28,
                            color="black"
                           ),
                       # margin_t=300,
                        # legend={'itemwidth': '60',},
                       # legend_itemwidth=50,                  
                        width=700,
                        height=550,
                      plot_bgcolor='whitesmoke', 
    )
    
    pio.write_html(fig, file=r'/path_to_output_save_folder/Supp_Fig_5B.html') #produces an interactive .html graph. You may want to comment out lines "autosize=False, width=some_number, height=some_number" in fig.update_layout to allow the size to be automatically made for this file.
    fig.write_image('/path_to_output_save_folder/Supp_Fig_5B.png', scale=8) #produces a high res .png image
    print("Finished!")
        
make_fig()  
