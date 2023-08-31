# -*- coding: utf-8 -*-
"""
Created on Tue Sep 22 23:20:01 2020

@author: thsia 
"""

import plotly.io as pio
import plotly.graph_objects as go

subjects_dict_by_age = {   
                'Pt20_R' : 1.7,
                'Pt23_R' : 2.2,
                'Pt14_R' : 2.3,
                'Pt21_R' : 2.6,
                'Pt19_R' : 3.3,
                'Pt17_R' : 5.4,     
                'Pt25_R': 9.32,

                'D207' : 23,               
                'D181' : 46,
                'D182' : 46,
                'D149' : 55,
                'D168' : 56,
                'D145' : 58,
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
        

def make_age_legend():

    fig = go.Figure()

    for patient_name, age in subjects_dict_by_age.items():
                
        fig.add_trace(go.Scatter(
            x=[None],
            y=[None],
            name="<b>" + patient_name + ": </b>" + str(age),
            marker_size=8,
            mode='markers',
            marker_color='black',
            marker_symbol=subjects_dict_shapes[patient_name],
            marker_line_width=1,                
            ), 
            )
                           
    fig.update_layout(
                    font=dict(size=12, family='Calibri', color='black'),
                    title='Subjects by age legend.',
                    plot_bgcolor='whitesmoke',
                    barmode='group',
                    height=500,
                    width=800,
                    # margin_t=270,
                    )
    
    pio.write_html(fig, file=r'/path_to_output_save_folder/Fig4_name_symbol_age_legend.html') #Saves legend as .html graph
    fig.write_image('/path_to_output_save_folder/Fig4_name_symbol_age_legend.png', scale=8) #Saves legend as .png image
    print("Finished!")
    
make_age_legend()

