"""
This script gets all of the clone_ids that have a minimum number of nodes, from node_min 1 to node_min 50.

Node_min 1 means all clones, because all clones have at least 1 node.

Input: The clones_table filtered for a specific subject. If you used the preprocessing script to download and save the tables as pickle files, this step has already been done and you simply need to iterate through the subjects to load each table one by one.
Output: dict[node_min] = set(clone_ids that meet the node_min)

"""

import json
import pickle
import pandas as pd
from pathlib import Path, PureWindowsPath

lp16_subjects_dict = {
                2: 'D145',
                3: 'D149',
                5: 'D168',
                6: 'D181',
                7: 'D182',
                8: 'D207', 
                    }
    
def get_clones_per_node_min():
        
    for patient_ID, patient_name in lp16_subjects_dict.items():
        
        print("Making list of clones for: " + lp16_subjects_dict[patient_ID])
        
        clones_table = pd.read_pickle(r"/path_to_lp16_data_folder/clones_table/{}_clones_table.pkl".format(patient_ID))
            
        results_1 = {}
        for node_min in range(1, 51):
            results_1[node_min] = []
            
        for clone, clone_id in zip(clones_table['tree'], clones_table['id']):
            if clone != None:
                clone = json.loads(clone)    
                node_counter = {'count': 0}
                def traverse_tree(tree):
                    if len(tree['data']['seq_ids']) > 0:
                        node_counter['count'] += 1
                    for child in tree['children']:
                        traverse_tree(child)                
                traverse_tree(clone['tree'])                   
    
                num_nodes = node_counter['count']
                for node_min in range(1, 51):
                    if num_nodes >= node_min:
                        results_1[node_min].append(clone_id)
         
        results_2 = {}
        for node_min, clone_list in results_1.items():
            temp = set(clone_list)
            results_2[node_min] = temp
            
        # =============================================================================
        # Save results    
        # =============================================================================          
        path = Path(r"/path_to_lp16_data_folder/{}_clones_per_node_min_using_tree.pkl".format(lp16_subjects_dict[patient_ID]))
        pickle_out = open(path,"wb")
        pickle.dump(results_2, pickle_out)
        pickle_out.close() 
         
    print("Finished!")
                      
get_clones_per_node_min()



