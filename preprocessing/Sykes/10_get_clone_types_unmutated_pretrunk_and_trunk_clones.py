"""
Summary:
    This script simplifies the results of the trunk length script. 
    This script removes the extra information from the trunk length script and simply returns the clones that belong to the following groups:
        result['unmutated'] = set(list of clone_ids...)
        result['pretrunk'] = set(list of clone_ids...)
        result['trunk'] = set(list of clone_ids...)

Prerequisite scripts:
    Script that calculates trunk length.
"""

import json
import pickle
import pandas as pd
from pathlib import Path, PureWindowsPath
import os



# =============================================================================
# Load data.
# =============================================================================
path = Path(r"/path_to_sykes_children_data_folder/trunk_length_simple_method_all_clones_with_tissues_combined.pkl")
pickle_in = open(path,"rb")
data = pickle.load(pickle_in)
pickle_in.close()
    
def get_clone_types(df):
    
    subjects_dict = {   
                    2: 'Pt19_R',
                    3: 'Pt21_R',
                    4: 'Pt23_R',
                    8: 'Pt14_R',
                    9: 'Pt20_R',
                    10: 'Pt17_R',
                    17: 'Pt25_R',
                    }

    for patient in subjects_dict.values():
        print("Getting clone types for: " + patient)
        patient_data = data[patient]
    
        result = {}
        for clone_type, clones in patient_data.items():
            clones = set(clones.keys())
            result[clone_type] = clones
            
        # =============================================================================
        # Save the result as a pickle    
        # =============================================================================
        # path = Path(r"/home/user/Tommy/sykesAnchoring2020_updated1/{}_unmutated_mutated_and_trunk_clones.pkl".format(patient))
        path = Path(r"/path_to_sykes_children_data_folder/{}_unmutated_pretrunk_and_trunk_clones.pkl".format(patient))
        pickle_out = open(path,"wb")
        pickle.dump(result, pickle_out)
        pickle_out.close() 
            
    print("Finished!")
            
get_clone_types(data)




