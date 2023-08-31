"""

NOTE:
    You will have to install the Haskell program find-clumpiness on your machine before running this script.
    For more info, see: https://github.com/GregorySchwartz/find-clumpiness
    Also, this script calls the find-clumpiness program using the terminal via linux commands. The commands may not work if you are on Windows.
    Lastly, this script requires you to make a temp folder, where you will not care if all of the files in that folder get deleted.
    The script uses the temp folder to store the JSON strings of the lineages and clumpiness results, which are deleted after each set of clones is done.
    
    8_get_clones_in_tissue_pairs.py is a prerequisite preprocessing script for this script. Its output is loaded on line 68.

The script then feeds in a list of clones into a function that reformats their trees into a format that can be used by find-clumpiness' "JSON" parser. The end
result of this step produces a file for each clone's tree that can now be used with find-clumpiness.

Next, find-clumpiness is run on the files that were just made, outputting a .txt file for each corresponding file with the results.

Finally, the script builds a dictionary that collects all of the find-clumpiness results by iterating through all of the .txt result files.
Some of the clones produce empty results or do not produce results for the selected criteria, such as clones with only a single node or where the
tissues being compared are only found in a single node. The script then removes all of these empty values from the dictionary so that the end result
has a dictionary with each key (clone) having meaningful results that can be used for plotting. This dictionary is then saved as a .JSON file.

The original find-clumpiness results are in the following format:

    property1,property2,value
    PBMC,PBMC,1.0
    PBMC,MLN,1.0
    MLN,PBMC,1.0
    MLN,MLN,1.0

and this script turns it into:

    results[94] : {'PBMC,PBMC': 1.0, 'PBMC,MLN': 1.0, 'MLN,PBMC': 1.0, 'MLN,MLN': 1.0}    
    
with 94 being the 95th (starts from 0) clumpiness result file. There is 1 file for each clone.
    
"""

import json
import pickle
import pandas as pd
from pathlib import Path, PureWindowsPath
import os
import re
import subprocess

subjects_dict = {   
                2: 'Pt19_R',
                3: 'Pt21_R',
                4: 'Pt23_R',
                8: 'Pt14_R',
                9: 'Pt20_R',
                10: 'Pt17_R',
                17: 'Pt25_R',
                }
    
def get_clumpiness():
         
    # =============================================================================
    # Load tables          
    # =============================================================================
    def load_data_pretransplant_version(patient_name, patient_ID):
        if patient_name in subjects_dict.values():

            clones_table = pd.read_pickle(r"/path_to_sykes_children_data_folder/clones_table/{}_clones_table.pkl".format(patient_ID))
            clone_stats = pd.read_pickle(r"/path_to_sykes_children_data_folder/clone_stats/{}_clone_stats.pkl".format(patient_ID))
            functional_clones = clone_stats['clone_id'].loc[clone_stats['functional'] == 1].unique().tolist()
            clones_table = clones_table.loc[clones_table['id'].isin(functional_clones)]
            tissue_pair_clones_dict = pd.read_pickle(r"/path_to_sykes_children_data_folder/{}_clones_per_tissue_pairs_new_1.pkl".format(patient_name))

        return clones_table , tissue_pair_clones_dict       
    # =============================================================================
    # This function returns two dicts:
    # samples_per_tissue: dict[tissue] : [list of sample_ids in that tissue]    
    # tissues_per_sample: dict[sample_id] : [tissue in that sample_id] ###this will be a list with a single element (the tissue)
    # =============================================================================
    def get_tissue_sample_ids(patient_name):
      
        tissues_per_sample = {}
        samples_per_tissue = {}
        
        if patient_name in subjects_dict.values():
            samples_table = pd.read_pickle(r"/path_to_sykes_children_data_folder/samples.pkl")
            metadata_table = pd.read_pickle(r"/path_to_sykes_children_data_folder/sample_metadata.pkl")
            subject_samples = samples_table['id'].loc[samples_table['subject_id'] == patient_ID].unique() 
            metadata_table = metadata_table.loc[metadata_table['sample_id'].isin(subject_samples)]
            metadata_table['value'] = metadata_table['value'].replace("Colon ", "Colon")
            metadata_table['value'] = metadata_table['value'].replace("MLN_allograft", "AxLN")
            metadata_table['value'] = metadata_table['value'].replace("POD37 44", "POD37")
            ###Change Ileum to Ileum_allograft for Pt21_R POD1145
            POD1145_sample_ids = metadata_table['sample_id'].loc[metadata_table['value'] == 'POD1145']
            metadata_table['value'].loc[metadata_table['sample_id'].isin(POD1145_sample_ids)] =  metadata_table['value'].loc[metadata_table['sample_id'].isin(POD1145_sample_ids)].replace("Ileum", "Ileum_allograft")
            sample_ids = metadata_table['sample_id'].unique()
            tissues = metadata_table['value'].loc[metadata_table['key'] == 'sample_origin'].unique()
                     
            for sample_id in sample_ids:
                tissues_per_sample[sample_id] = list(metadata_table['value'].loc[(metadata_table['key'] == 'sample_origin') & (metadata_table['sample_id'] == sample_id)])[0]
            for tissue in tissues:
                samples_per_tissue[tissue] = list(metadata_table['sample_id'].loc[metadata_table['value'] == tissue].unique())
                
            return tissues_per_sample, samples_per_tissue, tissues
        
    for patient_ID, patient_name in subjects_dict.items():
            
        ###IMPORTANT: PATH TO TEMP FOLDER WHERE ALL FILES WILL BE DELETED
        ###change directory to the directly where find-clumpiness results will be saved
        temp_folder_where_files_will_be_stored_and_deleted = "/path_to_sykes_children_data_folder/temp"
        
        os.chdir(temp_folder_where_files_will_be_stored_and_deleted) 
        subprocess.run(r'find . -type f -delete', shell=True, check=True, stdout=subprocess.PIPE, universal_newlines=True)
    
        print("Making list of clones for: " + patient_name)
    
        clones_table, tissue_pair_clones_dict = load_data_pretransplant_version(patient_name, patient_ID)
        tissues_per_sample, samples_per_tissue, tissues = get_tissue_sample_ids(patient_name)
      
        results_1 = {}
        
        ###Make list of clones_ids to be analyzed
        for tissue1 in range(0, len(tissues) - 1):
            
            for tissue2 in range(tissue1 + 1, len(tissues)):
                tissue_pair_clones = tissue_pair_clones_dict[tissues[tissue1] + "+" + tissues[tissue2]]
                list_of_clones_result = tissue_pair_clones
                
                print("Finished making list of clones. Now making files: " + str(len(list_of_clones_result)) + " total files.")
            
                # =============================================================================
                # Iterate through each tree, reformatting it and making files for find-clumpiness' "JSON" parser            
                # =============================================================================    
                clones_table_1 = clones_table.loc[clones_table['id'].isin(list_of_clones_result)]
    
                for clone, clone_id in zip(clones_table_1['tree'], clones_table_1['id']):
                        
                    if clone == None:
                        continue
                    
                    def reformat_trees_for_clumpiness(input):
                    	result = []
                    	for key, value in input.items():
                    		if key == "data":
                    			result.append(parse_data(value))
                    		if key == "children":
                    			if len(value) > 0:
                    				childrenResult = []
                    				for child in value:
                    					childrenResult.append(reformat_trees_for_clumpiness(child))
                    				result.append(childrenResult)
                    			else:
                    				result.append([])
                    	return result
                    
                    def parse_data(data):
                        result = {"nodeID" : "ID", "nodeLabels" : []}
                        seq_ids = data['seq_ids']
                        for key, seq_id in seq_ids.items():
                            metadata = seq_id['metadata']
                            if 'Colon ' in metadata['sample_origin']:
                                metadata['sample_origin'] = metadata['sample_origin'].replace("Colon ", "Colon")
                            if 'Ileum' in metadata['sample_origin'] and metadata['pod'] == 'POD1145':
                                metadata['sample_origin'] = metadata['sample_origin'].replace("Ileum", "Ileum_allograft")
                            result['nodeLabels'].append(metadata['sample_origin'])
                        return result
                    
                    json_data = json.loads(clone)
                    result = reformat_trees_for_clumpiness(json_data["tree"])

                    path = Path(r"{}/{}.JSON".format(temp_folder_where_files_will_be_stored_and_deleted, clone_id))
                    with open(path, 'w') as outfile:
                        json.dump(result, outfile)
                                          
                print("Finished making files. Now running find-clumpiness.")
        
                # =============================================================================
                # This command runs find-clumpiness from the terminal on the files produced above        
                # =============================================================================
                subprocess.run(r'for file in *.JSON; do cat ${file%%.*}.JSON | find-clumpiness --format "JSON" -e Majority -o ${file%%.*}.txt; done', shell=True, check=True, stdout=subprocess.PIPE, universal_newlines=True)
        
                print("find-clumpiness finished. Now collecting results as a dictionary.")
               
                ###this directory contains the find-clumpiness results files
                directory = temp_folder_where_files_will_be_stored_and_deleted
                
                # =============================================================================
                # Collect the results from the find-clumpiness output files as a dictionary      
                # =============================================================================
                file_results = {}
                for entry in os.scandir(directory):
                    if entry.path.endswith(".txt") and entry.is_file():
                        with open(entry) as reader:
                            file_results[entry.name.split(".")[0]] = {}
                            line_count = 1
                            for line in reader.readlines():
                                line_count += 1
                                if 'property1,property2,value' not in line:
                                    key_value = line.replace("," + str(re.findall("\d+\.\d+", line)[0]) + '\n', "")
                                    file_results[entry.name.split(".")[0]][key_value] = float(re.findall("\d+\.\d+", line)[0])
                
                ###remove empty dictionaries, which can be due to trees with only a single node not returning any values from find-clumpiness or the labels all being in a single node. 
                final_results = file_results.copy()
                node_removed = 0
                current_tissues = str(tissues[tissue1] + "," + tissues[tissue2])
                for key in file_results:
                    if current_tissues not in file_results[key]:
                        del final_results[key]
                        node_removed += 1
                        
                print("Finished with " + tissues[tissue1] + "+" + tissues[tissue2] + ". Removed " + str(node_removed) + " results where there were no values.")
                print("Now deleting find-clumpiness files to prepare for next iteration.")
        
                # =============================================================================
                # Delete all of the files that were created this iteration via the terminal
                # =============================================================================
                subprocess.run(r'find . -type f -delete', shell=True, check=True, stdout=subprocess.PIPE, universal_newlines=True)
              
                results_1[tissues[tissue1] + "+" + tissues[tissue2]] = final_results
                
        # =============================================================================
        # Save results    
        # =============================================================================          
        path = Path(r"/path_to_sykes_children_data_folder/{}_clumpiness_all_tissue_pair_clones_updated_node_min_script.pkl".format(patient_name))
        pickle_out = open(path,"wb")
        pickle.dump(results_1, pickle_out)
        pickle_out.close() 
                 
    print("Finished!")
          
get_clumpiness()



