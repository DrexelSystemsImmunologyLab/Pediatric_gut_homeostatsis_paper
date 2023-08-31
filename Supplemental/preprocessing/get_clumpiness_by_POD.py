"""

NOTE:
    You will have to install the Haskell program find-clumpiness on your machine before running this script.
    For more info, see: https://github.com/GregorySchwartz/find-clumpiness
    Also, this script calls the find-clumpiness program using the terminal via linux commands. The commands may not work if you are on Windows.
    Lastly, this script requires you to make a temp folder, where you will not care if all of the files in that folder get deleted.
    The script uses the temp folder to store the JSON strings of the lineages and clumpiness results, which are deleted after each set of clones is done.

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
    
def get_clumpiness_by_POD():
    
    # =============================================================================
    # Load tables          
    # =============================================================================
    def load_data(patient_ID, patient_name):
        if patient_name in subjects_dict.values():
            samples_table = pd.read_pickle(r"/path_to_sykes_children_data_folder/samples.pkl")
            metadata_table = pd.read_pickle(r"/path_to_sykes_children_data_folder/sample_metadata.pkl")    
            metadata_table['value'] = metadata_table['value'].replace("Colon ", "Colon")
            metadata_table['value'] = metadata_table['value'].replace("MLN_allograft", "AxLN")
            metadata_table['value'] = metadata_table['value'].replace("POD37 44", "POD37")
            ###Change Ileum to Ileum_allograft for Pt21_R POD1145
            POD1145_sample_ids = metadata_table['sample_id'].loc[metadata_table['value'] == 'POD1145']
            metadata_table['value'].loc[metadata_table['sample_id'].isin(POD1145_sample_ids)] =  metadata_table['value'].loc[metadata_table['sample_id'].isin(POD1145_sample_ids)].replace("Ileum", "Ileum_allograft")
            
            df_of_trees = pd.read_pickle(r"/path_to_sykes_children_data_folder/clones_table/{}_clones_table.pkl".format(patient_ID))
            seq_table = pd.read_pickle(r"/path_to_sykes_children_data_folder/sequences_table/{}_sequences_table.pkl".format(patient_ID))
            seq_table.dropna(subset=["clone_id"], inplace=True)       
            seq_table['clone_id'] = seq_table['clone_id'].astype(int)
            
            # =============================================================================
            # Load rejection samples and remove these samples from data
            # =============================================================================
            path = Path(r"/path_to_sykes_children_data_folder/{}_samples_with_rejection.pkl".format(patient_name))
            pickle_in = open(path,"rb")
            rejection_samples = pickle.load(pickle_in)
            pickle_in.close()
            seq_table = seq_table.loc[~seq_table['sample_id'].isin(rejection_samples)]
            metadata_table = metadata_table.loc[~metadata_table['sample_id'].isin(rejection_samples)]
            samples_table = samples_table.loc[~samples_table['id'].isin(rejection_samples)]
            
            clone_stats = pd.read_pickle(r"/path_to_sykes_children_data_folder/clone_stats/{}_clone_stats.pkl".format(patient_ID))
            functional_clones = clone_stats['clone_id'].loc[clone_stats['functional'] == 1].unique().tolist()
            seq_table = seq_table.loc[seq_table['clone_id'].isin(functional_clones)]
            
            subject_samples = samples_table['id'].loc[samples_table['subject_id'] == patient_ID].unique()   
            subject_metadata_table = metadata_table.loc[metadata_table['sample_id'].isin(subject_samples)]        
            tissues = subject_metadata_table['value'].loc[subject_metadata_table['key'] == 'sample_origin'].unique()
            node_min_clone_dict = pd.read_pickle(r"/path_to_sykes_children_data_folder/{}_clones_per_node_min_using_tree.pkl".format(patient_name))
  
        return df_of_trees, seq_table, tissues, node_min_clone_dict, rejection_samples
            
    
    # =============================================================================
    # POD grouping function to group PODs together that are within 2 days of each other
    # Input is name of patient, such as "Pt19_R"
    # Output is list of lists, with each sublist containing the grouped PODs (or single POD if there are no PODs that are grouped with it)
    # And a second output in the same format, except each sublist is a list of sample_ids that belong to the corresponding POD
    # Example: [[0], [20, 22], [100], [110, 112, 114]], [[sample_ids that belong to 0], [sample_ids that belong to 20, 22], [etc], [etc]]
    # =============================================================================
    def group_PODs(patient_ID, patient_name):
        samples_table = pd.read_pickle(r"/path_to_sykes_children_data_folder/samples.pkl")
        metadata_table = pd.read_pickle(r"/path_to_sykes_children_data_folder/sample_metadata.pkl")
        metadata_table['value'] = metadata_table['value'].replace("POD37 44", 'POD37')
        metadata_table['value'] = metadata_table['value'].replace("Colon ", "Colon")
        metadata_table['value'] = metadata_table['value'].replace("MLN_allograft", "AxLN")
        metadata_table['value'] = metadata_table['value'].replace("POD37 44", "POD37")
        ###Change Ileum to Ileum_allograft for Pt21_R POD1145
        POD1145_sample_ids = metadata_table['sample_id'].loc[metadata_table['value'] == 'POD1145']
        metadata_table['value'].loc[metadata_table['sample_id'].isin(POD1145_sample_ids)] =  metadata_table['value'].loc[metadata_table['sample_id'].isin(POD1145_sample_ids)].replace("Ileum", "Ileum_allograft")
             
        # =============================================================================
        # Load rejection samples and remove these samples from data
        # =============================================================================
        path = Path(r"/path_to_sykes_children_data_folder/{}_samples_with_rejection.pkl".format(patient_name))
        pickle_in = open(path,"rb")
        rejection_samples = pickle.load(pickle_in)
        pickle_in.close()
        metadata_table = metadata_table.loc[~metadata_table['sample_id'].isin(rejection_samples)]
        samples_table = samples_table.loc[~samples_table['id'].isin(rejection_samples)]
                        
        subject_samples = samples_table['id'].loc[samples_table['subject_id'] == patient_ID].unique()
        metadata_table = metadata_table.loc[metadata_table['sample_id'].isin(subject_samples)]
        pod_list_strings = metadata_table['value'].loc[metadata_table['key'] == 'pod'].unique()
        pod_list = [int(''.join(i for i in x if i.isdigit())) for x in pod_list_strings]     
        ###sort both using the int version as the guide
        pod_list, pod_list_strings = (list(t) for t in zip(*sorted(zip(pod_list, pod_list_strings))))
        
        result = []
        result_samples = []
        current_pod = pod_list[0]
        current_samples = metadata_table['sample_id'].loc[metadata_table['value'] == pod_list_strings[pod_list.index(current_pod)]].unique().tolist()
        temp = []
        temp_samples = []
        for next_pod in pod_list[1:]: 
            next_samples = metadata_table['sample_id'].loc[metadata_table['value'] == pod_list_strings[pod_list.index(next_pod)]].unique().tolist()
            if len(temp) == 0:
                temp.append(current_pod)
                for sample in current_samples:
                    temp_samples.append(sample)
            if next_pod - current_pod <= 2:       
                temp.append(next_pod)
                for sample in next_samples:
                    temp_samples.append(sample)                
            else:
                result.append(temp)
                result_samples.append(temp_samples)
                temp_samples = []
                temp = []            
            current_pod = next_pod
            current_samples = next_samples
        if pod_list[-1] not in result[-1]:
            result.append([pod_list[-1]])
            result_samples.append(metadata_table['sample_id'].loc[metadata_table['value'] == pod_list_strings[pod_list.index(pod_list[-1])]].unique().tolist())

        return result, result_samples
        

    # =============================================================================
    # This function returns two dicts:
    # samples_per_tissue: dict[tissue] : [list of sample_ids in that tissue]    
    # tissues_per_sample: dict[sample_id] : [tissue in that sample_id] ###this will be a list with a single element (the tissue)
    # =============================================================================
    def get_tissue_sample_ids(patient_ID, patient_name):
      
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
            
            # =============================================================================
            # Load rejection samples and remove these samples from data
            # =============================================================================
            path = Path(r"/path_to_sykes_children_data_folder/{}_samples_with_rejection.pkl".format(patient_name))
            pickle_in = open(path,"rb")
            rejection_samples = pickle.load(pickle_in)
            pickle_in.close()
            metadata_table = metadata_table.loc[~metadata_table['sample_id'].isin(rejection_samples)]
            samples_table = samples_table.loc[~samples_table['id'].isin(rejection_samples)]
            
            sample_ids = metadata_table['sample_id'].unique()
            tissues = metadata_table['value'].loc[metadata_table['key'] == 'sample_origin'].unique()
                     
            for sample_id in sample_ids:
                tissues_per_sample[sample_id] = list(metadata_table['value'].loc[(metadata_table['key'] == 'sample_origin') & (metadata_table['sample_id'] == sample_id)])[0]
            for tissue in tissues:
                samples_per_tissue[tissue] = list(metadata_table['sample_id'].loc[metadata_table['value'] == tissue].unique())
                
        return tissues_per_sample, samples_per_tissue
    
    def get_timepoint_sample_ids(patient_ID, patient_name):
      
        timepoints_per_sample = {}
        samples_per_timepoint = {}
        
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
            
            # =============================================================================
            # Load rejection samples and remove these samples from data
            # =============================================================================
            path = Path(r"/path_to_sykes_children_data_folder/{}_samples_with_rejection.pkl".format(patient_name))
            pickle_in = open(path,"rb")
            rejection_samples = pickle.load(pickle_in)
            pickle_in.close()
            metadata_table = metadata_table.loc[~metadata_table['sample_id'].isin(rejection_samples)]
            samples_table = samples_table.loc[~samples_table['id'].isin(rejection_samples)]
            
            sample_ids = metadata_table['sample_id'].unique()
            timepoints = metadata_table['value'].loc[metadata_table['key'] == 'pod'].unique().tolist()
            timepoints_original = metadata_table['value'].loc[metadata_table['key'] == 'pod'].unique().tolist()
            timepoints = [int(''.join(i for i in x if i.isdigit())) for x in timepoints]
            timepoints, timepoints_original = (list(t) for t in zip(*sorted(zip(timepoints, timepoints_original))))
            
            for sample_id in sample_ids:
                result = list(metadata_table['value'].loc[(metadata_table['key'] == 'pod') & (metadata_table['sample_id'] == sample_id)])[0]
                result = int(''.join(i for i in result if i.isdigit()))
                timepoints_per_sample[sample_id] = result                
            for timepoint, timepoint_original in zip(timepoints, timepoints_original):
                samples_per_timepoint[str(timepoint)] = list(metadata_table['sample_id'].loc[metadata_table['value'] == timepoint_original].unique())
                
        return timepoints_per_sample, samples_per_timepoint
        
    # =============================================================================
    # Body of script starts here     
    # =============================================================================
    for patient_ID, patient_name in subjects_dict.items():
        
        """
        ###IMPORTANT: BELOW IS THE PATH TO THE TEMP FOLDER WHERE ALL FILES WILL BE DELETED
        """
        ###change directory to the directly where find-clumpiness results will be saved
        temp_folder_where_files_will_be_stored_and_deleted = "/path_to_sykes_children_data_folder/temp"
        os.chdir(temp_folder_where_files_will_be_stored_and_deleted) 
        subprocess.run(r'find . -type f -delete', shell=True, check=True, stdout=subprocess.PIPE, universal_newlines=True)
    
        print("Making list of clones for: " + patient_name)
    
        df_of_trees, seq_table, tissues, node_min_clone_dict, rejection_samples = load_data(patient_ID, patient_name)
        clones_per_tissue_pairs = pd.read_pickle(r"/path_to_sykes_children_data_folder/Clumpiness/{}_clones_per_tissue_pairs_new_1.pkl".format(patient_name))
        timepoints_per_sample, samples_per_timepoint = get_timepoint_sample_ids(patient_ID, patient_name)
        tissues_per_sample, samples_per_tissue = get_tissue_sample_ids(patient_ID, patient_name)
        pod_groups, _ = group_PODs(patient_ID, patient_name)
    
        results_1 = {}
        node_min = 3 #No reason to get clumpiness for clones with fewer than 3 nodes because their lineages are not complex enough.
        df_node_min = seq_table.loc[seq_table['clone_id'].isin(node_min_clone_dict[node_min])]
        df_node_min['pod'] = df_node_min['sample_id'].apply(lambda x: timepoints_per_sample[x])
        df_node_min['tissue'] = df_node_min['sample_id'].apply(lambda x: tissues_per_sample[x])
        
        for pod_group in pod_groups:
            print("Currently analyzing: " + str(pod_group) + " for " + subjects_dict[patient_ID])
            clones_in_pod_group = set(df_node_min['clone_id'].loc[df_node_min['pod'].isin(pod_group)].unique().astype(int)) 
            tissues_per_pod_group = df_node_min['tissue'].loc[(df_node_min['pod'].isin(pod_group))].unique()
            df_by_pod_and_node_min = df_node_min.loc[df_node_min['clone_id'].isin(clones_in_pod_group)]
            
                   
            results_1["POD" + str(pod_group)] = {}
            ###Make list of clones_ids to be analyzed
            for tissue1 in range(0, len(tissues) - 1):
                
                for tissue2 in range(tissue1 + 1, len(tissues)):
  
                    clones_in_tissue_pair = set(clones_per_tissue_pairs[tissues[tissue1] + "+" + tissues[tissue2]])           
                    clones_per_tissue_and_pod_group = clones_in_tissue_pair.intersection(clones_in_pod_group)
                    
                    ###Only add clones that have both tissues at this pod_group
                    list_of_clones_result = []
                    if tissues[tissue1] in tissues_per_pod_group and tissues[tissue2] in tissues_per_pod_group:
                        for clone in clones_per_tissue_and_pod_group:
                            tissues_per_clone_and_pod = df_by_pod_and_node_min['tissue'].loc[(df_by_pod_and_node_min['clone_id'] == clone)].unique()
                            if tissues[tissue1] in tissues_per_clone_and_pod and tissues[tissue2] in tissues_per_clone_and_pod:
                                list_of_clones_result.append(clone)

                    if len(list_of_clones_result) > 0:
                        print("Finished making list of clones. Now making files: " + str(len(list_of_clones_result)) + " total files.")
                        print("POD group: " + str(pod_group) + " Tissues: " + tissues[tissue1] + "+" + tissues[tissue2])
                        results_1["POD" + str(pod_group)][tissues[tissue1] + "+" + tissues[tissue2]] = {}
                        # =============================================================================
                        # Iterate through each tree, reformatting it and making files for find-clumpiness' "JSON" parser            
                        # =============================================================================    
                        df_of_trees_1 = df_of_trees.loc[df_of_trees['id'].isin(list_of_clones_result)]
            
                        for clone, clone_id in zip(df_of_trees_1['tree'], df_of_trees_1['id']):   
                            def parse(input):
                                result = []
                                for key, value in input.items():
                                    if key == "data":
                                        # test1 = value
                                        result.append(parse_data(value))
                                    if key == "children":
                                        if len(value) > 0:
                                            childrenResult = []
                                            for child in value:
                                                childrenResult.append(parse(child))
                                            result.append(childrenResult)
                                        else:
                                            result.append([])
                                return result
                            
                            def parse_data(data):
                                result = {"nodeID" : "ID", "nodeLabels" : []}
                                seq_ids = data['seq_ids']
                                for key, seq_id in seq_ids.items():
                                    metadata = seq_id['metadata']
                                    if len(metadata) > 0:
                                        metadata_pod = int(''.join(i for i in metadata['pod'] if i.isdigit()))
                                        ###Only count nodes that are in the current pod_group
                                        if metadata_pod in pod_group: 
                                            if 'Colon ' in metadata['sample_origin']:
                                                metadata['sample_origin'] = metadata['sample_origin'].replace("Colon ", "Colon")
                                            if 'Ileum' in metadata['sample_origin'] and metadata['pod'] == 'POD1145':
                                                metadata['sample_origin'] = metadata['sample_origin'].replace("Ileum", "Ileum_allograft")                                                
                                            result['nodeLabels'].append(metadata['sample_origin'])
                                return result
                            
                            json_data = json.loads(clone)
                            result = parse(json_data["tree"])                            
                            
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
                        results_1["POD" + str(pod_group)][tissues[tissue1] + "+" + tissues[tissue2]] = final_results
      
                    print("Now deleting find-clumpiness files to prepare for next iteration.")
                    # =============================================================================
                    # Delete all of the files that were created this iteration via the terminal
                    # =============================================================================
                    subprocess.run(r'find . -type f -delete', shell=True, check=True, stdout=subprocess.PIPE, universal_newlines=True)
                                
        # =============================================================================
        # Save results    
        # =============================================================================          
        path = Path(r"/path_to_sykes_children_data_folder/{}_clumpiness_POD_Specific.pkl".format(patient_name))
        pickle_out = open(path,"wb")
        pickle.dump(results_1, pickle_out)
        pickle_out.close() 
                 
    print("Finished!")
    
get_clumpiness_by_POD()


