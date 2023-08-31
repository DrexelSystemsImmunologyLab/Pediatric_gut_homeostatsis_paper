To run these scripts, you must first run the "preprocessing" scripts. Some of the preprocessing scripts are prerequisites for others, so run them in the order they are numbered to avoid running into issues.

The first preprocessing script, 1_Download_tables_as_pickle.py, queries an ImmuneDB mysql database and saves the tables locally as pickle files. Most tables are saved individually for each subject in the database.

Notes: 

1. You will need to change the file paths in each script to point to whereever you want to save the data.

2. The dictionaries at the top of each script that show the subjects will need to be replaced with the subjects from your database (unless you used the exact database used in this paper). If your database is a remade version, the subject_ids will most likely be different from what is shown in the scripts. Therefore, you will need to get {subject_id : subject_name} for each subject in your database and then replace the subjects_dict at the top of the script with your subjects. You can easily find this information for your database in the 1_Download_tables_as_pickle.py script, where a dictionary of {subject_id : subject_name} is made, containing each subject in the database.

3. 12_get_clumpiness.py and 6_get_diveristy_of_clone_size_instances_by_POD_brackets.py require the corresponding Haskell program to be installed on your machine. More information is provided at the top of those scripts.

4. These scripts require python version 3+ and Plotly version 5+.

