# -*- coding: utf-8 -*-
"""
Created on Wed Jan 22 11:45:57 2020

@author: thsia

This script queries the mysql database and downloads the ImmuneDB tables and saves them locally as pickle files.
For the sequences, clones, clone_stats, and selection_pressure tables, the tables are downloaded for each subject individually.
The file name of the table will have the subject_id at the front. For example, subject_id 1 will have a clones_table named 1_clones_table.pkl.

Note: you must define the save paths and create the folders on your system before running the script.
Scroll down and update the paths to save the pickle files where you want them to be saved.

Also, this script connects to a remote machine via ssh. If you don't need to do that, you can remove the ssh part.

ImmuneDB was used to process and annotate the raw sequencing data. For more info visit: https://immunedb.readthedocs.io/en/latest/

"""


import json
from sqlalchemy import create_engine
import pandas as pd
from sshtunnel import SSHTunnelForwarder

 
# ssh variables
host = 'XXX.XX.XXX.XX' #IP of the machine hosting the database.
localhost = '127.0.0.1' #Your local machine's IP. Probably don't have to change this.
ssh_username = 'username' #Your ssh username
ssh_password = 'ssh_password'
port = 22 #ssh port. Probably won't have to change this but if you can't connect via ssh you may have to change this.


# mysql database variables
user='user' #username for the mysql database.
password='password' #password for the database.
database='database_name' #name of the database.

server = SSHTunnelForwarder(
    (host, port),
    ssh_username=ssh_username,
    ssh_password=ssh_password,
    remote_bind_address=(localhost, 3306) #This is the default port. You may have to change this if you have issues connecting.
    )

server.start()
local_port = str(server.local_bind_port)
connect_string = 'mysql+pymysql://{}:{}@{}:{}/{}'.format(user, password, localhost, local_port, database)
sql_engine = create_engine(connect_string)

# =============================================================================
# Download subjects for all subjects
# =============================================================================

print("Saving subjects...")
 
query = "select * from subjects"
result = pd.read_sql_query(query, sql_engine)
result.to_pickle(r"/path_to_lp16_data_folder/subjects.pkl")

def make_subject_key_dict(subjects_table):
    
    result = {}
    
    for subject_id, subject_name in zip(subjects_table['id'], subjects_table['identifier']):
        
        result[subject_id] = subject_name
        
    return result

subjects_dict = make_subject_key_dict(result)
     
# =============================================================================
# Download sequence_collapse_table for all subjects
# =============================================================================  
print("Saving sequence_collapse_table...") 
 
query = "select * from sequence_collapse"
result = pd.read_sql_query(query, sql_engine)
result.to_pickle(r"/path_to_lp16_data_folder/sequence_collapse.pkl")
     
# =============================================================================
# Download sample_metadata for all subjects
# =============================================================================
print("Saving sample_metadata...")
 
query = "select * from sample_metadata"
result = pd.read_sql_query(query, sql_engine)
result.to_pickle(r"/path_to_lp16_data_folder/sample_metadata.pkl")
      
# =============================================================================
# Download sample_stats for all subjects
# =============================================================================
print("Saving sample_stats...")
 
query = "select * from sample_stats"
result = pd.read_sql_query(query, sql_engine)
result.to_pickle(r"/path_to_lp16_data_folder/sample_stats.pkl")

# =============================================================================
# Download samples for all subjects
# =============================================================================
print("Saving samples...")
 
query = "select * from samples"
result = pd.read_sql_query(query, sql_engine)
result.to_pickle(r"/path_to_lp16_data_folder/samples.pkl")

# =============================================================================
# Download studies for all subjects
# =============================================================================
print("Saving studies...")
 
query = "select * from studies"
result = pd.read_sql_query(query, sql_engine)
result.to_pickle(r"/path_to_lp16_data_folder/studies.pkl")

# =============================================================================
# Download clones_table for all subjects
# =============================================================================
for subject_id in subjects_dict.keys():  
    print("Saving clones_table for id: ", subjects_dict[subject_id] + "...")
    
    query = "select * from clones where subject_id = {}".format(subject_id)
    result = pd.read_sql_query(query, sql_engine)
    result.to_pickle(r"/path_to_lp16_data_folder/clones_table/{}_clones_table.pkl".format(subject_id))

# =============================================================================
# Download sequences_table for all subjects
# =============================================================================

for subject_id in subjects_dict.keys():
    print("Saving sequences_table for id: ", subjects_dict[subject_id] + "...")
        
    query = "select * from sequences where subject_id = {} and clone_id IS NOT NULL".format(subject_id)
    result = pd.read_sql_query(query, sql_engine)
    result.to_pickle(r"/path_to_lp16_data_folder/sequences_table/{}_sequences_table.pkl".format(subject_id))

# =============================================================================
# Download clone_stats for all subjects
# =============================================================================

for subject_id in subjects_dict.keys():
    print("Saving clone_stats for id: ", subjects_dict[subject_id] + "...")
        
    query = "select * from clone_stats where subject_id = {} and id IS NOT NULL".format(subject_id)
    result = pd.read_sql_query(query, sql_engine)
    result.to_pickle(r"/path_to_lp16_data_folder/clone_stats/{}_clone_stats.pkl".format(subject_id))
     
# =============================================================================
# Download selection_pressure for all subjects
# =============================================================================
for subject_id in subjects_dict.keys():
    print("Saving selection_pressure for id: ", subjects_dict[subject_id] + "...")
    
    ###selection_pressure table does not have subject_id column so must filter by clone_id for each subject
    subject_clone_stats = pd.read_pickle(r"/path_to_lp16_data_folder/clone_stats/{}_clone_stats.pkl".format(subject_id))
    subject_clones = subject_clone_stats['clone_id'].unique().tolist()
    subject_clones = tuple(subject_clones) 
    placeholders = ", ".join(["%s"] * len(subject_clones))

    query = "select * from selection_pressure where clone_id in ({})".format(placeholders)
    
    result = pd.read_sql_query(query, sql_engine, params=subject_clones)
    result.to_pickle(r"/path_to_lp16_data_folder/selection_pressure/{}_selection_pressure.pkl".format(subject_id))
    print(result.shape)

print("Finished!")
