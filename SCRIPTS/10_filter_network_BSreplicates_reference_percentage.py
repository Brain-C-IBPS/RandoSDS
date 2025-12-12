 #!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""

            This script is part of the Rando-SDS v1.0 suite

------------------------------------------------------------------------------

Python script to filter SDS network files in order to keep  edges that resist
a certain number of bootstrap replicate analysis, i.e. to keep edges that appear
e.g. on 80% of the bootstrap replicate networks. 

The SDS network of choice is filtered based on a CSV file produced from running
the Python script 9_calculate_BSreplicates_edge_presence.py.

Python script to modify gene names in input BS replicates LFC files from 
ENSEMBL IDs to gene symbols.

------------------------------------------------------------------------------

Contributors (alphabetic order): Maialen Arrieta-Lobo (1) 
                                 Lucile Mégret
                                 Christian Neri (1)
    
Laboratories:   
(1) Institut de Biologie Paris Seine (IBPS), UMR CNRS 8256, 
    Team Brain-C (Compensation systems in neurodegenerative diseases and aging)

Affiliations:
(1) Sorbonnes Université
(1) Centre National de la Recherche Scientifique (CNRS)
(1) Institut National de la Santé et de la Recherche Médicale (INSERM)

==============================================================================
|  Copyright (C) 2023 Maialen Arrieta-Lobo, Lucile Mégret and Christian Neri            | 
==============================================================================
|                                                                            |
|  This work is licensed under the Creative Commons                          | 
|  Attribution-NonCommercial-NoDerivatives 4.0                               |
|  International License. To view a copy of this license, visit              |
|  http://creativecommons.org/licenses/by-nc-nd/4.0/ or send a letter to     |
|  Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.              | 
|                                                                            |                                                                          |
============================================================================== 


Input:
-----
    + SDS network to be filtered, in SIF format
    + CSV file outputted from $calculate_BSreplocates_edge_presence.py

The CSV file for filtering contains two columns: 
      
    + edge_in_networks_counter: integer column representing the number of 
                                networks the edge was found in
    + edge_in_networks_percentage: float column containing the percentage of 
                                   networks the edge was found in
                                   
                                   
Output:
------
    + SIF file containing the edges that were present in the percentage threshold 
      set by the user through the $PERCENTAGE_THRESHOLD variable


"""


# ---------------
# Import modules
# ---------------
import pandas as pd
import numpy as np
import networkx as nx
import os
from datetime import datetime
import argparse



# ---------------
#      MAIN
# ---------------
time_start = datetime.now()


# Initialise the parser and read arguments:
# -----------------------------------------
parser = argparse.ArgumentParser('Script to filter SDS network')

required_named = parser.add_argument_group('Required arguments')

# Add long and short arguments
required_named.add_argument("--edge_presence_threshold", "-e", 
                    type = float, 
                    help = ("State edge presence filtering threshold, float between" 
                    "0 and 100 (e.g. 75.5 for edges that resist 75.5% of bootstrapping)"),
                    required = True)

required_named.add_argument("--sds_deregulation", "-deregulation", 
                    type = str,
                    choices =['ABS', 'U', 'D'],
                    help = "Select SDS deregulated network",
                    required = True)

args = parser.parse_args()   

# Print out inputted arguments:
# -----------------------------
print("Edge presence threshold for network filtering:  %s" % args.edge_presence_threshold)
PERCENTAGE_THRESHOLD = args.edge_presence_threshold

print("Deregulated networks to be filtered:  %s" % args.sds_deregulation)
SDS_DEREGULATION = args.sds_deregulation 


# Create output folder if non-existent:
# -------------------------------------
output_folder = './ANALYSIS_RESULTS/FILTERED_NETWORKS/'
CHECK_FOLDER = os.path.isdir(output_folder)

if not CHECK_FOLDER:
    os.makedirs(output_folder)

# Iterate over disease condition
# -------------------------------
disease_condition_array = ['Q50', 'Q111', 'Q170', 'Q175'] # In the example file, the original experiment was performed for 4 mutation lengths


for disease_condition in disease_condition_array:


    print(' -----------------------------------')
    print('       Filtering {} {} network' .format(disease_condition, SDS_DEREGULATION))  
    

    # Read file containing the edge presence counter and percentage
    # -----------------------------------------------------------------
    reference_edge_file = ('./ANALYSIS_RESULTS/BS_REPLICATES_EDGE_PRESENCE/'
                           'BSreplicates_edge_presence_scores_{}{}.csv'
                           .format(disease_condition, SDS_DEREGULATION))
    
    # Convert to df
    reference_edge_df = pd.read_csv(reference_edge_file, header = 0)
    
    # Filter reference dataframe to retain edges above the selected percentage threshold
    reference_edges_above_threshold_df = reference_edge_df[reference_edge_df.edge_in_networks_percentage >= PERCENTAGE_THRESHOLD]
    reference_edges_above_threshold_df = reference_edges_above_threshold_df[['First', 'Last']]
    
    # Sort dataframe
    reference_edges_above_threshold_df_sorted = reference_edges_above_threshold_df.sort_values(by = ['First', 'Last'], 
                                                                                                inplace = False, 
                                                                                                ascending = True)
    
    # Remove duplicates if any
    reference_edges_above_threshold_df = reference_edges_above_threshold_df_sorted[~reference_edges_above_threshold_df_sorted.duplicated()]
    
    # Delete dataframes to release memory space
    del reference_edge_df, reference_edges_above_threshold_df_sorted
    
    
    # Read current network SIF file
    # -----------------------------
    file_to_filter = ('./ANALYSIS_RESULTS/SDS_NETWORKS/randoSDS_original_experiment/'
                      '{}_LFC_geneSymbols/{}/{}_10.sif'
                      .format(disease_condition, SDS_DEREGULATION, SDS_DEREGULATION))
    
    # Convert to df
    df_to_filter = pd.read_csv(file_to_filter, sep = ' ', header = None) 
    df_to_filter = df_to_filter.drop(1, axis = 1)
    df_to_filter.columns = ['First', 'Last']
    df_to_filter["First"] = df_to_filter["First"].astype(str)
    df_to_filter["Last"] = df_to_filter["Last"].astype(str)
    
    # Remove duplicates if any
    df_to_filter = df_to_filter[~pd.DataFrame(np.sort(df_to_filter.values), 
                                                     columns = df_to_filter.columns, 
                                                     index = df_to_filter.index).duplicated(keep = 'first')]
    
    
    # Use NetworkX lnbrary for edge filtering
    # ---------------------------------------
    
    # Initialise array to store edges from reference edge presence network
    reference_edges_list = []    
    
    # Iterate over rows in filtered reference edge df and store to array
    for index, rows in reference_edges_above_threshold_df.iterrows(): 
        my_list = [rows.First, rows.Last]          
        reference_edges_list.append(my_list) 
        
    # Generate networkX object from edge list
    reference_net = nx.from_edgelist(reference_edges_list)
    reference_edges = set(reference_net.edges())
    
    # Create networkX object from network to be filtered
    edges_to_check_list = []    
    
    for index, rows in df_to_filter.iterrows(): 
        my_list = [rows.First, rows.Last]          
        edges_to_check_list.append(my_list) 
    
    net_to_filter = nx.from_edgelist(edges_to_check_list)
    edges_to_filter = set(net_to_filter.edges())
    
    # Iterate over edges to filter:
    edges_to_keep = set()
    
    for v in net_to_filter.edges():
        
        if (v[1], v[0]) in edges_to_keep or (v[0], v[1]) in edges_to_keep: 
            continue
        if ((v[0], v[1]) in reference_net.edges()) or ((v[1], v[0]) in reference_net.edges()):
            edges_to_keep.add(v)
       
          
    # Convert to dataframe and add 'IN' interaction column for each edge
    edges_to_keep_df = pd.DataFrame({'edge' : list(edges_to_keep)})
    edges_to_keep_df[['g1', 'g2']] = pd.DataFrame(edges_to_keep_df['edge'].tolist())
    edges_to_keep_df = edges_to_keep_df.drop('edge', axis = 1)
    edges_to_keep_df = edges_to_keep_df[~pd.DataFrame(np.sort(edges_to_keep_df.values), 
                                                     columns = edges_to_keep_df.columns, 
                                                     index = edges_to_keep_df.index).duplicated(keep = 'first')]
    
    # Are the two methods yielding the same filtered network??  
    tmp = ['IN'] * len(edges_to_keep_df)
    edges_to_keep_df.insert(1, 'Interaction', tmp)  
    
    
     
    # Filtered network output filename
    filename = ('filtered_network_{0}percentageBSreps_{1}{2}.sif'
                .format(PERCENTAGE_THRESHOLD, disease_condition, SDS_DEREGULATION))
    
    edges_to_keep_df.to_csv(os.path.join(output_folder, filename), 
                           header = False, index = False, sep = '\t')    
    
# -----------------------------------------------------------------------------

time_finish = datetime.now()

print('''
\t==================================================================== 

\t              Script ended after {0} 
      
\t====================================================================         

'''.format( time_finish - time_start))
