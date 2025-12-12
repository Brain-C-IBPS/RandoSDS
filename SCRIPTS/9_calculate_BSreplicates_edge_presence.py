#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

            This script is part of the Rando-SDS v1.0 suite

------------------------------------------------------------------------------

Python script to quantify the presence of a given SDS network edge across 
the 'original' network and/or several bootstrap replicates of the network. 

Python script to modify gene names in input BS replicates LFC files from 
ENSEMBL IDs to gene symbols

------------------------------------------------------------------------------

Contributors (alphabetic order): Maialen Arrieta-Lobo (1) 
                                 Lucile Mégret (1)
                                 Christian Neri (1)
    
Laboratories:   
(1) Institut de Biologie Paris Seine (IBPS), UMR CNRS 8256, 
    Team Brain-C (Compensation systems in neurodegenerative diseases and aging)

Affiliations:
(1) Sorbonnes Université
(1) Centre National de la Recherche Scientifique (CNRS)
(1) Institut National de la Santé et de la Recherche Médicale (INSERM)

==============================================================================
| Copyright (C) 2023 Maialen Arrieta-Lobo, Lucile Mégret and Christian Neri  | 
==============================================================================
|                                                                            |
|  This work is licensed under the Creative Commons                          | 
|  Attribution-NonCommercial-NoDerivatives 4.0                               |
|  International License. To view a copy of this license, visit              |
|  http://creativecommons.org/licenses/by-nc-nd/4.0/ or send a letter to     |
|  Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.              |
|                                                                            |
============================================================================== 


Input:
------
    + Path to SDS network SIF files (original and/or boostrap replicates)

Output:
-------
    + CSV file containing each unique edge in the network combination per row
    + CSV file format: 2-columns
        => First column contains the count of the number of networks the edge 
           was present in
        => Second column contains the corresponding percentage for the total 
           number of networks taken into account

Requires previously running SDS analyses on the original and/or boostrapped LFC
data to produce SIF files. 

The script uses the networkX module to convert the SIF file into a network and
consider the edges as a tuple. 

The output file generated with this script can be used to filter out edges in 
the original network based on a presence score threshold set by the user by
running the script '10_filter_network_BSreplicates_reference_percentage.py'


"""

# -------------------------------------
#       Import modules
# -------------------------------------
import pandas as pd
import numpy as np
import networkx as nx
from datetime import datetime
import os
import argparse



# ---------------
#      MAIN
# ---------------
time_start = datetime.now()


# Initiate the parser and read arguments:
# ---------------------------------------
parser = argparse.ArgumentParser('Script to calculate in how many BS replicates an edge occurs')

required_named = parser.add_argument_group('Required arguments')

# Add long and short argument
required_named.add_argument("--n_bs_replicates", 
                            "-n2", 
                            type = int,
                            help = "Set number of bootstrap replicates generated in step 1",
                            required = True)

# Add long and short argument
required_named.add_argument("--deregulation", 
                            "-d", 
                            type = str,
                            choices =['ABS', 'U', 'D'],
                            help = "Select deregulation type ",
                            required = True)
args = parser.parse_args()

# Check for --n_bs_replicates
print("Number of bootstrap replicates generated for the experiment = %s" % args.n_bs_replicates)
N_BS_REPLICATES = int(args.n_bs_replicates)

# Check for --deregulation
print("Type of deregulated network to analyse = %s" % args.deregulation)
DEREGULATION = str(args.deregulation)

# --------------------------------------------------
# Initialize dictionary to store edge presence count
# --------------------------------------------------
edge_dictionary = {}

effect = DEREGULATION
disease_conditions = ['Q50', 'Q111', 'Q170', 'Q175']

output_folder = './ANALYSIS_RESULTS/BS_REPLICATES_EDGE_PRESENCE/'
CHECK_FOLDER = os.path.isdir(output_folder)

if not CHECK_FOLDER:
    os.makedirs(output_folder)

# ----------------------------------------------------------------------
#  Get SDS networks for original experiment and add edges to dictionary 
# ----------------------------------------------------------------------
for condition in disease_conditions:
    
    original_data_network_file = ('./ANALYSIS_RESULTS/SDS_NETWORKS/'
                                  'randoSDS_original_experiment/{}_LFC_geneSymbols/'
                                  '{}/{}_10.sif'.format(condition, effect, effect))
    
    original_net_df = pd.read_csv(original_data_network_file, sep = ' ', header = None)    
    original_net_df = original_net_df.drop(1, axis = 1)
    original_net_df.columns = ['First', 'Last']
        
    # Sort columns
    df = pd.DataFrame(np.sort(original_net_df.values, axis = 1), 
                        columns = original_net_df.columns)
     
    # Remove duplicate edges if any
    original_unique_edges_df = df[~df.duplicated()]
    
    # Create an empty list to store edges as list
    original_net_row_list = [] 
    
    # Iterate over each row 
    for index, rows in original_unique_edges_df.iterrows(): 
        # Create list for the current row 
        my_list = [rows.First, rows.Last] 
          
        # append the list to the final list 
        original_net_row_list.append(my_list) 
    
    # Create networkX graph from list of edges
    original_net_df = nx.from_edgelist(original_net_row_list)
    
    # Store edges as list
    net_edges = list(original_net_df.edges())
    
    # Iterate over edges to add them to dictionary
    for edge in net_edges:
        
        # Store edge in edge dictionary for eventual edge filtering:
        if edge not in edge_dictionary.keys():
            
            edge_dictionary[edge] = 1
            
        else: 
            
            edge_dictionary[edge] += 1 
        
        
    # ----------------------------------------
    #  Repeat process for bootstrap replicates
    # ----------------------------------------
        
    for i in range(N_BS_REPLICATES):
            
        print(' -------------------------------------------')
        print('  Analysing BS network replicate {} of {}'.format(i + 1, condition))  
        
        path = ('./ANALYSIS_RESULTS/SDS_NETWORKS/randoSDS_BS{}/{}_LFC_geneSymbols/'
                '{}/{}_10.sif'.format( i + 1, condition, effect, effect))
        
        net_df = pd.read_csv(path, sep = ' ', header = None)    
        net_df = net_df.drop(1, axis = 1)
        net_df.columns = ['First', 'Last']
        
        df1 = pd.DataFrame(np.sort(net_df.values, axis = 1), 
                           columns = net_df.columns)
        
        unique_edges_df = df1[~df1.duplicated()]
     
        net_row_list = []    
    
        for index, rows in unique_edges_df.iterrows(): 
            my_list = [rows.First, rows.Last]          
            net_row_list.append(my_list) 
            
        net = nx.from_edgelist(net_row_list)
        
        net_edges = list(net.edges())
        
        for edge in net_edges:
            
            if edge not in edge_dictionary.keys():
                
                edge_dictionary[edge] = 1
                
            else: 
                
                edge_dictionary[edge] += 1
        
    
    
    # -----------------------------------------------------
    #  Convert edge presence count dictionary to a DataFrame
    # -----------------------------------------------------
    edge_df = pd.DataFrame.from_dict(data = edge_dictionary,
                                     orient = 'index', 
                                     columns = ['edge_in_networks_counter'])
    
    edge_df['Edge'] = edge_df.index
    
    edge_df[['First', 'Last']] = pd.DataFrame(edge_df['Edge'].tolist(), index = edge_df.index)
    
    # ------------------------------------------------------------------------
    # Add column reporting the percentage of networks each edge is present in
    # ------------------------------------------------------------------------
    total_networks = N_BS_REPLICATES
    edge_df['edge_in_networks_percentage'] = edge_df['edge_in_networks_counter'] * 100. / total_networks    
    
    # ---------------------------
    # Save dataframe to CSV file
    # ---------------------------
    edge_df.to_csv(os.path.join(output_folder, 'BSreplicates_edge_presence_scores_{}{}.csv'.format(condition, effect)), 
                   header = True)    
    
    
    
    
    # Create SIF file too
    
    sif_df = edge_df[['First', 'Last']]
    
    tmp = ['IN'] * len(sif_df)
    sif_df.insert(1, 'Interaction', tmp)    
    sif_df.to_csv(os.path.join(output_folder, 'BSreplicates_edge_presence_scores_{}{}.sif'.format(condition, effect)), 
                   header = False, index = False, sep = '\t')    
    
time_finish = datetime.now()


print('''
\t==================================================================== 

\t              Script ended after {0} 
      
\t====================================================================         

'''.format( time_finish - time_start))





















