#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""

            This script is part of the Rando-SDS v1.0 suite

------------------------------------------------------------------------------

Python script to generate bootstrap replicates of experimental read count data 
that must more than one replicate. 

------------------------------------------------------------------------------

Contributors (alphabetic order): Maialen Arrieta-Lobo (1) 
                                 Christian Neri (1)
    
Laboratories:   
(1) Institut de Biologie Paris Seine (IBPS), UMR CNRS 8256, 
    Team Brain-C (Compensation systems in neurodegenerative diseases and aging)

Affiliations:
(1) Sorbonnes Université
(1) Centre National de la Recherche Scientifique (CNRS)
(1) Institut National de la Santé et de la Recherche Médicale (INSERM)

==============================================================================
|      Copyright (C) 2023 Maialen Arrieta-Lobo and Christian Neri            | 
==============================================================================
|                                                                            |
|  This work is licensed under the Creative Commons                          | 
|  Attribution-NonCommercial-NoDerivatives 4.0                               |
|  International License. To view a copy of this license, visit              |
|  http://creativecommons.org/licenses/by-nc-nd/4.0/ or send a letter to     |
|  Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.              |
|                                                                            |
==============================================================================      

------------------------------------------------------------------------------

INPUT
-----

    File containing at least the following columns
    
    1 - gene name column
    2 - one column per replicate and disease condition, column names having the 
    following format: ${CELLTYPE}_${POLYQLENGTH}_${REPLICATENUMBER}_counts
    
    By default, it runs on the example count data file provided with this notebook where:
        
    ${CELLTYPE} = CP73
    ${POLYQLENGTH} = Q20, Q50, Q111, Q170 or Q175
    ${REPLICATENUMBER} = 1, 2, 3, 4, 5, 6, 7, 8, 9 or 10
    
    We want to resample those mouse replicates to reproduce Boostrap replicates of the experiment
    from which we can calculate the corresponding LFC files and apply SDS to them
    
    To adapt it to another input count file, you will need to change the line of code 105

Steps of script:
----------------
    
    1 - Resample the N_EXPERIMENTAL_REPLICATES replicates of the experiment with replacement to 
        generate N_BS_REPLICATES bootstrap replicates of the experiment
        
    2 - For each BS replicate of the experiment:
        
        2.1 - Extract read data columns corresponding to the replicates within the BS replicate
        
        2.2 - Generate a new Data Frame (BS data frame) containing those columns 
    
        2.3 - Store the BS data frame as a CSV file in the Bootstrap_replicates folder
        
        
OUTPUT:
-------

    Folder containing:
    
        1 - bootstrap_replicates_ids.txt  ===> file containing the list of replicates
                                                used in each BS replicate
                                                The file contains N_BS_REPLICATES rows, 
                                                and each row features N_EXPERIMENTAL_REPLICATES
                                                elements.
                                                
        2 - N_BS_REPLICATES CSV files containing the replicates of the experiments                                           


"""


# -------------------------------------
#       Import modules
# -------------------------------------
import numpy as np
import pandas as pd
import os
import collections
import argparse

# -------------------------------------
#       Function definition
# -------------------------------------
def draw_bs_sample(data):
    """Draw a bootstrap sample from a 1D data set."""
    return np.random.choice(data, size = len(data))

# -------------------------------------
#             Main
# -------------------------------------

# Create directory to store bootstrap replicates of the experimental data
# -----------------------------------------------------------------------
outdir = './ANALYSIS_RESULTS/BOOTSTRAP_REPLICATES_OF_EXPERIMENT/'

CHECK_FOLDER = os.path.isdir(outdir)
if not CHECK_FOLDER:
    os.makedirs(outdir)


# Initialise parser and read arguments
# -------------------------------------
parser = argparse.ArgumentParser()

# Add long and short argument
parser.add_argument("--n_experimental_replicates", "-n1", type = int, help = "State number of replicates in experiment")
parser.add_argument("--n_bs_replicates", "-n2", type = int,help = "Set number of bootstrap replicates to generate")
args = parser.parse_args()

# Check for --n_experimental_replicates
if args.n_experimental_replicates:
    print("Number of experimental replicates = %s" % args.n_experimental_replicates)
N_EXPERIMENTAL_REPLICATES = int(args.n_experimental_replicates)  # Number of replicates in actual experiment

# Check for --n_bs_replicates
if args.n_bs_replicates:
    print("Number of bootstrap replicates of the experiment to generate = %s" % args.n_bs_replicates)
N_BS_REPLICATES = int(args.n_bs_replicates)


# Create array to identify each replicate to be generated:
# -------------------------------------------------------
replicate_ids = [ i + 1 for i in range(N_EXPERIMENTAL_REPLICATES) ]

# Create dictionary containing replicate ID as string
replicate_ids_dict = {}

for identifier in replicate_ids:
    
    replicate_ids_dict[identifier] = "_{}_".format(identifier)
    
    
# Read experimental data file
# ---------------------------
read_count_file = 'EXAMPLE_INPUT_DATA/example_read_count_file.csv'
read_count_df = pd.read_csv(read_count_file, sep = '\t', header = 0)

#  Get cfile olumns containing generic information other than read counts
base_cols = read_count_df.columns.tolist()[:2]

    
# Generate boostrap replicates of the identifiers
# -----------------------------------------------
# Initialise array for further use:
bs_replicates = []

# Create file to keep track of BS replicate IDs
bs_replicate_ids_file = open(os.path.join(outdir, 'bootstrap_replicates_ids.txt'), 'a')

# Iterate over BS replicate numbers and create replicates:
for j in range(N_BS_REPLICATES):
    
    # Draw a BS replicate of the mice in experiment
    replicate = draw_bs_sample(replicate_ids)
    
    # Save the BS replicate in an array
    bs_replicates.append(replicate)
    
    # Save mouse IDs in BS replicate to keep track
    np.savetxt(bs_replicate_ids_file, [replicate],  fmt ='%s')

# Close file
bs_replicate_ids_file.close()


# For each replicate, generate a CSV file corresponding to the replicates of the bootstrap
# -----------------------------------------------------------------------------------------
i = 1

for bs_replicate in bs_replicates:
    
    print(' ---------------------------------------')
    print('    Generating BS replicate number ', i)  
    print(' ---------------------------------------')
    
    # Initialize array that will contain the df columns containing the data of 
    # the replicate in current bs replicate
    bs_replicate_df_cols = []
    
    # Get the columns for each mice present in the bootstrapped experiment  
    replicate_str_array = []
    
    for replicate in bs_replicate:
        
        replicate_str = replicate_ids_dict[replicate]
        replicate_str_array.append(replicate_str)
        
        # Get columns for replicate in question
        cols = [ c for c in read_count_df.columns if replicate_str in c ]
        
        bs_replicate_df_cols.append(cols)
        
    # Add generic columns to df columns
    bs_replicate_df_cols = sum(bs_replicate_df_cols, [])
    bs_replicate_df_cols = base_cols + bs_replicate_df_cols
    bs_replicate_df =  read_count_df.loc[:, bs_replicate_df_cols] 
    
    # Identify replicate IDs duplicated in the bs replicates to avoid naming conflicts
    duplicates = [item for item, count in collections.Counter(replicate_str_array).items() if count > 1]
   
    # Loop to change column names if column name contains a duplicated replicate    
    new_cols = []    
    dict_dup = {}
    
    for dup in duplicates:
        
        dict_dup[dup] = 0
    
    for col in bs_replicate_df.columns:
        
        # Column contains duplicate:
        if any(dup in col for dup in duplicates):
            
            # Find what duplicate            
            for dup in duplicates:
                         
                if dup in col:
                    
                    col = col.replace( dup, ''.join([dup, 'AA', str(dict_dup[dup]), 'AA_']))
                    dict_dup[dup] += 1 
                
        new_cols.append(col)
        
    bs_replicate_df.columns = new_cols
    
    # Save replicate_df to CSV file
    bs_replicate_df.to_csv( os.path.join(outdir, 'BSreplicate_{}.csv'.format(i)), 
                           sep = '\t', index = False)
    
    i += 1
    del bs_replicate_df








