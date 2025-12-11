#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""

            This script is part of the Rando-SDS v1.0 suite

------------------------------------------------------------------------------

Python script to modify gene names in input LFC files from ENSEMBL IDs 
to gene symbols

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

INPUT
-----

    LFC files corresponding to the original experiment (not the BS replicates of the 
    experiment) containing ENSEMBL gene IDs in first column and LFC values on second colum.
    
    By default, input LFC files are taken from the /DESEq2_RESULTS/ folder within the 
    /ANALYSIS_RESULTS/PRE_PROCESSING_ORIGINAL_EXPERIMENT/ directory, and have a naming
    structure that follows 'lfc_.txt'.
    
    The user should modify either the path to LFC files and the naming structure thereof
    if not concurrent with the default settings. 

OUTPUT
------

    'lfc_*_geneSymbols.txt' files inside the input data folder containing gene symbols
    in the first column and corresponding LFC values in the second column. 


"""

# -------------------------------------
#       Import modules
# -------------------------------------

import pandas as pd
import mygene
import os
from datetime import datetime
import glob
from pathlib import Path

# ---------------
#      MAIN
# ---------------

time_start = datetime.now()

#  ----------------------------------------------------------------------------   

# Set path and structure of LFC files to modify:
lfc_files_path_string = './ANALYSIS_RESULTS/PRE_PROCESSING_ORIGINAL_EXPERIMENT/DESEq2_RESULTS/lfc_*.txt'
    
# Retrieve LFC files for which we want to change ENSEMBL IDs to GeneSymbols:
lfc_files_to_modify = []

# Loop over files and query mygene library for translation:
for lfc_file in glob.glob(lfc_files_path_string):   
    
    # Read as dataframe and modify column names
    lfc_df = pd.read_csv( lfc_file, header = None, sep = '\t')
    lfc_df.columns = ['gene', 'lfc']
    
    # Get gene list:
    genes = lfc_df.loc[:, 'gene'].tolist()
    
    # Initialise dictionary to store translated names
    gene_dict = {}
    
    # Launch mygene instance
    mg = mygene.MyGeneInfo()
    
    # Make a query to convert gene list names:    
    results = mg.querymany(genes, 
                           scopes = ["ensembl.gene"], 
                           fields = ["symbol"], 
                           species = "mouse", 
                           verbose = False)
    
    # Start looping over query results:
    not_found_array = []
    i = 0
    for res in results:
        q = res['query']
        s = 'NaN'
        if 'symbol' in res:
            s = res['symbol']
        else:
            not_found_array.append(i)
        i += 1 
        gene_dict[q] = s
    
    # Store conversion to dataframe    
    new_lfc_df = lfc_df.replace({"gene": gene_dict})    
    new_lfc_df = new_lfc_df[ ~new_lfc_df.index.isin( not_found_array )]
    new_lfc_df = new_lfc_df.sort_values(by = ['gene'])
    
    # Filter out genes with 0 counts
    new_lfc_df = new_lfc_df[ new_lfc_df['lfc'] != 0.0]
    
    # Save df with gene symbols instead of ensembl ids to txt file
    outdir = os.path.dirname(lfc_file)
    old_file_name = Path(lfc_file).stem
    new_file_name = ''.join([old_file_name, '_geneSymbols.txt'])
    new_lfc_df.to_csv( os.path.join(outdir, new_file_name), 
                               sep = '\t', index = False, header = False)
        

#  ----------------------------------------------------------------------------   

time_finish = datetime.now()

print('''
\t==================================================================== 

\t              Script ended after {0} 
      
\t====================================================================         

'''.format( time_finish - time_start))


