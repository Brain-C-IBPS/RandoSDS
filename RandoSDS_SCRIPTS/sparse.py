#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""

             This module is part of the Rando-SDS v1.0 suite

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
| Copyright (C) 2025 Maialen Arrieta-Lobo, Lucile Mégret and Christian Neri  | 
==============================================================================
|                                                                            |
|  This work is licensed under the Creative Commons                          | 
|  Attribution-NonCommercial-NoDerivatives 4.0                               |
|  International License. To view a copy of this license, visit              |
|  http://creativecommons.org/licenses/by-nc-nd/4.0/ or send a letter to     |
|  Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.              | 
|                                                                            |                                                                          |
============================================================================== 

"""

# ----------------------------------------------------------------------------
#
#                           Import modules
#
# ----------------------------------------------------------------------------

import glob, os, re, pandas as pd
import numpy as np
import time
import datetime

# ----------------------------------------------------------------------------
#
#                           Function definition
#
# ----------------------------------------------------------------------------



def get_sparse_graph_edges( output_folder_path ):
    
    """
    
    Function that fetches and processes the network files to convert them
    into pandas DataFrames
    
    """
    
    # Select first of output folders, since the graph is the same for all of them, 
    # and thus we cover for the case were there is a single output file path
    sparse_graph_file = os.path.join( output_folder_path, 'graph.sparse' )
     
    sparse_graph_df = pd.read_csv( sparse_graph_file, 
                                  sep = "\t", header = 0,
                                  dtype = {'start_gene_idx' : int,
                                           'end_gene_idx' : int,
                                           'edge_weights' : float} )   
    
    edge_gene_matrix = sparse_graph_df.values
    
    return sparse_graph_df, len(sparse_graph_df)
 
# -----------------------------------------------------------------------------
    
def process_network_files( path_to_bionetwork_folder, edge_score_thresh = None):
    
    """
    
    Function that fetches and processes the network files to convert them
    into pandas DataFrames
    
    """
    
    # Read genes and edges files from bionetwork folder path:
    path_to_network_genes_file =  os.path.join( path_to_bionetwork_folder,  'geneNames.txt')   
    path_to_network_edges_file = os.path.join( path_to_bionetwork_folder,  'networkNames.txt')    
    
    # Convert network gene name file into pandas dataframes and read genes:
    graph_genes_df = pd.read_csv( path_to_network_genes_file, sep = "\t", header = None, index_col = False) 
    graph_genes = graph_genes_df[0].values # Convert form df column to list
     
    # Convert network edge gene names file into pandas dataframes and read genes in edges:
    graph_edges_df = pd.read_csv( path_to_network_edges_file, sep = "\t", header = None, index_col = False)  
    
    
     # Check if score threhold was selected for edge filtering:
    if edge_score_thresh is not None:
        print('\n'
               '\t           ****** WARNING ****** \n'
               '\n Only edges with score >= {0} will be kept in the network \n'
               '\n'
               '\t           ********************* \n'
               '\n'
               .format(edge_score_thresh))
               
        graph_edges_df = graph_edges_df[graph_edges_df[2].values >= edge_score_thresh ]
        
          
    return graph_genes, graph_edges_df #graph_edges_startGenes, graph_edges_endGenes, graph_edge_weights

# ----------------------------------------------------------------------------    

def process_lfc_files(  lfc_file_path ):
    
    """
    
    Function that fetches and processes the LFC files to convert them
    into pandas DataFrames
    
    """
    
    # Get the LFC data from the input file and the expressed genes from its 1st column
    expression_df = pd.read_csv( lfc_file_path, sep = '\t', header = None)
    expression_df.columns = ['gene', 'enrichment']
    expressed_genes = expression_df['gene'].values
    
    return expression_df, expressed_genes

# ----------------------------------------------------------------------------

def check_directory_exists( directory_name ):
    
    '''Check if directory exists, if not, create it'''
   
    MYDIR = directory_name
    CHECK_FOLDER = os.path.isdir(MYDIR)
    
    # If folder doesn't exist, create it.
    if not CHECK_FOLDER:
        os.makedirs(MYDIR)
    
    
        
# ----------------------------------------------------------------------------    

def create_output_directory( lfc_file_path, path_to_input_folder, output_prefix):
    
    """
    
    Funtion to create output directory based on the LFC file name
    
    """

    # Get lfc file name from lfc file path:
    lfc_file_name = os.path.basename(lfc_file_path)
        
    # Set output directory name and path
    output_directory_name = lfc_file_name.replace( ".txt", "").replace( "lfc", output_prefix)    
    output_directory_path =  os.path.join(path_to_input_folder, output_prefix)
    output_directory_path = os.path.join(output_directory_path, output_directory_name)
    
    # Call function to check whether output directory exists
    check_directory_exists(output_directory_path)
    
    return output_directory_path




# ----------------------------------------------------------------------------


def matched_index_finder_dict( l1, l2 ):
    
    """
    
    Function that finds the indices at which the elements of l2 occur in l1
    
    Input: two arrays

    Output: list containing the l1 indices at which the elements in l2 are located
   
    
    """
    
    dict1 =  { k: v for v, k in enumerate(l1)}
    
    match = [ dict1[k] for k in l2 if k in dict1.keys()]
   
    return match


# ----------------------------------------------------------------------------


def missing_gene_finder_dict( l1, l2 ):
    
    """
    
    Function that finds the elements in data2 missing in data1
    
    Input: two arrays
    Output: list containing the indices in data1 at which the elements in data2 
            are located
            
    TODO: OPTIMIZE THIS METHOD! Takes a very long time
    
    """
    miss = [ k for k in l2 if k not in l1]
   
    return miss

# ----------------------------------------------------------------------------

def sparse_df_to_file( sparse_df, output_directory_path, prefix_string ):
    
    """
    
    Funtion to store a sparse df into a .sparse file 
    
    """
    sparse_file =  os.path.join(output_directory_path, prefix_string + ".sparse")
    sparse_df.to_csv(sparse_file, index = False, header = True, sep = '\t') 
       
       
# ----------------------------------------------------------------------------

def default_smooth(genes_in_network_but_not_in_LFC_list, graph_edges_df, expressed_genes, expression_df):
    
    '''
      set default lfc value for a gene not present in the LFC file but present in the bionetwork

      the default value is calculated from the gene's network neighbours as a (weighted) mean
     
    '''
        
    start_time = time.time()
    
    # Rename columns of the LFC dataframe
    expression_df.columns = ['gene', 'enrichment'] 

    # Initialize dictionary to contain the default LFC value for the 
    # genes not in LFC but present in network
    missing_gene_expression_dict = {}
    
    number_of_genes_with_no_neighbours_in_lfc_gene_list = 0
    
    print('# of genes not in LFC but present in network = ', len(genes_in_network_but_not_in_LFC_list))
    
    test_count = 0
    
    # For each missing gene, search for its neighbours in order to calculate a weighted average expression that
    # is asigned to the missing gene
    
    
    for gene in genes_in_network_but_not_in_LFC_list:

        test_count += 1
        
        # Get subnetwork where the missing gene is the edge start gene
        tmp1 = graph_edges_df[graph_edges_df[0] == gene].iloc[:, [1,2]]  # Neighbour and edge weight
       
        # Get subnetwork where the missing gene is the edge end gene
        tmp2 = graph_edges_df[graph_edges_df[1] == gene].iloc[:, [0,2]] # Neighbour and edge weight
        
        # Concatenate the neighbours and edge weights of the missing gene into single dataframe 
        neighbours = pd.concat( [tmp1, tmp2.rename(columns = { 0 : 1})], ignore_index=True )
        neighbours.columns = ['gene', 'weight']
        neighbours = neighbours.drop_duplicates(subset = 'gene')


        # Find the neighbour genes in the original LFC expressed gene list
        # Note that not all the neighbours will necessarily be in the expressed gene list
        # so len(neighbours) != len(indices)
        indices_of_neighbours_in_expressed_gene_list = matched_index_finder_dict(expressed_genes, neighbours['gene'].tolist())

        # We need to deal with the case where no neighbour is present in the input LFC file
        if len(indices_of_neighbours_in_expressed_gene_list) == 0:
            # None of the neighbouring genes are present in the LFC file 
            number_of_genes_with_no_neighbours_in_lfc_gene_list += 1
            continue
        # Calculate averaged and weighted expression from neighbours:
        else:
        
            # Names of the neighbours that are present in the lfc expressed gene list
            neighbour_genes_in_lfc_gene_list = expressed_genes[indices_of_neighbours_in_expressed_gene_list].tolist()

            # Get the edge weights of the neighbour genes found in the lfc expression genes:
            neighbours_subset = neighbours[neighbours['gene'].isin(neighbour_genes_in_lfc_gene_list)]

            # Ge the expression of the neighbours found in the LFC file:
            expression_subset = expression_df[expression_df['gene'].isin(neighbour_genes_in_lfc_gene_list)]

            # Some genes in the original LFC file have the same HGNC symbol
            # For now, we remove them, but this is a patch, need to find a true solution
            if (len(neighbours_subset) != len(expression_subset)):

               print('DIFFERENT LENGTHS')
               from collections import Counter

               genes_that_are_repeated_in_lfc_nomenclature = [item for item, count in Counter(expression_subset['gene']).items() if count > 1]

               for repeated_gene in genes_that_are_repeated_in_lfc_nomenclature:
                   neighbour_genes_in_lfc_gene_list.remove(repeated_gene)

            
            expression_for_neighbours_in_LFC = expression_df[expression_df['gene'].isin(neighbour_genes_in_lfc_gene_list)]['enrichment']
            neighbour_subset = neighbours[neighbours['gene'].isin(neighbour_genes_in_lfc_gene_list)]
            weight_for_neighbours_in_LFC =  neighbour_subset['weight']

            val = np.average(expression_for_neighbours_in_LFC, weights = (weight_for_neighbours_in_LFC))
        
            # Store in dictionary
            missing_gene_expression_dict[ gene ] = val 
            
                
    missing_gene_expression_df = pd.DataFrame(list(missing_gene_expression_dict.items()), 
                                              columns = ['gene', 'enrichment'])
    missing_gene_expression_df.to_csv( 'missing_gene_expression.csv', sep = '\t' )
    print("--- %s seconds ---" % (time.time() - start_time)) 
    
    return missing_gene_expression_df

# ----------------------------------------------------------------------------    
def sparse_expression( lfc_data_files_paths, 
                      path_to_input_folder, 
                      output_prefix, 
                      graph_genes, 
                      graph_edges_df, 
                      analysis_version  ):
    
    """
    
    Funtion to convert the LFC expression file into a sparse file
    
    """    
       
    output_folder_path_array = []
    
    for lfc_file_path in lfc_data_files_paths:
        
        path_to_output_folder = create_output_directory(lfc_file_path, path_to_input_folder, output_prefix)
        output_folder_path_array.append( path_to_output_folder)       
        
        lfc_df, lfc_genes_list = process_lfc_files( lfc_file_path )   
              
        # Find the gene names that are both in the LFC list and in the network
        # m <- expr[,1] %in% genes # match against bionetwork genes
        indices_of_lfc_genes_present_in_graph = matched_index_finder_dict(graph_genes, lfc_genes_list)
        #matched_graph_gene_names = graph_genes[indices_of_lfc_genes_present_in_graph]
        names_of_lfc_genes_present_in_graph = graph_genes[indices_of_lfc_genes_present_in_graph]

        number_of_lfc_genes_not_in_graph = len(lfc_genes_list) - len(names_of_lfc_genes_present_in_graph)

        print('\t {0} genes in {1} file could not be mapped to the bionetwork'.format(number_of_lfc_genes_not_in_graph, os.path.basename(lfc_file_path)))       

        # Create a data frame containing the names of the LFC genes present in the network
        # and the expression level of the corresponding gene
        lfc_genes_present_in_graph_df = lfc_df[lfc_df['gene'].isin(names_of_lfc_genes_present_in_graph)]
              
        # Filter the wrongly named genes
        bad_pattern = re.compile("\\..+\\.")
        
        #indices_containing_bad_orf_pattern = [i for i in range(len(lfc_genes_list)) if bad_pattern.match(lfc_genes_list[i])]
        indices_containing_bad_orf_pattern = [i for i in range(len(names_of_lfc_genes_present_in_graph)) if bad_pattern.search(names_of_lfc_genes_present_in_graph[i])]
        

        # I think that what we need is to
        for i in indices_containing_bad_orf_pattern:
            wrong_string = lfc_genes_list[i]
            right_string = wrong_string[0 : i + (len(wrong_string) - 2)]
            lfc_genes_list[i] = right_string
            #substr(orfs, 1, a + attributes(a)$match.length - 2)[badnames]
        
        
        # Case where the graph contains weighted edges:
        if ( len(graph_edges_df.columns) == 3 and analysis_version == 2 ):
            
            # Extract gene names that are present in the network but missing
            # in the LFC expression file:
            genes_in_network_but_not_in_lfc = missing_gene_finder_dict(names_of_lfc_genes_present_in_graph, graph_genes)
            
            
            if len(genes_in_network_but_not_in_lfc) > 0:
                
                print('\n'
                      '\t           ****** WARNING ****** \n'
                      '\n {0} bionetwork gene IDs not found in {1} data file \n'
                      '\n Calculating default expression values based on neighbour expression'
                      ' for missing genes'.format(len(genes_in_network_but_not_in_lfc), lfc_file_path))
               
                d = default_smooth( genes_in_network_but_not_in_lfc, 
                                   graph_edges_df, 
                                   lfc_genes_list, 
                                   lfc_df)
                
        
                lfc_df = lfc_df.append(d, ignore_index = True)
        
        # Add the genes for which the default expression was added to the
        # expressed_genes array:
                lfc_genes_list = lfc_df['gene'].to_list()
                
        
        expressed_genes_list = lfc_genes_list
        expression_df = lfc_df

        # Find the gene names that are both in the LFC list and in the network data
        # m <- expr[,1] %in% genes # match against bionetwork genes
        gene_index_in_graph_genes = matched_index_finder_dict(graph_genes, expressed_genes_list)
        
        # Extract the graph gene names corresponding to the indices
        # expr<-expr[m,]
        matched_graph_gene_names = graph_genes[gene_index_in_graph_genes]
        
         # Create a data frame containing the names of the LFC genes present in the network
        # and the expression level of the corresponding gene
        sparse_expression_df = expression_df[expression_df['gene'].isin(matched_graph_gene_names)]
        
        # Substitute the expression_df gene name with the indices of the graph_gene data
        sparse_expression_df = sparse_expression_df.assign( gene = gene_index_in_graph_genes)
        
        
        # Save to output directory
        sparse_df_to_file( sparse_expression_df , path_to_output_folder, 'expression' )
        
        
        
    return output_folder_path_array


    
        
# ----------------------------------------------------------------------------    
    
#def sparse_graph( output_folder_path_array, graph_genes, graph_edges_startGenes, graph_edges_endGenes, graph_edges_weights ):
def sparse_graph( output_folder_path_array, 
                     graph_genes,
                     graph_edges_df, 
                     analysis_version):
      
    """
    
    Funtion to convert the network file into a sparse file
    
    """    
    
    # Get edge start and end gene columns
    graph_edges_startGenes = graph_edges_df[0].values
    graph_edges_endGenes = graph_edges_df[1].values
    
    # Get the indices of the genes in the network as they appear in the geneNames file
    start_gene_matches_idx = matched_index_finder_dict(graph_genes, graph_edges_startGenes)
    end_gene_matches_idx = matched_index_finder_dict(graph_genes, graph_edges_endGenes)
    
    # Index dictionary for unweighted network:
    idx_dict = {'start_gene_idx'   : start_gene_matches_idx, 
                'end_gene_idx'     : end_gene_matches_idx}
    
    sparse_graph_df = pd.DataFrame(idx_dict)
    
    # Case where the graph contains weighted edges:
    if ( len(graph_edges_df.columns) >= 3 and analysis_version == 2 ):
        
        # Get edge scores from third column:
        graph_edge_weights = graph_edges_df[2].values
        
        weighted_graph_dict = { 'start_gene_name' : graph_edges_startGenes,
                                'end_gene_name' : graph_edges_endGenes, 
                                'edge_weights' : graph_edge_weights}
        
        weighted_graph_df = pd.DataFrame(weighted_graph_dict)  
                
        #Get the names of the grap genes that matched the geneNames file
        start_gene_matches = [ graph_genes[m] for m in start_gene_matches_idx ]   
        end_gene_matches = [ graph_genes[m] for m in end_gene_matches_idx ]
        
        # Create a dictionary and the corresponding data frame
        idx_dict['start_gene_name']  = start_gene_matches 
        idx_dict['end_gene_name']    = end_gene_matches
        idx_dict['gene_name_tuples'] = list(zip(start_gene_matches, end_gene_matches))
        
        idx_df = pd.DataFrame(idx_dict)
        
        # Merge the matched index dataframe with the graph dataframe that contains
        # edge weights so that we obtain the matched-index-filtered weights
        sparse_graph_df = pd.merge(idx_df, weighted_graph_df, 
                             on = ['start_gene_name','end_gene_name'], 
                             how = 'left')
        
        sparse_graph_df.drop(['start_gene_name', 'end_gene_name', 'gene_name_tuples'], 
                       inplace = True, axis = 1)
        
         
                
    for path in output_folder_path_array: 
    # Save to output directory
        sparse_df_to_file( sparse_graph_df , path, 'graph' )

                  


# ----------------------------------------------------------------------------    

def convert_to_sparse( path_to_bionetwork, 
                      path_to_input_folder,  
                      output_prefix, 
                      analysis_version, 
                      edge_score_threshold):
    
    '''
    
    Function to write sparse versions of the gene expression levels and the network 
    
    Input: 
        
        => path_to_input_folder -- Path to directory containing LFC dataset to analyse
                                LFC file must be a 2-column .txt file of the 
                                format:   gene name | enrichment 
                                
        => bionetwork -- Name of the bionetwork to be used in the analysis. 
                         This should correspond to the directory where the 
                         two following files are located:
        
                             => graph_genes -- Single column file containing the 
                                               names of the genes present in the network
                                               [ NETWORK GENE NAME FILE ]
                        
                            => graph_edges -- 2-column file containing the names 
                                              of the genes involved in each network edge.
                                              Each row represents an edge ant the 
                                              start and end genes involved
                                             [ NETWORK EDGE GENES NAMES FILE ]
                        
        => output_prefix - Prefix of the directory were the sparse files will be stored
        
        
    Output:
        
        => expression.sparse -- 2-column file containing the index of the LFC 
                                expressed gene as found in the network gene name 
                                file and its corresponding expression level
        
        => graph.sparse -- 2-column file containing the indices of the network 
                           edge genes as
                           found in the network gene name file
    
    '''
    
    # Get processed genes & edges from network
    #graph_genes, graph_edges_startGenes, graph_edges_endGenes, graph_edge_weights = process_network_files( path_to_bionetwork )
    graph_genes, graph_edges_df =   process_network_files( path_to_bionetwork, edge_score_threshold )      
    # Find the files containing the LFC data (should be stored in .txt format)
    # Note: glob.glob() returns an array of pattern matches. Need to iterate!
    lfc_data_files_paths = glob.glob( path_to_input_folder + '*.txt')
    
    # Run write_sparses.py for each of the input LFC.txt files found within the
    # input folder and store the output folder name containing the sparse files
    
    t1 = time.time()
    
    output_folder_path_array = sparse_expression( lfc_data_files_paths, 
                                                path_to_input_folder, 
                                                output_prefix, 
                                                graph_genes,
                                                graph_edges_df, 
                                                analysis_version) 

    t2 = time.time()
    
    print(  
       '\t===============================================\n'
       '\n' 
       '\t Sparse expression created after {:4.2f} seconds \n'.format( t2 - t1 ) )
    
    del t1, t2
    
    # Create sparse graph and, since it is independent of the LFC expression
    # file, store the same sparse graph in all sparse expression output directories:       
    t1 = time.time() 
    sparse_graph(output_folder_path_array, 
                 graph_genes, 
                 graph_edges_df,
                 analysis_version)
    
    t2 = time.time()
    
    
    print(  
       '\t===============================================\n'
       '\n' 
       '\t Sparse graph created after {:4.2f} seconds \n'
       '\n' 
       '\t===============================================\n'.format( t2 - t1 ) )
      
    return output_folder_path_array
    
    
    
    
    
