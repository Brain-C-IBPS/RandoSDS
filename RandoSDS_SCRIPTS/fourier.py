#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""

             This module is part of the Rando-SDS v1.0 suite

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
|                                                                            |                                                                          |
============================================================================== 

"""

# ----------------------------------------------------------------------------
#
#                           Import modules
#
# ----------------------------------------------------------------------------

import os
import numpy as np
import scipy
import scipy.sparse
import scipy.linalg
import datetime
import sparse
import pandas as pd

# ----------------------------------------------------------------------------
#
#                           Function definition
#
# ----------------------------------------------------------------------------

 
def read_txt_file( path_to_file_directory, file_name ):
    
    path_to_txt_file = os.path.join( path_to_file_directory, file_name)
    txt_file = open( path_to_txt_file , 'r')
    txt_file = np.loadtxt(txt_file)
     
    return txt_file 

# -----------------------------------------------------------------------------


def diagonalize_graph( output_folder_path ):
                
    '''
    Function that diagonalises the graph within the input data directory
    '''
    
    graph_edge_matrix, n_edges = sparse.get_sparse_graph_edges( output_folder_path  )
    
    #Adjacency and Laplacian matrices
    edges = graph_edge_matrix.iloc[:, 0:2].astype(int)    
    
    # Get the max index of the genes present in the matrix
    n_genes = np.amax( edges.values + 1 )
    
    k = np.empty( len(edges) )
    k.fill(1)
    
    n_edge_matrix_cols = len(graph_edge_matrix.columns)
    
    # Create sparse matrix in Coordinate format.    
    if n_edge_matrix_cols == 2: # Non-weighted network graph
        
        print('\n'
              '\t *************************************************\n'
              '\t         The provided bionetwork is a \n'
              '\t non-weighted network graph (2 column bionetwork) \n'
              '\t ************************************************ \n')
        A = scipy.sparse.coo_matrix( (k, (edges.iloc[:, 0] ,  edges.iloc[:, 1] ) ), 
                                    shape = (n_genes , n_genes )).toarray()
    
    elif n_edge_matrix_cols == 3: # Weighted network graph
        
        print('\n'
              '\t *************************************************\n'
              '\t         The provided bionetwork is a \n'
              '\t weighted network graph (3 column bionetwork) \n'
              '\t ************************************************ \n')
        weights = graph_edge_matrix.iloc[:, 2].astype(float)
        
        #Adjacency and Laplacian matrices
        A = scipy.sparse.coo_matrix( (weights, (edges.iloc[:, 0], edges.iloc[:, 1] ) ), 
                                     shape = (n_genes  , n_genes  )).toarray()
        

    # convert from triangular to interaction matrix and removes parasite same-node interactions if any
    A = A + np.transpose(A) - 2 * np.diag(np.diag(A))
    
    L = np.diag(np.sum(A, axis = 0)) - A # creates Laplacian, weight of node (sum of edge weights)
    (D, eigenvector_matrix) = scipy.linalg.eigh(L)
    
    # eigenvector_matrix is an orthogonal base fot fourier transform, meaning V.t(V).e=e
    D = np.diag(D)
  
    return (eigenvector_matrix, n_genes) # (eigenvector_matrix, number_of_genes)
        
# -----------------------------------------------------------------------------

def smoothe_expression( output_folder_path, eigenvector_matrix, n_genes, gene_extraction_method ):
            
    '''
    Function that smoothes expression
    
    Input:  - Absolute path to expression data
            - Eigenvector matrix
            - Number of genes from diagonalization
            - Number of smoothing iterations
        
    Output: Tuple containing smoothed expression and means
    
    
    '''
    #Load expression
    expression_df = pd.read_csv( os.path.join( output_folder_path, 'expression.sparse'), 
                                    sep = '\t', header = 0)
    
    expression_matrix = expression_df.values
    
    # Check this! Remember that not all the genes that are present in the graph are also present in 
    # the smoothed expression
    #n_genes = len(expression_matrix)
    #print(n_genes)
        
    e = scipy.sparse.csc_matrix( ( expression_matrix[:, 1], ( expression_matrix[:, 0]  , np.zeros(len(expression_matrix) ) ) ), 
                             shape = (n_genes  , 1)).toarray()
   
    
    e[e == 0] = np.mean(e[e != 0])
    
    ft = np.dot(np.transpose(eigenvector_matrix), e) # fourier transform
    
    #Extraction of 100 smoothed expressions
    n_iterations = 40
    
    smoothed_expression = np.zeros((n_genes , n_iterations))
    
    for i in range(n_iterations - 1, -1, -1):
#        ft[round(ngenes*i/N):] = 0
        ft[int(round( (n_genes + 1) * (i + 1) / n_iterations)) : ] = 0
        #ft[int(round( ngenes       * (N-i)   /  N))           : ] = 0
            
        smoothed_expression[:, (n_iterations - 1) - i] = np.dot(eigenvector_matrix, ft)[:,0]
        #smoothexpression[:,i] = sp.dot(P,ft)[:,0]
            
    if ( gene_extraction_method == 1 ) or ( gene_extraction_method == 3):
        means = np.zeros((n_iterations + 1, 1))    
        means[0, 0] = np.mean(smoothed_expression[:, 0])
    
        for i in range(n_iterations):
            means[i + 1, 0] = np.std(smoothed_expression[:,i])
            # MAI MAI means[ i + 1 , 0] = sum((smoothed_expression[:,i] - smoothed_expression[:,i]) / len(smoothed_expression))
    
    if gene_extraction_method == 2:   
        
        means = np.zeros((n_iterations + 1, 2))
        means[0, 0] = np.mean(smoothed_expression[:,0])
        means[0, 1] = np.std(smoothed_expression[:,0])
        
        for i in range(n_iterations):
            means[i + 1, 0] = np.mean(smoothed_expression[:,i])
            means[i + 1, 1] = np.std(smoothed_expression[:,i])
            
        
        
    return (smoothed_expression, means)
    
# -----------------------------------------------------------------------------    
        
def run_smoothing( input_folder_path, output_folder_path_array, same_graph, 
                  gene_extraction_method ):
            
    '''
    Function to diagonalise and smoothe the data
    
    '''
    time_start = datetime.datetime.now()
    
    counter = 0   
    
    for output_folder_path in output_folder_path_array:
        
        # Call diagonalisation:
        if same_graph == 0:
            (V, ngenes) = diagonalize_graph( output_folder_path )
        elif counter == 0:
            (V, ngenes) = diagonalize_graph( output_folder_path )
         
        # Perform smoothing    
        (smoothed_expression, means) = smoothe_expression(output_folder_path, V, ngenes, gene_extraction_method)
        np.savetxt(output_folder_path + '/smoothexpression.txt', smoothed_expression)
        np.savetxt(output_folder_path + '/means.txt', means)        
        counter = counter + 1
    
    time_end = datetime.datetime.now()
    print('\t Diagonalisation and smoothing duration = ' + str(time_end - time_start) )
    

