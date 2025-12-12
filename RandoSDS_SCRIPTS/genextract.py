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

import numpy as np
import os
import itertools
import fourier
import sparse 
from scipy.stats import norm
from textwrap import dedent


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
           
def multidim_intersect( arr1, arr2):
    
    '''
    Calculate intersection of two arrays of different dimensions and return 
    unique values in each of them
    '''
    return  np.array([x for x in set(tuple(x) for x in arr1) & set(tuple(x) for x in arr2)])  

    
# -----------------------------------------------------------------------------
           
def save_edges_subset( output_folder_path,
                      edges_subset, 
                      graph_genes,
                      selected_option, 
                      i ): 
    
    '''
    Save edges subset to .txt and .sif files
    
    '''
    
    save_to_folder = os.path.join(output_folder_path, selected_option)
    
    save_to_file = os.path.join(save_to_folder, 
                                    selected_option +
                                    '_' + str(i))
    
    #np.savetxt(  save_to_file + '.txt' , edges_subset, fmt = '%i')

    f = open( save_to_file + '.sif', 'w')
    f2 = open(save_to_file + '.txt', 'w')
    
    for e in range(len(edges_subset)):
        e0, e1 = int(edges_subset[e, 0]), int(edges_subset[e, 1])
        f.write(str(graph_genes[e0]) + ' IN ' + str(graph_genes[e1]) +'\n')
        f2.write(graph_genes[e0] + '\n')
        f2.write(graph_genes[e1] + '\n')
    f.close()
    f2.close()

#f.write(geneNames_list[graph[e, 0]].split("\t")[0] + ' IN ' + geneNames_list[graph[e,1]].split("\t")[0]  +'\n')

# -----------------------------------------------------------------------------       

def extract_genes_v1( output_folder_path, 
                     n_genes, 
                     graph_genes, 
                     graph_edge_tuples,
                     smoothexpression, 
                     means, 
                     gene_expr_options, 
                     n_fourier_iteration ):
    
    
    '''
    Function that extracts the expression of the genes from input data 
    
    '''
     
    for selected_option in gene_expr_options:
        
        genes_dereg = np.zeros((n_genes + 1, 40)) # This should be n_iterations, the number of Fourier transforms?
        
        threshold_pos = means[0] + means[1:] 
        threshold_pos_arr = np.ones((n_genes + 1, 40)) * threshold_pos
        threshold_neg = means[0] - means[1:]
        threshold_neg_arr = np.ones((n_genes + 1, 40)) * threshold_neg
                    
        if selected_option == 'U':
           genes_dereg[smoothexpression > threshold_pos_arr] = 1
           
        elif selected_option == 'D':
           genes_dereg[smoothexpression < threshold_neg_arr] = 1   
           
        elif selected_option == 'ABS':
           genes_dereg[smoothexpression > threshold_pos_arr] = 1
           genes_dereg[smoothexpression < threshold_neg_arr] = 1        
                
        genes_dereg_index = [list(np.where(genes_dereg[:,i])[0]) for i in range(40)]
                       
        for i in [n_fourier_iteration]:
            
            
            ''' Code to store the 3-column format similar to V2'''
            genes_dereg_tosave = np.array(["NS"] * (n_genes + 1))            
            genes_dereg_tosave[smoothexpression[:, i] > threshold_pos[i]] = "UP"
            genes_dereg_tosave[smoothexpression[:, i] < threshold_neg[i]] = "DOWN"
                
                
            f = open(os.path.join( output_folder_path,'smoothexpression_' + str(i) + '.txt'), 'w')
            
            for j in range(n_genes):
                
                f.write("\t".join( (graph_genes[j], str(smoothexpression[j,i]), genes_dereg_tosave[j]) ) + "\n")
                
            f.close()
            
            l = genes_dereg_index[i]
            
            # Extract dereulated genes:
            if len(l) > 0:
                candidateGenePairs = np.array([np.array([p ,q ]) for (p,q) in itertools.permutations(l, 2)])
                
                
                # R based 1 indexing                
                graph_edge_tuples_index_R = graph_edge_tuples + 1
                edges_subset = multidim_intersect(candidateGenePairs, graph_edge_tuples_index_R)
                
               
                msg = dedent('''\t====================================================================
                             \t       {0} NETWORK                                                                                                    
            Number of deregulated genes = {1}                             
            Number of edges in interaction network (IN) subset = {2}      
                                                                          
        ==================================================================== '''
    .format( selected_option, len(l),  len(edges_subset) ))
                
                print(msg)
                
                # """""""""""""""""""""""""""""""             
                # Save the edges subset to files:
                #               
                save_edges_subset( output_folder_path,
                              edges_subset, 
                              graph_genes,
                              selected_option,
                              i )
                
            else:
                print('No deregulated genes')    
            
           
           
                
# -----------------------------------------------------------------------------       
    
def extract_genes_v2( output_folder_path, 
                                n_genes, 
                                graph_genes, 
                                graph_edge_tuples, 
                                smoothexpression, 
                                means, 
                                gene_expr_options, 
                                n_fourier_iteration,
                                pval = 0.05 ): # default pvalue set to 0.05
    

    '''
    Function that extracts the expression of the genes from input data 
    for a threshold based a normal distribution standard deviation multiplier 
    
    '''
    
    for i in [n_fourier_iteration]:
        
        # deregulation for positive and negative deregulation // normal distribution
        '''
        The method norm.ppf() takes a percentage and returns a standard deviation 
        multiplier for what value that percentage occurs at. 
        
        For instance, 
            
            norm.ppf(0.95, loc=172.7815, scale=4.1532)

        will return a value (that functions as a 'standard-deviation multiplier') 
        marking where 95% of data points would be contained if our data is a 
        normal distribution.

        To get the exact number, we take the norm.ppf() output and multiply it 
        by our standard deviation for the distribution in question.
        
        interval_value = std * ppf
        '''
        
        threshold_pos = norm.ppf(1.0 - pval, means[i, 0], means[i, 1])
        threshold_neg = norm.ppf(pval, means[i, 0], means[i, 1])
        
        
            
        genes_dereg = np.array(["NS"] * (n_genes + 1))
        genes_dereg[smoothexpression[:,i] > threshold_pos] = "UP"
        genes_dereg[smoothexpression[:,i] < threshold_neg] = "DOWN"    
        
        print('UP', len(np.where(smoothexpression[:, i] > threshold_pos)[0]))
        print('DOWN', len(np.where(smoothexpression[:, i] < threshold_neg)[0]))
        print('ABS', len(np.where(np.logical_or((smoothexpression[:,i] > threshold_pos),(smoothexpression[:,i] < threshold_neg)))[0]))

            
        f = open(os.path.join( output_folder_path,'smoothexpression_' + str(i) + '.txt'), 'w')
        
        for j in range(n_genes):
            
            f.write("\t".join( (graph_genes[j], str(smoothexpression[j,i]), genes_dereg[j]) ) + "\n")
            
        f.close()
        
        for selected_option in gene_expr_options:        
            
            print('i = ' + str(i) + ', ' + selected_option)
            
            if selected_option == 'U': 
                
                l = np.where(smoothexpression[:, i] > threshold_pos)[0]
                
            elif selected_option == 'D':
                
                l = np.where(smoothexpression[:, i] < threshold_neg)[0]
                
            elif selected_option == 'ABS':
                
                l = np.where(np.logical_or((smoothexpression[:,i] > threshold_pos),(smoothexpression[:,i] < threshold_neg)))[0]
            
                
            print(len(l))
            
            if len(l) > 0:
                
                    
                # V2:
                # l is python-based indexed so we can directly search for those
                # positions in our graph:
                edges_subset_idx_python = np.where(np.logical_and(np.isin(graph_edge_tuples[:, 0], l), 
                                                       np.isin(graph_edge_tuples[:, 1], l)))[0] # graph edges between two deregulated genes
                
                  
                msg = dedent('''\t====================================================================
                                 \t       {0} NETWORK                                 
            Number of deregulated genes = {1}                             
            Number of edges in interaction network (IN) subset = {2}      
                                                                          
        ==================================================================== '''
    .format( selected_option, len(l),  len(edges_subset_idx_python) ))
                print(msg)
                
                # """""""""""""""""""""""""""""""""""            
                # Save the edges subset to files:
                #    
                f = open(os.path.join( output_folder_path, selected_option, selected_option + '_' + str(i) + '.sif'),'w')
                f2 = open(os.path.join( output_folder_path, selected_option, selected_option + '_' + str(i) + '.txt'),'w')
                
                for e in edges_subset_idx_python:
                    
                     f.write(graph_genes[graph_edge_tuples[e,0]].split("\t")[0] + ' IN ' + graph_genes[graph_edge_tuples[e,1]].split("\t")[0]  +'\n')
                     f2.write(graph_genes[graph_edge_tuples[e,0]] + '\n')
                     f2.write(graph_genes[graph_edge_tuples[e,1]] + '\n')
                    
                f.close()     
                f2.close()
        
# -----------------------------------------------------------------------------       
    
def extract_genes_v3( output_folder_path, 
                     n_genes, 
                     graph_genes, 
                     graph_edge_tuples,
                     smoothexpression, 
                     means, 
                     gene_expr_options, 
                     n_fourier_iteration ):
    
    
    '''
    Function that extracts the expression of the genes from input data 
    for a threshold equivalent to the mean
    
    '''
     
    for selected_option in gene_expr_options:
        
        genes_dereg = np.zeros((n_genes + 1, 40)) # This should be n_iterations, the number of Fourier transforms?
        
        threshold = means[0]
                    
        if selected_option == 'U':
           genes_dereg[smoothexpression > threshold] = 1
           
        elif selected_option == 'D':
           genes_dereg[smoothexpression < threshold] = 1   
           
        elif selected_option == 'ABS':
           genes_dereg[smoothexpression > threshold] = 1
           genes_dereg[smoothexpression < threshold] = 1        
                
        genes_dereg_index = [list(np.where(genes_dereg[:,i])[0]) for i in range(40)]
                 
        for i in [n_fourier_iteration]:
            
            l = genes_dereg_index[i]
            
            
            # Extract deregulated genes:
            if len(l) > 0:
                
                candidateGenePairs = np.array([np.array([p ,q ]) for (p,q) in itertools.permutations(l, 2)])
                   
                # So that the indices correspond to those in V1 (R based 1 indexing)                   
                graph_edge_tuples_index_R = graph_edge_tuples + 1
                edges_subset = multidim_intersect(candidateGenePairs, graph_edge_tuples_index_R)
                    
                   
                msg = dedent('''\t====================================================================
                                                 
                             \t       {0} NETWORK                                 
            Number of deregulated genes = {1}                             
            Number of edges in interaction network (IN) subset = {2}      
                                                                          
        ==================================================================== '''
    .format(selected_option, len(l),  len(edges_subset) ))
                
                print(msg)
                
                # """""""""""""""""""""""""""""""""""            
                # Save the edges subset to files:
                #               
                save_edges_subset( output_folder_path,
                                  edges_subset, 
                                  graph_genes,
                                  selected_option,
                                  i )
                
            else:
                print('No deregulated genes')    
                
                
                
          
# -----------------------------------------------------------------------------       

def extract_genes( output_folder_path_array,  
                          path_to_bionetwork, 
                          same_graph, 
                          gene_expr_options, 
                          n_fourier_iteration, 
                          extraction_method, 
                          *pval):
    

    print('\n \t Starting gene extraction \n')

     # Get the gene names from the bionetwork directory:
    graph_genes, _ = sparse.process_network_files( path_to_bionetwork )  
                
    counter = 0
    #Loads the graph (only once if possible)
    if (( same_graph == 0 )  or ( same_graph == 1 and counter == 0)):               
        sparse_graph_df, _ = sparse.get_sparse_graph_edges( output_folder_path_array[0])
    
    
    sparse_graph_tuples = sparse_graph_df.to_numpy( dtype = int )
    
    
    # Get the max index of the genes present in the matrix
    n_genes = int(np.amax( sparse_graph_tuples ))
    
    for output_folder_path in output_folder_path_array:
                
        #Loads smoothexpression and means files
        smoothexpression = fourier.read_txt_file( output_folder_path, 'smoothexpression.txt' )
        means = fourier.read_txt_file( output_folder_path, 'means.txt')
        
        print('\n \t Running for files in {0}\n'.format(output_folder_path ))

        # """""""""""""""""""""""""""""""""""
        
        # Create folders to store the extracted genes depending on the
        # chosen analysis type:
        for selected_option in gene_expr_options: # empty folder creation only
            if selected_option == 'U':
                sparse.check_directory_exists( os.path.join(output_folder_path, 'U') )                   
                    
            elif selected_option == 'D':
                sparse.check_directory_exists( os.path.join(output_folder_path, 'D') )               
                
            elif selected_option == 'ABS':
                sparse.check_directory_exists( os.path.join(output_folder_path, 'ABS') )      
  
    
        # """""""""""""""""""""""""""""""""""
        # V1 threshold:
            
        if extraction_method == 1 : 
            
            extract_genes_v1( output_folder_path,
                                n_genes, 
                                graph_genes, 
                                sparse_graph_tuples, 
                                smoothexpression, 
                                means, 
                                gene_expr_options, 
                                n_fourier_iteration )
            
        # """""""""""""""""""""""""""""""""""
        # V2 threshold:
            
        elif extraction_method == 2 : 
            
            extract_genes_v2( output_folder_path, 
                                n_genes, 
                                graph_genes, 
                                sparse_graph_tuples, 
                                smoothexpression, 
                                means, 
                                gene_expr_options, 
                                n_fourier_iteration )
            
        # """""""""""""""""""""""""""""""""""
        # V3 threshold:
            
        elif extraction_method == 3 : 
            
            extract_genes_v3( output_folder_path,
                                n_genes, 
                                graph_genes, 
                                sparse_graph_tuples, 
                                smoothexpression, 
                                means, 
                                gene_expr_options, 
                                n_fourier_iteration )

    
    
# ----------------------------------------------------------------------------
#
#                           Main function
#
# ----------------------------------------------------------------------------
    
    
if __name__ == '__extract_genes__':
    
   extract_genes( )
    
        
    
    
    
    
    
    
    
    
    
    
    
