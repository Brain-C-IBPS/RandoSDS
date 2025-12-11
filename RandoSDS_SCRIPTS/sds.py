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
#                           Import Rando-SDS modules
#
# ----------------------------------------------------------------------------

import sparse
import fourier
import genextract
import os

# ----------------------------------------------------------------------------
#
#                           Create Rando-SDS class
#
# ----------------------------------------------------------------------------

class spectralDecomposition( object ):
    
    def __init__(self, params, same_graph):
        
        self.params = params
        self.bionetwork = params[0] 
        self.species = params[1]  
        self.input_folder_name = params[2]  
        self.input_folder_path = params[3]  
        self.output_prefix = params[4] 
        self.preSparse = params[5]
        self.preFourier = params[6] 
        self.n_fourier_iteration = params[7] 
        self.deregulation_options = params[8] 
        self.analysis_version = params[9]
        self.pval = params[10]
        self.edge_score_thresh = params[11]
        
        self.path_to_main = os.path.dirname(os.path.realpath(__file__))
        self.path_to_bionetwork_folder =  self.path_to_main +'/' + self.bionetwork + '/'
        self.same_graph = same_graph
        
        
        # ---------------------------------------------------------------------

    def convert_input_to_sparse( self ): 
    
        '''
        Function to convert input LFC and network files to sparse format
        
        '''          
        self.output_folder_path_array = sparse.convert_to_sparse(self.path_to_bionetwork_folder, 
                                                                     self.input_folder_path, 
                                                                     self.output_prefix, 
                                                                     self.analysis_version,
                                                                     self.edge_score_thresh)
                                                                           
              
         # -----------------------------------------------------------------------------       
        
    def fourier_smoothe_input( self ): 
            
        '''
        Function to run Fourier smooting on sparse expression data
        
        If the -sp and/or -p preprocessing options are set to False, the code
        assumes that the sparse, smoothexpression and means files have all
        been already generated, and it tries to fetch them from the folders 
        within the input folder path that start with the output_prefix. 
        
        ''' 

        # First check if the output folder array attribute has been assigned
        # If not, perform the assignment
        try:
            
            self.output_folder_array
            
        except AttributeError:
            
            path_to_output_folder =  os.path.join(self.input_folder_path, self.output_prefix)
            
            directories_in_output_folder = [ d for d in os.listdir(path_to_output_folder) 
                            if os.path.isdir(os.path.join(path_to_output_folder, d))]
            
            
            self.output_folder_path_array = [ os.path.join(path_to_output_folder, folder) for folder in directories_in_output_folder]
            
        # Run smoothing
        fourier.run_smoothing( self.input_folder_path, 
                                        self.output_folder_path_array, 
                                        self.same_graph, 
                                        self.analysis_version)
        
        
        # -----------------------------------------------------------------------------
        
    def extract_genes( self ):
    
    
        '''
        Function that extracts the expression of the genes from input data 
                
        If the -sp and/or -p preprocessing options are set to False, the code
        assumes that the sparse, smoothexpression and means files have all
        been already generated, and it tries to fetch them from the folders 
        within the input folder path that start with the output_prefix. 
        
        '''
	# First check if the output folder array attribute has been assigned
        # If not, perform the assignment
        
        try:
            
            self.output_folder_array
            
        except AttributeError:
            
            path_to_output_folder =  os.path.join(self.input_folder_path, self.output_prefix)
            
            directories_in_output_folder = [ d for d in os.listdir(path_to_output_folder) 
                            if os.path.isdir(os.path.join(path_to_output_folder, d))]
            
            
            self.output_folder_path_array = [ os.path.join(path_to_output_folder, folder) for folder in directories_in_output_folder]
            
         

        genextract.extract_genes( self.output_folder_path_array, 
                                          self.bionetwork, 
                                          self.same_graph, 
                                          self.deregulation_options, 
                                          self.n_fourier_iteration, 
                                          self.analysis_version, 
                                          self.pval)	
 
                                                    
         
     
