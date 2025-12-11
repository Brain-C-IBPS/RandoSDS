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

import argparse
import os


# ----------------------------------------------------------------------------
#
#                           Function definition
#
# ----------------------------------------------------------------------------

def str2bool(v):
    if isinstance(v, bool):
       return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

# -----------------------------------------------------------------------------

def get_args():
    
    #if isinstance(args, list):
    parser = argparse.ArgumentParser(description = 'Please launch using the following format:\n\
                                                    python3 launch_SDS.py \n\
                                                    -b <bionetworkDirectoryName> \n\
                                                    -d ABS \n\
                                                    -s <Species> \n\
                                                    -i <InputFolder> \n\
                                                    -o <Outputfileprefix>\n\
                                                    -sp <SparsePreProcessing>\n\
                                                    -p <preFourier>\n\
                                                    -n <numberOfFourierIterations>\n\
                                                    -m <geneExtractionMethod>\n\
                                                    -pval <pvalThresholdV2>\n\
                                                    -e <edgeScoreThresholdV2>\n\
                                                    ')
    
    # Read the name of the folder were the Bionetwork is stored 
    parser.add_argument('-b', 
                        type = str, 
                        dest = 'bionetwork',
                        help = ('Specify the bionetwork: Custom bionetworks can be given. '
                                'Please refer to README.txt file before adding a custom network.'
                                'A example bionetwork (mousenetV2) is provided within the zip package'))
    
    parser.add_argument('-d', 
                        type = str, 
                        dest = 'deregulation_type', 
                        default = 'ABS',
                        help = ('Select option for analysis: U = up_regulated genes, '
                               'D = down_regulated, ABS = both [default = ABS]'))
           
    # Read the species corresponding to the analysis results
    parser.add_argument('-s', 
                        type = str, 
                        dest = 'species',
                        help = 'Species to consider for the analysus [Mouse or Human]')
    
    # Read the input folder name to read data from
    parser.add_argument('-i', 
                        type = str, 
                        dest = 'input_folder_name',
                        help = 'Name of input folder containing gene expression data')
    
    # Read the prefix for output Fourier transform files
    parser.add_argument('-o', 
                        type = str, 
                        dest = 'output_prefix',
                        help = ('Prefix string for the output .sif files '))
    
    # Read sparse preprocessing option
    parser.add_argument('-sp', 
                        type = str2bool, 
                        dest = 'preSparse', 
                        default = True,
                        help = 'Run sparse preprocessing [default = True].')
    
    # Read preFourier option
    parser.add_argument('-p', 
                        type = str2bool, 
                        dest = 'preFourier', 
                        default = True,
                        help = 'Run Prefourier [default = True].')
    
    # Read number of Fourier transform iterations
    parser.add_argument('-n', 
                        type = int, 
                        dest = 'n_fourier_iteration', 
                        default = 10,
                        help = ('Number of Fourier iterations for expression '
                                'smoothing [default = 10]'))
    
    # Read the gene extraction method
    parser.add_argument('-m', 
                        type = int, 
                        dest = 'gene_extraction_method',
                        default = 1,
                        help = ('Method to select threshold for deregulated genes: 1, 2  or 3 [default = 1]. '
                                '1: deregulation threshold = mean smoothed expression +/- std'
                                '2: deregulation threshold = ppf of  mean smoothed expressions'
                                '3: deregulation threshold  = mean smoothed expression'
                                'Please refer to README.txt file for more details on each method.'
                                ))
    
    # Read the gene extraction method
    parser.add_argument('-pval', 
                        type = float, 
                        dest = 'pval',
                        default = 0.05,
                        help = ('Threshold pval for gene_extraction_method = 2'
                                '[default = 0.05]'
                                ))
    
    # Read the gene extraction method
    parser.add_argument('-e', 
                        type = float, 
                        dest = 'edge_score_thresh',
                        default = None,
                        help = ('Threshold for weighted network edge score. Only edges'
                                'whose score is above the value will be considered'
                                'for the analysis'
                                ))
    
    args = parser.parse_args()
    
    return parser, args
   

# -----------------------------------------------------------------------------
    
def parameter_processor(parser, args):

    processed_params = []
    
    # ---------------------------
    # Process input arguments 
    # ---------------------------
        
    # Bionetwork:
    if args.bionetwork:
        bionetwork = args.bionetwork
        processed_params.append(bionetwork)
    else: 
        parser.error("PARSER ERROR -- Missing bionetwork argument")
        exit(parser.print_help())
    
    # Species:    
    if args.species:
        species = args.species
        processed_params.append(species)
    else: 
        parser.error("PARSER ERROR -- Missing species argument")
        exit(parser.print_help())
    
    # Input folder:    
    if args.input_folder_name:
        input_folder_name = args.input_folder_name
        processed_params.append(input_folder_name)    
                                       
        input_folder_path = os.path.dirname(os.path.realpath(__file__))  +'/' + input_folder_name + '/'         
        processed_params.append(input_folder_path)    
                                       
        if not os.path.isdir(input_folder_path):
            parser.error('PARSER ERROR -- Input folder does not exist. Please check path.\n')    
            exit(parser.print_help())
    else: 
        parser.error("PARSER ERROR -- Missing input folder name argument")
        exit(parser.print_help())        
    
    
    if args.output_prefix:
        output_prefix = args.output_prefix
        processed_params.append(output_prefix) 
        
        #if os.path.isdir(input_folder_path + '/' + output_prefix):
         #   parser.error(('PARSER ERROR -- Output folder {0} already exists. '
         #                'Please change output prefix.\n').format(output_prefix) )   
        #    exit(parser.print_help())
        
    else: 
        parser.error("PARSER ERROR -- Missing output prefix argument")
        exit(parser.print_help()) 
    
    preSparse = args.preSparse
    processed_params.append(preSparse) 
    
    preFourier = args.preFourier
    processed_params.append(preFourier) 
    
    n_fourier_iteration = args.n_fourier_iteration
    processed_params.append(n_fourier_iteration) 
   
    deregulation = args.deregulation_type   
    if deregulation == 'ABS':
        deregulation_options = ['ABS']
    elif deregulation == 'U':
        deregulation_options = ['U']
    elif deregulation =='D':
        deregulation_options = ['D']
    elif deregulation == 'ALL':
        deregulation_options = ['U','D','ABS']
    elif deregulation == 'UD':
        deregulation_options = ['U','D']
    processed_params.append(deregulation_options) 
    
    gene_extraction_method = args.gene_extraction_method
    processed_params.append(gene_extraction_method) 
    
    pval = args.pval
    processed_params.append(pval)
    
    edge_score_thresh = args.edge_score_thresh
    processed_params.append(edge_score_thresh)
    
    
    return processed_params 

# -----------------------------------------------------------

class confirmator( object ):
    
    def __init__( self, message, error_message ):
        
        self.message = message
        self.error_message = error_message
        self.confirmation_options = ('yes', 'y', 'YES', 'Yes')
        
    def confirm( self ):
        
        #confirmation = input(self.message)
        confirmation = 'yes'

        if confirmation in self.confirmation_options:    
            return True
        else:
            print(self.error_message)
            return False
            
# -----------------------------------------------------------------------------    

def confirm_parameter_choice(processed_params):
    
    
    confirm_parameters = False
    
    print(  
       '\t===============================================\n'
       '\n' 
       '\t Rando-SDS will run with the following parameters: \n' 
       '\n' 
       '\t===============================================\n'
       '\n'
       '\t Bionetwork = {0}     \n' 
       '\t Species = {1}        \n'
       '\t Input Folder = {2}   \n'
       '\t Output Prefix = {3}  \n'     
       '\t Sparse preprocessing = {4}        \n'
       '\t Run preFourier = {5}  \n'
       '\t Number of Fourier iterations for smoothing = {6}  \n'       
       '\t Deregulation type = {7}        \n'     
       '\t Gene extraction method = {8}        \n'  
       '\n'
       '\t================================================\n'    
        .format(processed_params[0], 
                processed_params[1], 
                processed_params[2], 
                processed_params[4], 
                processed_params[5], 
                processed_params[6], 
                processed_params[7],
                processed_params[8],
                processed_params[9]))
    
    message = ('\n'
               '\t Please confirm above options to move forward: \n'
              ' \t [yes / no] or [y / n] '   )

    error_message = '\t Discarding input parameters and aborting!'
    
    parameter_confirmator = confirmator(message, error_message)
    parameter_confirmation = parameter_confirmator.confirm()
             
    if parameter_confirmation:
        confirm_parameters = True
    
    return confirm_parameters
    

# -----------------------------------------------------------------------------
        
def check_directory_occupancy(processed_params):
    
    # Generates the file names in the input directory tree:
    # New in V3: if the input directory contains a folder whose name is the same
    # as the output prefix parsed, abort
    # Otherwise, continue

    input_folder_path = processed_params[3]
    
    if os.listdir(input_folder_path) :
        
        output_prefix = processed_params[4]
        analysis_type = processed_params[7]
        
        for i in analysis_type:        
            
            output_folder_name = output_prefix + i
            
            if output_folder_name in os.listdir(input_folder_path):
                
                print('FATAL ERROR -- Output folder already exists. '
                     'Relaunch analysis for unique output prefix')
                
                return False
                
            else:
                continue
    else: 
        return True
  
        
# -----------------------------------------------------------------------------
  
    
def process_parameters():
    
    
    # ---------------------------
    # Get arguments from command line
    # ---------------------------       
    parser, args = get_args()
    
    # ---------------------------
    # Process arguments 
    # ---------------------------       
    processed_parameters = parameter_processor(parser, args) 
    
    if not confirm_parameter_choice(processed_parameters):
        print('ABORTED -- Parameters not confirmed. '
              'Relaunch analysis for correct set of parameters')
        exit()
    else:
        return processed_parameters
    
        
# ----------------------------------------------------------------------------
#
#                           Main function
#
# ----------------------------------------------------------------------------
    
    
if __name__ == '__process_parameters__':
    
      process_parameters()
    
    
     
