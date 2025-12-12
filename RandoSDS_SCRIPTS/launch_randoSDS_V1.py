#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""             

        This script is the main script of the Rando-SDS v1.0 suite

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

from datetime import datetime
import parameters
import sds

            
# ----------------------------------------------------------------------------
#
#                           Function definition
#
# ----------------------------------------------------------------------------

def print_info():
    
    '''    
    Function that prints the general information about the code, credits and 
     of the authors, as well as copyright specifications.

    '''

    print('''
    
    
======                                           ======== ========     ========                                   
=======                                          ======== ========     ========                                    
==    ==                         ==              ==       ==     ==    ==
==     ==                        ==              ==       ==      ==   ==                                                         
==     ==                        ==              ==       ==       ==  ==  
==    ==                         ==              ==       ==        == ==   
========   =====     =====    =====  =====       ======== ==        == ========                  
========   =====     =====   ======  =====       ======== ==        == ========                
==   ==   ==   ==   ==   == ==   == ==   ==            == ==       ==        ==
==    ==  ==   ==   ==   == ==   == ==   == ====       == ==       ==        ==       
==     == ==   ==   ==   == ==   == ==   == ====       == ==      ==         ==           
==     == ==   ==   ==   == ==   == ==   ==            == ==     ==          ==                
==     ==  ====  == ==   == =======  =====       ======== ========      ======           
==     ==   ===  == ==   ==  ======  =====       ======== ========      ======    
            
            
    ==========================================================================
    | Source code for performing spectral decomposition of the signal (SDS)  | 
    | against pre-existing biological networks for analysis of bootstrapped  | 
    | genome-wide data sets                                                  |
    ==========================================================================
    
    Version: Rando-SDS v1.0
    Date: 11 March 2023
    
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
    | Copyright (C) 2023 Maialen Arrieta-Lobo and Christian Neri                | 
    ==============================================================================
    |                                                                            |
    |  This work is licensed under the Creative Commons                          | 
    |  Attribution-NonCommercial-NoDerivatives 4.0                               |
    |  International License. To view a copy of this license, visit              |
    |  http://creativecommons.org/licenses/by-nc-nd/4.0/ or send a letter to     |
    |  Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.              |
    |                                                                            |
    ==============================================================================      
       

                        **** DISCLAIMER ****
                   
        The Rando-SDS v1.0 workflow comes with ABSOLUTELY NO WARRANTY
        See the GNU General Public License version 2 for more details
        
        
                        ********************      
   
        
    ''')
    
    # ---------------------------------------------------------------------
    
def print_help():
    
    '''    
    Function that prints the detailed information about the code, its usage options
    and additional specifications about structure and calls. 
    
    '''

    print('''    
          
     ==========================================================================
     |                                                                        |
     |                        Rando-SDS HELP                                  |
     |                                                                        |
     ==========================================================================
     
     Rando-SDS generates Fourier .sif files from the logfoldchange output from 
     any differential analysis that has been bootstrapped for robustness. 
     See Usage section for launching options. 
     
   
     ------------------------   
     Usage
     ------------------------
                         
     python launch_randoSDS_V1.py -b <BIONETWORK> 
                                  -d <deregulation(ALL, ABS, UP, DOWN)> 
                                  -s <species(Mouse/Human)> 
                                  -i <InputFolderName> 
                                  -o <OutputPrefix>
                                  -m <gene_extraction_method> 
                                  -sp <sparsePreprocessing(True/False)>
    
     ------------------------ 
     Additional options
     ------------------------ 
     
       1) It is possible to run the matrix diagonalization only once by adjusting 
          the 'sameGraph' parameter to 1.
       2) The matrix 'V' giving the eigenvectors is not saved (because it is not 
          used further in the workflow)
    
    ------------------------
    Structure
    ------------------------
    
    The parent folder is the folder in which the launch_randoSDS_V1.py code is located. 
    Modules necessary for running the script are also in the parent folder
    The code retrieves the path to the script file to run the analysis.
    The following parent folder structure is mandatory:

    Parent Folder:
     - launch_randoSDS_V1.py
     - fourier.py
     - genextract.py
     - sds.py
     - parameters.py
     - sparse.py   
     - bionetwork(folder)
     - dataExample(folder)
    
    ------------------------
    Storing the bionetwork
    ------------------------
    
    Download the bionetwork and save the network to a folder (e.g. bionetwork). 
    Create 2 files inside the folder as follows
      - A file with the unique set of gene list from the bionetwork and name 
        as 'geneNames.txt'
      - A file with the edge list of the bionetwork and named as 'networkNames.txt'
      
    ------------------------
    INPUT files
    ------------------------
    * The input files are the LFC files generated after performing a differential 
    expression analysis.
    * These LFC files have to be stored inside a folder (e.g. DataExample) 
    within the parent folder.
    * The DataExample folder should only contain LFC files in the format (Gene <tab> LFC) 
    

    ''')     

# ----------------------------------------------------------------------------- 
            
def sparse_condition_checker( processed_params ):
    
    if ( processed_params[5] == False):        
         return False
    
    else :
        return True

# ----------------------------------------------------------------------------- 
            
def fourier_condition_checker( processed_params ):
    
    if ( processed_params[6] == False):        
         return False
    
    else :
        return True

# -----------------------------------------------------------------------------       
    
def main( same_graph ):
    
    print_info()
    print_help()
    
    time_start = datetime.now()
       
    # ---------------------------
    # Processed arguments 
    # ---------------------------
    processed_params = parameters.process_parameters()  
     
    # ---------------------------
    # Instantiate SDS Object 
    # ---------------------------
    sds_obj = sds.spectralDecomposition(processed_params, same_graph)
         
    # ---------------------------------
    # Check if preFourier and smoothing 
    # ---------------------------------
    if sparse_condition_checker( processed_params ):
        sds_obj.convert_input_to_sparse()

        
    if fourier_condition_checker( processed_params ):
        sds_obj.fourier_smoothe_input()
        
    
    # ---------------------------------
    # Perform gene extraction
    # ---------------------------------
    sds_obj.extract_genes()
        
    time_finish = datetime.now()
    
    print('''
    \t              Script ended after {0} 
          
    \t====================================================================         
    
    '''.format( time_finish - time_start))


#*****************************


if __name__ == '__main__':
    
    #sameGraph = 0 if matrix diagonalization is to be run for each folder in the list
    #sameGraph = 1 if matrix diagonalization is to be run only once (i.e. graph.sparse is the same for each folder in the list)
   
    same_graph = 1
    main(same_graph)


#*****************************
