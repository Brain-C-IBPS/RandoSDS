#!/bin/bash

#/  ------------------------------------------------------------------------------
#/
#/             This script is part of the Rando-SDS v1.0 suite
#/
#/  ------------------------------------------------------------------------------

#/ Usage: 7_randoSDS_run_SDS_multiple_files.sh 
#/
#/ 
#/ OPTIONS
#/   -h, --help
#/       Script to run Rando-SDS
#/       Modify server access password and file paths in order to match computer and server setup

#/ Contributors (alphabetic order): Maialen Arrieta-Lobo (1) 
#/                                  Lucile Mégret (1)
#/                                  Christian Neri (1)
    
#/ Laboratories:   
#/ (1) Institut de Biologie Paris Seine (IBPS), UMR CNRS 8256, 
#/     Team Brain-C (Compensation systems in neurodegenerative diseases and aging)

#/ Affiliations:
#/ (1) Sorbonnes Université
#/ (1) Centre National de la Recherche Scientifique (CNRS)
#/ (1) Institut National de la Santé et de la Recherche Médicale (INSERM)

#/ ==============================================================================
#/ | Copyright (C) 2023 Maialen Arrieta-Lobo, Lucile Mégret and Christian Neri  | 
#/ ==============================================================================
#/ |                                                                            |
#/ |  This work is licensed under the Creative Commons                          | 
#/ |  Attribution-NonCommercial-NoDerivatives 4.0                               |
#/ |  International License. To view a copy of this license, visit              |
#/ |  http://creativecommons.org/licenses/by-nc-nd/4.0/ or send a letter to     |
#/ |  Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.              |
#/ |                                                                            |
#/ ============================================================================== 


#/  ------------------------------------------------------------------------------


declare_variables() {

	# Set common variables for Rando-SDS script
	species='Mouse'
	bionetwork='mousenet_V2_filtered'
	output_prefix='randoSDS'
	gene_extraction_method='1'
	sparse_preprocessing='True'
	prefourier='True'

	# Number of BS replicates	
	N=3

}



run_randoSDS_bs_replicates() {

	for ((i=1; i<=$N; i++)); do

		lfc_folder="randoSDS_BS$i"
		echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
		echo ""
		echo "        Running Rando-SDS for "$lfc_folder" BS replicate"
		echo ""
		echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
		echo ""
		echo ""
		python launch_randoSDS_V1.py -b $bionetwork -d ALL -s $species -i $lfc_folder -o $output_prefix -m $gene_extraction_method -sp $sparse_preprocessing
	done

}

run_randoSDS_original_experiment() {

	lfc_folder="randoSDS_original_experiment/"
	echo "+++++++++++++++++++++++++++++++++++++++++++++++++"
	echo ""
	echo "        Running Rando-SDS for "$lfc_folder""
	echo ""
	echo "+++++++++++++++++++++++++++++++++++++++++++++++++"
	echo ""
	echo ""
	python launch_randoSDS_V1.py -b $bionetwork -d ALL -s $species -i $lfc_folder -o $output_prefix -m $gene_extraction_method -sp $sparse_preprocessing

}

main() {

	declare_variables
	run_randoSDS_bs_replicates
	run_randoSDS_original_experiment
    
}



main "${@}"


