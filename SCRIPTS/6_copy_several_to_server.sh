#!/bin/bash

#/  ------------------------------------------------------------------------------
#/
#/             This script is part of the Rando-SDS v1.0 suite
#/
#/  ------------------------------------------------------------------------------

#/ Usage: 6_copy_several_to_server.sh 
#/
#/ 
#/ OPTIONS
#/   -h, --help
#/       Script to scp original LFC files and BS replicate LFC files from local to server
#/       Modify server access password and file paths in order to match computer and server setup

#/ Contributors (alphabetic order): Maialen Arrieta-Lobo (1) 
#/                                  Christian Neri (1)
    
#/ Laboratories:   
#/ (1) Institut de Biologie Paris Seine (IBPS), UMR CNRS 8256, 
#/     Team Brain-C (Compensation systems in neurodegenerative diseases and aging)

#/ Affiliations:
#/ (1) Sorbonnes Université
#/ (1) Centre National de la Recherche Scientifique (CNRS)
#/ (1) Institut National de la Santé et de la Recherche Médicale (INSERM)

#/ ==============================================================================
#/ |      Copyright (C) 2023 Maialen Arrieta-Lobo and Christian Neri            | 
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

	# Modify your server access password
	PASSWORD_SSH="yourpassword"
	SERVER_ADRESS="user@XXX.XXX.XXX.XXX"

	# Declare path where to Rando-SDS is
	path_to_rando_sds="/path/to/RandoSDS"
	
	# Declare disease condition array
	disease_condition_array=("Q50" "Q111" "Q170" "Q175")
	
	# Declare number of BS replicates calculated on step 1 of the workflow
	N=3	
	
	# Declare paths to folders containing LFC files to be copied
	parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )
	grandparent_path=$( cd "$(dirname "$parent_path")" ; pwd -P )
	path_to_folder_containing_original_data="$grandparent_path/ANALYSIS_RESULTS/PRE_PROCESSING_ORIGINAL_EXPERIMENT"
	path_to_folder_containing_bs_replicates="$grandparent_path/ANALYSIS_RESULTS/PRE_PROCESSING_BS_REPLICATES"
}

copy_original_lfc_files_to_server(){
	echo "==> Copying to server"
	for j in "${disease_condition_array[@]}"; do

		file="$path_to_folder_containing_original_data/DESEq2_RESULTS/lfc_"$j"_geneSymbols.txt"
		echo $file
		
		# The path in $destination needs to be exactly the same as the one defined in script number 5 (i.e. 5_mkdir_several_to_server.sh)
		destination="$SERVER_ADRESS:$path_to_rando_sds/randoSDS_original_experiment/"$j"_LFC_geneSymbols.txt"	
		sshpass -p $PASSWORD_SSH scp -r $file $destination 
	done

}

copy_bs_lfc_files_to_server(){
	echo "==> Copying to server"
	for ((i=1; i<=$N; i++)); do

		for j in "${disease_condition_array[@]}"; do
			file="$path_to_folder_containing_bs_replicates/BSreplicate_"$i"/DESEq2_RESULTS/lfc_"$j"_geneSymbols.txt"
			echo $file
			destination="$SERVER_ADRESS:$path_to_rando_sds/randoSDS_BS"$i"/"$j"_LFC_geneSymbols.txt"
			sshpass -p $PASSWORD_SSH scp -r $file $destination 
		done
		
		
	done

}



main() {

	declare_variables
	copy_original_lfc_files_to_server
	copy_bs_lfc_files_to_server
    
}



main "${@}"


