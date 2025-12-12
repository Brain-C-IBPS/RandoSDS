#!/bin/bash
#
#/  ------------------------------------------------------------------------------
#/
#/             This script is part of the Rando-SDS v1.0 suite
#/
#/  ------------------------------------------------------------------------------

#/ Usage: 8_copy_several_from_server.sh 
#/
#/ 
#/ OPTIONS
#/   -h, --help
#/        Script to retrieve post Rando-SDS files, i.e. Rando-SDS results to local computer
#/        Modify server access password and file paths in order to match computer and server setup

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

	# Set server access password		
	PASSWORD_SSH="yourpassword"
	SERVER_ADRESS="user@XXX.XXX.XXX.XXX"

	# Declare path where to Rando-SDS is
	path_to_rando_sds="/path/to/RandoSDS"
	
	# Declare path to server folder storing Rando-SDS result
	path_to_server_folder_sds_results="$SERVER_ADRESS:$path_to_rando_sds"
	
	# Declare paths to folders containing LFC files to be copied
	parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )
	grandparent_path=$( cd "$(dirname "$parent_path")" ; pwd -P )
	path_to_folder_to_store_sds_results="$grandparent_path/ANALYSIS_RESULTS/SDS_NETWORKS"

	if [ ! -d $path_to_folder_to_store_sds_results ]; then
		  mkdir -p $path_to_folder_to_store_sds_results;
	fi

	# Declare disease condition array
	disease_condition_array=("Q50" "Q111" "Q170" "Q175")
	
	# Declare number of BS replicates calculated on step 1 of the workflow
	N=3
		
		
}


copy_original_experiment_results_from_server(){


	for j in "${disease_condition_array[@]}"; do
		
		file="$path_to_server_folder_sds_results/randoSDS_original_experiment/randoSDS/"$j"_LFC_geneSymbols/"
		echo $file
		
		destination="$path_to_folder_to_store_sds_results/randoSDS_original_experiment/"
		
		echo $destination
		if [ ! -d $destination ]; then
		  mkdir -p $destination;
		fi

		# Copy file to destination   
		sshpass -p $PASSWORD_SSH scp -r $file $destination 
		
	done
}

copy_bs_replicates_results_from_server(){

	for ((i=1; i<=$N; i++)); do

		for j in "${disease_condition_array[@]}"; do
		
			file="$path_to_server_folder_sds_results/randoSDS_BS"$i"/randoSDS/"$j"_LFC_geneSymbols/"
			echo $file

			destination="$path_to_folder_to_store_sds_results/randoSDS_BS"$i"/"
			echo $destination
			if [ ! -d $destination ]; then
			  mkdir -p $destination;
			fi

			# Copy file to destination   
			sshpass -p $PASSWORD_SSH scp -r $file $destination
	       done 

	done
}


main() {

	declare_variables
	echo "==> Copying from server"
	copy_original_experiment_results_from_server
	copy_bs_replicates_results_from_server
    
}



main "${@}"



