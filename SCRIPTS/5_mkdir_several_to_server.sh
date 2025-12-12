#!/bin/bash

#/  ------------------------------------------------------------------------------
#/
#/             This script is part of the Rando-SDS v1.0 suite
#/
#/  ------------------------------------------------------------------------------
#
#/ Usage: 5_mkdir_several_to_server.sh 
#/
#/ 
#/ OPTIONS
#/   -h, --help
#/       Script to create folders in server to store LFC files to be copied in future steps
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
#/ | Copyright (C) 2025 Maialen Arrieta-Lobo, Lucile Mégret and Christian Neri  | 
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
		
}

create_folder_for_original_files(){

	# Create a directory where original LFC data file will be copied to
	original_dir_to_create="$path_to_rando_sds/randoSDS_original_experiment"
	echo $original_dir_to_create
	sshpass -p $PASSWORD_SSH ssh $SERVER_ADRESS "mkdir -p $original_dir_to_create"	
	

}

create_folder_for_bs_files(){

	# Iterate over BS replicates
	for ((i=1; i<=$N; i++));do
		
		# Create a directory for each BS replicate and copy it to the server
		dir_to_create="$path_to_rando_sds/randoSDS_BS"$i"/"
		echo $dir_to_create
		sshpass -p $PASSWORD_SSH ssh $SERVER_ADRESS "mkdir -p $dir_to_create"
	done
}


main() {

	declare_variables
	echo "==> Creating folders in  remote server "
	create_folder_for_original_files
	create_folder_for_bs_files
    
}

main "${@}"


