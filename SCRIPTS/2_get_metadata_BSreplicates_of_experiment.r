#-----------------------------------------------------------------------------
#
#            R script of the Rando-SDS V1.0 workflow to
#            preprocess BS replicates of read count data
#
#-----------------------------------------------------------------------------
#
# Contributors (alphabetic order): Maialen Arrieta-Lobo (1) 
#                                  Christian Neri (1)
#   
# Laboratories:   
# (1) Institut de Biologie Paris Seine (IBPS), UMR CNRS 8256, 
#     Team Brain-C (Compensation systems in neurodegenerative diseases and aging)
#
# Affiliations:
# (1) Sorbonnes Université
# (1) Centre National de la Recherche Scientifique (CNRS)
# (1) Institut National de la Santé et de la Recherche Médicale (INSERM)

# ==============================================================================
# |      Copyright (C) 2023 Maialen Arrieta-Lobo and Christian Neri            | 
# ==============================================================================
# |                                                                            |
# |  This work is licensed under the Creative Commons                          | 
# |  Attribution-NonCommercial-NoDerivatives 4.0                               |
# |  International License. To view a copy of this license, visit              |
# |  http://creativecommons.org/licenses/by-nc-nd/4.0/ or send a letter to     |
# |  Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.              |   
# |                                                                            |                                                                         |
# ============================================================================== 

#-----------------------------------------------------------------------------

# Loading R libraries
library("stringr")

# Get script directory
this.dir <- dirname(sys.frame(1)$ofile)
print(this.dir)

# Path to parent directory
parent_dir <- dirname(this.dir)

# Path to BS replicates of experiment from parent directory
bs_replicates_dir <- file.path(paste0(parent_dir, '/ANALYSIS_RESULTS/BOOTSTRAP_REPLICATES_OF_EXPERIMENT'))

# Get CSV files of bootstrap replicates of experiment
bs_replicates_of_experiment_files <- list.files(bs_replicates_dir, pattern = "\\.csv$")
print(bs_replicates_of_experiment_files)  ## list all files in path

# Create folder to store results
bs_preprocessing_dir <- file.path(paste0(parent_dir, '/ANALYSIS_RESULTS/PRE_PROCESSING_BS_REPLICATES/'))
dir.create(bs_preprocessing_dir)

# Now we need to loop over each of the BS replicates of the experiment to create count
# and meta datafiles
k = 1

for(iFile in bs_replicates_of_experiment_files) {  ## loop through the files
  
  print(paste0('Processing BS replicate file ', iFile))
  
  # Read csv file
  setwd(bs_replicates_dir)
  tot.data <- read.csv(file =  iFile, sep ="\t", stringsAsFactors = FALSE )
  dat.count <- as.matrix(tot.data)
  
  # Create a folder for each BS replicate within the folder we just created
  file_sans_ext <- tools::file_path_sans_ext(iFile)
  bs_out_dir <- paste0(bs_preprocessing_dir, file_sans_ext)
  dir.create(bs_out_dir)
  setwd(bs_out_dir)
  
  # Set first column as data frame index
  rownames(dat.count) <- dat.count[,1]
  
  # Get remaining columns as data, but only those containing normalised counts
  dat.count <- dat.count[,-(1:2)]
  
  # We only want coluns containing 'counts'
  dat.count<- dat.count[ , grepl( "counts" , colnames(dat.count) ) ]
  write.table(dat.count, "count_data.txt", row.name = T, col.name = T, quote = F, sep = "\t")
  
  # Create metadata from sample names
  meta <- c()
  toto <- colnames(dat.count) 
  l <- dim(dat.count)[2]
  coco <- vector(mode = "character", length = l)
  
  for(i in 1:l){
    
    sample_id = toto[i]
    
    temp2 <- str_split(sample_id, "_" )
    temp2 <- temp2[[1]]
    
    polyQ <- temp2[1]
    replicate <- temp2[2]
    
    if (length(temp2) == 6){
      last <- temp2[6]
    }
    else {
      last <- temp2[5]
    }
    
    meta <- rbind(meta, c( sample_id, polyQ, replicate))
  }
  
  rownames(meta) <- meta[, 1]
  meta <- meta[, -1]
  colnames(meta) <- c("polyQ.length", "replicate")
  write.table(meta, "meta_data.txt", row.name = T, col.name = T, quote = F, sep = "\t")

  # Convert character to numeric
  dat.count <- type.convert(dat.count)
  
  # Keep the rows whose counts are larger than 0
  dat <- dat.count[rowSums(dat.count) > 0.0,]
  dat <- dat[rowSums(dat) > 0,]
  print(paste0("Number of rows with counts > 0: ", dim(dat)))
  write.table(dat, "count_dat.txt", row.name = T, col.name = T, quote = F, sep = "\t")
  
  # Move on to next loop
  k = k + 1
  
}


