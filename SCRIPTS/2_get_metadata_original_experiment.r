#-----------------------------------------------------------------------------
#
#            R script of the Rando-SDS v1.0 workflow to
#            preprocess original experiment read count data
#
#-----------------------------------------------------------------------------
#
# Contributors (alphabetic order): Maialen Arrieta-Lobo (1) 
#                                  Lucile Mégret (1)
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
# |  Copyright (C) 2023 Maialen Arrieta-Lobo, Lucile Mégret and Christian Neri | 
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

# Get CSV files of bootstrap replicates of experiment
experimental_data_file <-  file.path(paste0(parent_dir, "/EXAMPLE_INPUT_DATA/example_read_count_file.csv"))

# Create folder to store results
original_experiment_preprocessing_dir <- file.path(paste0(parent_dir, '/ANALYSIS_RESULTS/PRE_PROCESSING_ORIGINAL_EXPERIMENT/'))
dir.create(original_experiment_preprocessing_dir)
setwd(original_experiment_preprocessing_dir)

# Read csv file
tot.data <- read.csv(file = experimental_data_file, sep ="\t", stringsAsFactors = FALSE )
dat.count <- as.matrix(tot.data)
View(head(dat.count))

# Set first column as data frame index
rownames(dat.count) <- dat.count[,1]

# Get remaining columns as data, but only those containing normalised counts
dat.count <- dat.count[,-(1:2)]

# We only want coluns containing 'norm'
dat.count<- dat.count[ , grepl( "counts" , colnames(dat.count) ) ]
write.table(dat.count, "count_data.txt", row.name = T, col.name = T, quote = F, sep = "\t")

#meta <- matrix(data = NA,nrow = dim(dat.count)[2],ncol = 2)

# Create metadata from sample names
meta <- c()
toto <- colnames(dat.count) 
print(toto)
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
View(meta)
write.table(meta, "meta_data.txt", row.name = T, col.name = T, quote = F, sep = "\t")

# Convert character to numeric
dat.count <- type.convert(dat.count)

# Keep the rows whose counts are larger than 0
dat <- dat.count[rowSums(dat.count) > 0.0,]
View(head( dat  ))
print(paste0("Number of rows with counts > 0: ", dim(dat)))

write.table(dat, "count_dat_norm_All.txt", row.name = T, col.name = T, quote = F, sep = "\t")

dat <- dat[rowSums(dat) > 0,]
print(paste0("Number of rows with counts > 0: ", dim(dat)))

write.table(meta, "meta_dat_norm.txt", row.name = T, col.name = T, quote = F, sep = "\t")
write.table(dat, "count_dat_norm.txt", row.name = T, col.name = T, quote = F, sep = "\t")



