#-----------------------------------------------------------------------------
#
#               R script of the Rando-SDS v1.0 workflow to
#       perform DESeq2 analysis and DEG extraction for BS replicates
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
# |  Copyright (C) 2025 Maialen Arrieta-Lobo, Lucile Mégret and Christian Neri | 
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

#-- Load R packages --#
library("DESeq2")
#library(flashClust)
library("ggplot2")

# Get script directory
this.dir <- dirname(sys.frame(1)$ofile)

# Path to parent directory
parent_dir <- dirname(this.dir)

# Create folder to store results
pre_processing_dir <- file.path(paste0(parent_dir, '/ANALYSIS_RESULTS/PRE_PROCESSING_BS_REPLICATES'))
print(pre_processing_dir)

# List of polyQ lengths in experiment
polyQ <- c("Q20", "Q50", "Q111", "Q170", "Q175")

# We need to iterate over each of the Boostrap replicate folders containing count and meta files 
bs_replicates_folders <- list.dirs(pre_processing_dir)
bs_replicates_folders <- bs_replicates_folders[ grepl("BSreplicate_", bs_replicates_folders) ]
print(bs_replicates_folders)

for(folder in bs_replicates_folders) {  ## loop through the bs replicate folders
  
  bs_folder <- tools::file_path_sans_ext(folder) 
  print(paste0('Starting evaluation of bootstrap:', folder  ))
  
  # Set data directory to bs folder
  work.dir = folder
  setwd(work.dir)
  
  # Get count file
  dat <- as.matrix(read.table("count_data.txt",header = T))
  
  # Convert character to numeric
  dat <- type.convert(dat)
  
  # Get metadata file
  lab <- as.matrix(read.table("meta_data.txt",header = T))
  print( paste0(  'Metadata matrix dimension = ', dim(lab)))
  
  # Get corresponding counts
  dat <- dat[, rownames(lab)]
  print( paste0( '6 months Dat matrix dimension = ', dim(dat)))
  
  # -- Define output directory --#
  out.dirB <- paste0(work.dir, "/DESeq2_RESULTS") 
  dir.create(out.dirB)
  setwd(out.dirB)
  
  #--------------------------------------------------------------------------------------------------------------------#
  
  filtering <- TRUE
  
  #-- remove genes with less than 10 read counts in all samples --#
  if (filtering == TRUE) dat <- dat[ rowSums(dat < 5) != 0.8*ncol(dat), ]
  
  for(i in 1:1){  
    out.dir <- paste0(out.dirB, "/")
    dir.create(out.dir)
    setwd(out.dir)
    
    # Loop over disease condition polyQ length  
    for( j in 2:5 ) {
      
      print(paste0( polyQ[j], 'loop'))
      
      #cat(polyQ[j],"-","\n", sep=" ")
      
      #-- Sample.ID: Control condition -> "Q20" --# 
      sublab1 <- lab[ lab[, "polyQ.length"] == "Q20" , ]
      
      #-- Sample.ID: HET condition -> polyQx --#
      sublab2 <- lab[ lab[, "polyQ.length"] == polyQ[j] , ]
      
      # Create metadata for control and current polyQ length
      sublab <- rbind(sublab1, sublab2)
      
      # Get row indices corresponding to selected samples
      ix <- rownames(sublab)
      
      #-- Extract 'countData' (matrix) and 'colData' (data.frame) --#
      countData <- as.matrix( dat[, ix] )  
      countData <- type.convert(countData)
      colData <- sublab
      colData <- as.data.frame(colData)
      
      #--------------------------------------------------------------------------------------------
      # The formula should be a tilde (~) followed by the variables with plus signs between them.
      # Note: In order to benefit from the default settings of the package, 
      # you should put the variable of interest at the end of the formula 
      # and make sure the control level is the first level.
      #--------------------------------------------------------------------------------------------
      
      #-- Truncate "read count decimal values" as non-negative integer values --#
      mode(countData) <- "integer"
      
      #-- Define a DESeq dataset object --#
      dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~polyQ.length)
      dds$polyQ.length <- relevel(dds$polyQ.length, "Q20") 
      
      png(filename = paste(out.dir, paste("Boxplots_pseudoreadcounts_",polyQ[j],"_",".png",sep=""), sep=""), width = 2000, height = 2500, units = "px", pointsize = 12, bg = "white",  res=300)
      names <- paste(rownames(colData),  colData[,"polyQ.length"], sep="_")
      par(mar=c(2+round(max(nchar(names))/2),4,2,1))
      title <- paste("pseudo read counts - ",polyQ[j],"m - ",sep="")
      bxp <- boxplot(log2(countData)+1, boxwex=0.7, notch=F, main=title, outline=FALSE, col = c("gray90","gray60")[as.numeric(dds$polyQ.length)], names=names, las=2)
      dev.off()
      
      #---------------------------------------------------------------------------------------------------------------------------------------
      # Results tables are generated using the function results, which extracts a results table with log2 FC, p values and adjusted p values. 
      # With no arguments to results, the results will be for the last variable in the design formula, 
      # and if this is a factor, the comparison will be the last level of this variable over the first level.
      # Details about the comparison are printed to the console. 
      # The text, condition treated vs untreated, tells you that the estimates are of the lFC log2(treated/untreated).
      #---------------------------------------------------------------------------------------------------------------------------------------
      
      #-- Perform differential expression analysis -- #
      
      # There are two reasons for so many genes being flagged as outlier: 
      # either the method for flagging outliers is not appropriate for the distribution of counts in your data 
      # and should be turned off (by setting minReplicatesForReplace=Inf and cooksCutoff=FALSE),
      
      dds <- DESeq(dds, test = "Wald", betaPrior=TRUE, minReplicatesForReplace = Inf)
      res <- results(dds, independentFiltering = FALSE, cooksCutoff = FALSE, pAdjustMethod = "fdr")
      res$padj[ is.na(res$padj) ] <- 1 
      res$pvalue[ is.na(res$pvalue) ] <- 1 
      res$log2FoldChange[ is.na(res$log2FoldChange) ] <- 0 
      
      #-- MA-plot --#
      png(filename = paste(out.dir, paste("MAPlot_DESeq2_",polyQ[j],"m_",".png",sep=""), sep=""), width = 2000, height = 2500, units = "px", pointsize = 12, bg = "white",  res=300)
      
      DESeq2::plotMA(dds, main=paste("DESeq2 - ",polyQ[j],"m - ",sep=""))
      
      dev.off()
      
      #-- Retrieve lFC and drg files --#
      lFC <- as.matrix(res$log2FoldChange)
      rownames(lFC) <- rownames(res) 
      
      png(filename = paste(out.dir, paste("VolcanoPlot_DESeq2_",polyQ[j],".png",sep=""), sep=""), width = 2000, height = 2500, units = "px", pointsize = 12, bg = "white",  res=300)
      par(las=1)
      
      if ( length(unique(res$padj)) > 200 ) { mycol <- densCols(res$log2FoldChange,-log10(res$padj),colramp = colorRampPalette(c("grey40", "grey0"))) 
      } else { mycol <- rep("grey40", length(res$padj)) }
      mycol[ res$padj < .05 & res$log2FoldChange > 0 ] <- "red"
      mycol[ res$padj < .05 & res$log2FoldChange < 0 ] <- "green4"
      
      plot(res$log2FoldChange,
           -log10(res$padj), 
           col="white", cex=.5, xaxt="n",yaxt="n", 
           xlab="log2(FC)", ylab="-log10(FDR)", 
           ylim=c(0,5), xlim=c(-1,1), 
           main=paste("DESeq2 - ",polyQ[j],"m - ",sep=""),
           font.main=1)
      
      #with(subset(res, rownames(res) %in% gene_list), points(log2FoldChange, -log10(pvalue), pch=20, col="orange"))
      
      axis(1,at=seq(-2,2,.5),seq(-2,2,.5))
      axis(2,at=seq(-0.5,7.5,.5),seq(-0.5,7.5,.5))
      abline(h=seq(-3,8,.5), lty=4, col="gray60")
      abline(v=seq(-60,60,.5), lty=4, col="gray60")
      abline(h=-log10(.05),lty=5,col="darkblue")
      points(res$log2FoldChange,-log10(res$padj), col=mycol, cex=.5, pch=16)
      with(subset(res, rownames(res) %in% c('ENSMUSG00000002324')), points(log2FoldChange,-log10(padj), pch=20, col="black"))
      legend("topleft", legend=c(paste("up = ",sum(res$padj < .05 & res$log2FoldChange > 0),sep=""),paste("down = ",sum(res$padj < .05 & res$log2FoldChange < 0),sep="")),
             pch=rep(16,2),col=c("red","green4"),bty="o",bg="mintcream")
      box(which="plot")
      
      dev.off()
      
      df <- data.frame(res)
      
      df$shape <- ifelse(df$padj < .05, "triangle", "circle")
      
      df$col <- df$shape
      df$col[ df$padj < .05 & df$log2FoldChange > 0 ] <- "triangle_U"
      df$col[ df$padj < .05 & df$log2FoldChange < 0 ] <- "triangle_D"
      df$col[ df$padj > .05 ] <- "circle"
      df$padj[ df$padj < .001] <- .001
      
      df$shape[(abs(df$log2FoldChange) > 1.5)] <- "triangle"
      df$log2FoldChange[df$log2FoldChange >  1.5] <- 1.5
      df$log2FoldChange[df$log2FoldChange < -1.5] <- -1.5
      
      highlight_df <- subset(df, rownames(df) %in% c('ENSMUSG00000069833'))
      highlight_df$shape <- "square"
      highlight_df$col <- "square"
      
      png(filename = paste(out.dir, paste("VolcanoGGPlot_DESeq2_",polyQ[j],"m_",".png",sep=""), sep=""), width = 2000, height = 2500, units = "px", pointsize = 12, bg = "white",  res=300)
      myplot <- ggplot(data=df, aes(x=log2FoldChange, y=-log10(padj), shape=shape)) +
        geom_point(alpha=0.5, aes(color=factor(df$col)), size=3) +
        geom_point(data=highlight_df, aes(x=log2FoldChange, y=-log10(padj), shape=shape)) + 
        scale_color_manual(values = c("circle" = "gray50", "triangle_D" = "green4", "triangle_U" = "red", "square" = "blue")) +
        scale_size_manual(values = c("circle" = 3, "triangle_D" = 3, "triangle_U" = 3, "square" = 10)) +
        theme(legend.position = "none") +
        xlim(c(-1.2, 1.2)) + ylim(c(0, -log10(0.001))) +
        xlab("log2 fold change") + ylab("-log10 p-value") + ggtitle(paste("DESeq2 - ",polyQ[j],"m - ", sep="")) +
        annotate("point",c(.80,.80),c(2.67,2.80), size=3,col=c("green4","red"), shape=17) +
        annotate("text",c(.85,.85),c(2.685,2.815),
                 label=c(paste("down = ",sum(res$padj < .05 & res$log2FoldChange < 0),sep=""),paste("up = ",sum(res$padj < .05 & res$log2FoldChange > 0),sep="")),
                 hjust=-0.1)  
        
        
      print(myplot)
      dev.off()
      
      drg <- rep(0, nrow(res))
      drg[ res$pvalue < .05 ] <- 1
      drg <- as.matrix(cbind(rownames(res),drg))
      out.file <- paste(polyQ[j],".txt",sep="")
      
      #-- Save lFC and dereg files --#
      
      write.table(lFC, paste(out.dir,"lfc_",out.file,sep=""), row.names=T, col.names=F, quote=F, sep="\t")
      write.table(drg, paste(out.dir,"drg_",out.file,sep=""), row.names=F, col.names=F, quote=F, sep="\t")
  
    }} # end_i_j
  

  # Loop over cell types  
  for(i in 1:1){  
    print('Looping')
    LFC.cible <- matrix(data = NA,nrow = nrow(dat) ,ncol = 4 )
    Pval.cible <- matrix(data = NA,nrow = nrow(dat) ,ncol = 4 )
    rownames(Pval.cible) <- rownames(dat)
    rownames(LFC.cible) <- rownames(dat)
    colna <- vector(mode = "character", length = 4)
    
    u=1
    for(j in 2:5){

      out.file <- paste(polyQ[j],".txt",sep="")
      lfc <- as.matrix(read.table(paste(out.dir,"lfc_",out.file,sep=""),header = F))
      pval <- as.matrix(read.table( paste(out.dir,"drg_",out.file,sep=""),header = F))
      colna[u]= paste0("at_",polyQ[j])
      for(s in 1:nrow(dat) ){
        LFC.cible[s,u] <- as.numeric(lfc[ which(lfc[,1]== rownames(dat)[s]),2]) 
        Pval.cible[s,u] <- pval[ which(pval[,1]== rownames(dat)[s]),2]
      }
      u <- u+1 
    }
    colnames(LFC.cible) <- colna
    colnames(Pval.cible) <- colna
    
    write.table(LFC.cible,paste0("Table_LFC.txt")  ,row.names = T,col.names = T, quote = F,sep = "\t")
    write.table(Pval.cible, paste0("Table_Pval.txt") ,row.names = T,col.names = T, quote = F,sep = "\t")
    
    
  }
  
  setwd(work.dir)
  print(work.dir)
  print(paste0('Ended evaluation of bootstrap:', bs_folder  ))
  
} # end bs_replicate folder loop
  
  
