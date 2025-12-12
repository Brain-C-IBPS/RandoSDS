
<img src="https://cdn-icons-png.flaticon.com/512/2103/2103633.png" width="48"> Rando-SDS user guide
==============


The present document is a user guide to the `Rando-SDS` workflow.

### What is Rando-SDS?

`Rando-SDS` is a network analysis software workflow that uses bootstrapping techniques for statistical robustness. It can be divided into two main steps: bootstrap replicate generation and `SDS` (Spectral Decomposition of the Signal) analysis. The output of the workflow are dysregulated networks, e.g. dysregulated gene expression networks. 

The `SDS` (Spectral Decomposition of the Signal) part of the `Rando-SDS` workflow has been adapted from the orignal `SDS` software tool developed by the Brain C Lab. `Rando-SDS` also comprises `python`, `R` and `shell` scripts to generate bootstrap replicates and for other tasks of the workflow. 

#####=================================================================#####

Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International license

#####=============================================================#####

Copyright (C) Maialen Arrieta-Lobo, Christian Neri, and Lucile Megret
Maialen Arrieta (maialen.arrieta_lobo@sorbonne-universite.fr) Christian Neri(christian.neri@inserm.fr) Lucile Mégret(lucile.megret@sorbonne-universite.fr) 2025

#####========================#####

### Rando-SDS workflow explained

1. Starting from read count data (e.g. RNAseq counts) containing several (at least 10) experimental replicates (e.g. 10 mice for which a disease condition was tested), `Rando-SDS` generates a number of bootstrap replicates of the experiment (the number of replicates being set by the user). If the user choses to generate 3 bootstrap replicates of the experiment, the workflow will create 3 "mock" experiments to which the experimental replicates are randomly assigned with replacement. 

2. Next, DESeq2 analysis is performed on read count data of bootstrap replicates and on read count data of the original experiment in order to convert the read counts into Log Fold Change (LFC) values.  

2. The obtained LFC files are used as the input to the `SDS` part of the  workflow, which returns dysregulated networks for the bootstrap replicates of the experiment and for the original experiment.

4. The edges in the dysregulated networks corresponding to the bootstrap replicates are counted across replicates to determine in how many bootstrap replicate networks a given edge appears.

5. Last, using a user-inputted threshold for edge presence, the dysregulated network corresponding to the original experiment is filtered to keep the edges that appear in at least the percentage of bootstrap replicate networks  given by the threshold. That is, if a 80% edge presence threshold is chosen, only edges that are present in 80% of the bootstrap replicate networks will be kept for the final network.




### Rando-SDS folder structure

The downloaded ZIP file contains this `README.md` file, a Jupyter notebook and several folders containing python, R and shella scripts and example datasets. 

* **EXAMPLE_INPUT_DATA**: folder containing an example CSV file of RNAseq read counts. More precisely, the example file corresponds to data from [Lee et al. 2020](https://pubmed.ncbi.nlm.nih.gov/32681824/) and contains 4 Huntington's disease conditions (4 different mutation lengths, namely Q2, Q50, Q111, Q170 and Q175) and 10 replicates per disease condition. 

* **mousenet_V2_filtered**: folder containing an example bionetwork, more precisely the [MouseNet v2](https://www.inetbio.org/mousenet/) filtered to keep interactions whose confidence is larger than 0.4.

* **SCRIPTS**: contains `python`, `R` and `shell` scripts necessary to run the `Jupyter` notebook of the workflow. 

* **RandoSDS_SCRIPTS**: contains 6 python scripts that make up the  core part of the `Rando-SDS` analysis workflow. 

* **`RandoSDS_workflow.ipynb`**: notebook to launch a `Rando-SDS` analysis instance for the example data file and bionetwork provided. 

* **`README.md`**: present user guide. 



### Setting up the Rando-SDS analysis environment

The `Rando-SDS` workflow is developped to partially run on a server with enough RAM. Parts of the workflow are run on the local computer and parts are run on a remote server. 

On the local computer, the user should have uncompressed the downloaded ZIP file, obtaining the **`RandoSDS`** folder where the present file is stored and whose structure is described above. Intermediate files necessary for running the whole analysis are generated and stored in the **`RandoSDS`** folder, as well as network analysis result files that are copied from the remote server after the analysis is therein performed. 

On the remote server, the user should create a folder named **`RandoSDS`** and copy two things inside:

* The **mousenet_V2_filtered** folder with all of its content
* The contents of the **RandoSDS_SCRIPTS** folder (6 python scripts)

### Running the Rando-SDS analysis through the Jupyter Notebook 

The **`RandoSDS_workflow.ipynb`** Jupyter Notebook can be followed in order to launch an `Rando-SDS` analysis. The Notebook includes step-by-step explanations so that the user can follow the workflow. 


### License

**Copyright (C) 2023 Maialen Arrieta-Lobo and Christian Neri**

This work is licensed under the Creative CommonsAttribution-NonCommercial-NoDerivatives 4.0  International License. To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-nd/4.0/ or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA. 


### Command to convert md file to html
pandoc -f markdown -t html README.md > README.html
