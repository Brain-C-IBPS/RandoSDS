# Rando-SDS

**Network analysis workflow using bootstrapping techniques for statistical robustness**

Rando-SDS is a network analysis software workflow that combines bootstrap replicates generation with Spectral Decomposition of the Signal (SDS) analysis to identify dysregulated networks, such as dysregulated gene expression networks.


#####=================================================================#####
#####       Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International license          #####
#####=============================================================#####
#####               Copyright (C) Maialen Arrieta-Lobo, Francesca Farina, Tamara Monteagudo Aboy, Megan Mair, Cloé Mendoza, Huy Tran, Jeff Aaronson, Jim Rosinski, Lisa Ellerby, Emmanuel Brouillet, Juan Botas, Christian Neri, and Lucile Megret #####
#####                     Maialen Arrieta (maialen.arrieta_lobo@sorbonne-universite.fr) Christian Neri(christian.neri@inserm.fr) Lucile Mégret(lucile.megret@sorbonne-universite.fr) 2024                               #####
#####========================#####
#      
#      This work is licensed under the Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 
#      International License. To view a copy of this license, visit 
#      http://creativecommons.org/licenses/by-nc-nd/4.0/ or send a letter to 
#      Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.
#      
#####======================#####
#####=====================#####


## Overview

The workflow uses bootstrapping techniques to enhance the statistical robustness of network analysis. Starting from read count data (e.g., RNA-seq), Rando-SDS generates bootstrap replicates, performs differential expression analysis, identifies dysregulated networks, and filters results based on edge presence across replicates.

### Key Features

- Bootstrap replicate generation from experimental data
- Integration with DESeq2 for differential expression analysis
- SDS-based network dysregulation detection
- Statistical validation through bootstrap edge counting
- Customizable edge presence thresholds

## Workflow Steps

1. **Bootstrap Generation**: Creates multiple bootstrap replicates from experimental data by randomly assigning samples with replacement
2. **Differential Expression**: Performs DESeq2 analysis on both bootstrap replicates and original data to obtain Log Fold Change (LFC) values
3. **Network Analysis**: Applies SDS methodology to identify dysregulated networks for each bootstrap replicate and the original experiment
4. **Edge Counting**: Counts edge appearances across bootstrap replicate networks
5. **Network Filtering**: Filters the original dysregulated network based on user-defined edge presence threshold (e.g., keep edges present in ≥80% of bootstrap networks)

## Repository Structure

```
RandoSDS/
├── EXAMPLE_INPUT_DATA/          # Example RNA-seq read count data
├── mousenet_V2_filtered/        # Example bionetwork (MouseNet v2)
├── SCRIPTS/                     # Python, R, and shell scripts
├── RandoSDS_SCRIPTS/            # Core Rando-SDS analysis scripts (6 Python scripts)
├── RandoSDS_workflow.ipynb      # Main workflow Jupyter notebook
└── README.md                    # User guide
```

## Installation

### Prerequisites

- Python (version TBD)
- R (version TBD)
- Jupyter Notebook
- DESeq2 (R package)
- Access to a server with sufficient RAM for network analysis

### Setup

1. Clone this repository:
   ```bash
   git clone https://github.com/yourusername/RandoSDS.git
   cd RandoSDS
   ```

2. Install dependencies:
   ```bash
   # Python dependencies
   pip install -r requirements.txt
   
   # R dependencies
   # Install DESeq2 and other required packages (details TBD)
   ```

3. Set up the remote server:
   - Create a folder named `RandoSDS` on your server
   - Copy the `mousenet_V2_filtered/` folder to the server
   - Copy all scripts from `RandoSDS_SCRIPTS/` to the server

## Usage

The analysis is performed through the Jupyter notebook `RandoSDS_workflow.ipynb`, which provides step-by-step instructions for:

1. Loading your RNA-seq read count data
2. Configuring bootstrap parameters (number of replicates, edge presence threshold)
3. Running the analysis pipeline
4. Interpreting results

### Quick Start

```bash
jupyter notebook RandoSDS_workflow.ipynb
```

Follow the notebook cells sequentially. The workflow is designed to run partially on your local machine and partially on a remote server.

## Example Data

The repository includes example data from [Lee et al. 2020](https://pubmed.ncbi.nlm.nih.gov/32681824/):
- RNA-seq read counts from Huntington's disease conditions
- 4 mutation lengths (Q2, Q50, Q111, Q170, Q175)
- 10 replicates per condition

The example bionetwork is [MouseNet v2](https://www.inetbio.org/mousenet/) filtered for interactions with confidence >0.4.

## Input Data Format

Your input data should be a CSV file containing:
- RNA-seq read counts
- Multiple experimental replicates (minimum 10 recommended)
- Appropriate column headers for samples and conditions

## Output

The workflow generates:
- Bootstrap replicate networks
- Edge count statistics across replicates
- Filtered dysregulated network based on edge presence threshold
- Visualization and summary statistics

## System Requirements

- **Local machine**: Standard desktop/laptop with Jupyter notebook capability
- **Remote server**: Server with sufficient RAM for network analysis (specific requirements TBD)

## Citation

If you use Rando-SDS in your research, please cite:

[Citation information to be added]

### Example Data Citation

If using the provided example data, please cite:
Lee H, et al. (2020). [Publication details from PMID: 32681824]

## Acknowledgments

The SDS (Spectral Decomposition of the Signal) component of Rando-SDS is adapted from the original SDS software tool developed by the Brain C Lab.

[Additional attribution details to be added]

## License

Copyright (C) 2025 Maialen Arrieta-Lobo, Lucile Megret and Christian Neri

This work is licensed under the Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International License. To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-nd/4.0/

## Support

For questions, issues, or feedback:
- Open an issue on GitHub
- [Contact information to be added]

## Contributing

[Contribution guidelines to be added]

## Changelog

See [CHANGELOG.md](CHANGELOG.md) for version history and updates.

---

**Note**: This repository is under active development. Some documentation sections are marked as "TBD" and will be updated as the project matures.
