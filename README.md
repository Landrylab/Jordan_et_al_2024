# Jordan_et_al_2024
## Residues Neighboring an SH3-Binding Motif Participate in Determining Affinity and Specificity In Vivo
All scripts and code for figures for Jordan et al 2024 https://doi.org/10.1101/2024.05.13.593936 


This repo is designed to be used with the data contained in the folowing Dryad repository: https://doi.org/10.5061/dryad.79cnp5j3z.
Both folders from Github and the Data folder from Dryad should be placed in the same folder, the paths in the scripts should work.

The [figures](https://github.com/Landrylab/Jordan_et_al_2024/tree/main/figures) folder contains R Markdown documents, that were used to generate all main and supplementary figures.

The [scripts](https://github.com/Landrylab/Jordan_et_al_2024/tree/main/scripts) folder contains all scripts used in the analysis of the data collected in the manuscript. It is seperated into 3 sections :

- sequencing_analysis contains the scripts involved in processing the sequencing data, and obtaining read counts for every mutant in every sequencing library
- DHFR_PCA_analysis contains scripts that process the data from the DHFR-PCA experiments, including the sequencing data. It also has a script for the DHFR-PCA on solid media and the DHFR-PCA growth curves
- AlphaFold contains the script used to predict structures using AlphaFold-Multimer and the package requirements folder used
