#!/bin/bash

#Originally obtained from the website of the Digital Research Alliance of Canada and run on the Alliance's Graham cluster. https://docs.alliancecan.ca/wiki/AlphaFold . The version below is a (slight) modification, and removed the information for running this using a SLURM system. You can add whatever SLURM parameters are necessary for your cluster. The results in the paper were run on CPU and GPU. 
#Make sure that the alphafold-requirements.txt file is in your home (otherwise change the path).
#The INPUT_DIR/$FILE should be a fasta file with the complete sequence of Sho1 and  of the Pbs2 mutant you are modeling between positions 71 and 126. Here is an example for I87W:
#>Sho1
#MSISSKIRPTPRKPSRMATDHSFKMKKFYADPFAISSISLAIVSWVIAIGGSISSASTNESFPRFTWWGIVYQFLIICSLMLFYCFDLVDHYRIFITTSIAVAFVYNTNSATNL#VYADGPKKAAASAGVILLSIINLIWILYYGGDNASPTNRWIDSFSIKGIRPSPLENSLHRARRRGNRNTTPYQNNVYNDAIRDSGYATQFDGYPQQQPSHTNYVSSTALAGFEN#TQPNTSEAVNLHLNTLQQRINSASNAKETNDNSNNQTNTNIGNTFDTDFSNGNTETTMGDTLGLYSDIGDDNFIYKAKALYPYDADDDDAYEISFEQNEILQVSDIEGRWWKAR#RANGETGIIPSNYVQLIDGPEEMHR
#>Pbs2_I87W
#SNQSEQDKGSSQSPKHWQQIVNKPLPPLPVAGSSKVSQRMSSQVVQASSKSTLKNV  
#                ^ note I87W mutation 


# Load modules dependencies.
module load gcc/9.3.0 openmpi/4.0.3 cuda/11.4 cudnn/8.2.0 kalign/2.03 hmmer/3.2.1 openmm-alphafold/7.5.1 hh-suite/3.3.0 python/3.8

DOWNLOAD_DIR=/datashare/alphafold   # set the appropriate path to your downloaded data
INPUT_DIR=$SCRATCH/alphafold/input     # set the appropriate path to your input data
OUTPUT_DIR=${SCRATCH}/alphafold/output # set the appropriate path to your output data

# Generate your virtual environment in $SLURM_TMPDIR.
virtualenv --no-download ${SLURM_TMPDIR}/env
source ${SLURM_TMPDIR}/env/bin/activate

# Install AlphaFold and its dependencies.
pip install --no-index --upgrade pip
pip install --no-index --requirement ~/alphafold-requirements.txt

# Edit with the proper arguments and run your commands.
run_alphafold.py \
   --fasta_paths=${INPUT_DIR}/$FILE \
   --output_dir=${OUTPUT_DIR} \
   --data_dir=${DOWNLOAD_DIR} \
   --db_preset=full_dbs \
   --model_preset=multimer \
   --bfd_database_path=${DOWNLOAD_DIR}/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt \
   --mgnify_database_path=${DOWNLOAD_DIR}/mgnify/mgy_clusters_2022_05.fa \
   --template_mmcif_dir=${DOWNLOAD_DIR}/pdb_mmcif/mmcif_files \
   --obsolete_pdbs_path=${DOWNLOAD_DIR}/pdb_mmcif/obsolete.dat \
   --pdb_seqres_database_path=${DOWNLOAD_DIR}/pdb_seqres/pdb_seqres.txt \
   --uniprot_database_path=${DOWNLOAD_DIR}/uniprot/uniprot.fasta \
   --uniref30_database_path=${DOWNLOAD_DIR}/uniref30/UniRef30_2021_03 \
   --uniref90_database_path=${DOWNLOAD_DIR}/uniref90/uniref90.fasta \
   --hhblits_binary_path=${EBROOTHHMINSUITE}/bin/hhblits \
   --hhsearch_binary_path=${EBROOTHHMINSUITE}/bin/hhsearch \
   --jackhmmer_binary_path=${EBROOTHMMER}/bin/jackhmmer \
   --kalign_binary_path=${EBROOTKALIGN}/bin/kalign \
   --max_template_date=2024-02-22 \
   --use_gpu_relax=True

