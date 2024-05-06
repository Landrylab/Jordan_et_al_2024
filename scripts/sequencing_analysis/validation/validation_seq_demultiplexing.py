# -*- coding: utf-8 -*-
"""
Created on Thu Apr 11 10:46:23 2024

@author: DAJOR4
"""

#From the raw reads, demultiplex into replicates of the DHFR-PCA screen based on barcodes. The results of this script are available in the Dryad repository

#%% Importing packages
import os

import subprocess

import pandas as pd
print(pd.__name__, pd.__version__)

import numpy as np
print(np.__name__, np.__version__)

import matplotlib.pyplot as plt
import matplotlib
print(matplotlib.__name__, matplotlib.__version__)

import scipy
print(scipy.__name__, scipy.__version__)

import seaborn as sns
print(sns.__name__, sns.__version__)

#%% Set working folder

#Set the folder containing the script, should contain script and RC_pcr_for.fa and RC_pcr_rev.fa files. common_sequences folder should be in the folder above the one containing the script
os.chdir("/your/path")

print(os.listdir())

empty_datasheet="../../../Data/demultiplexed_sequencing_reads/validation/validation_sample_datasheet_empty.csv"
demultiplexed_datasheet='../../../Data/demultiplexed_sequencing_reads/validation/validation_sample_datasheet_demultiplexed_merged.csv'

#%% Name stuff
experiment_name = 'Pbs2_validation'

#Path to the raw sequencing files
path_to_R1_file = '/path/to/R1/file.fastq.gz'

path_to_R2_file = '/path/to/R2/file.fastq.gz'



#%% Run FastQc on reads

fastqc_call= "fastqc "+ path_to_R1_file + " " + path_to_R2_file


print(fastqc_call)

subprocess.check_output(fastqc_call, shell=True)


#%%Filter short reads

if "temp" not in os.listdir("./"):
    os.mkdir("./temp")

#Here the path will have to be changed to your installation of trimmomatic
trim_call = 'java -jar /home/adminlandry/anaconda3/share/trimmomatic-0.39-2/trimmomatic.jar PE -threads 6 -trimlog log.txt '
trim_call += path_to_R1_file+' '+path_to_R2_file+' '+'-baseout ./temp/minlen298.fastq MINLEN:298'
#MinLEN at 298

print(trim_call)

subprocess.check_output(trim_call, shell=True)



#%% Demultiplexing with cutadapt

if "demultiplex" not in os.listdir("./"):
    os.mkdir("./demultiplex")

# Try with trimming 36 bp from the beginning of each read. This should begin at the RC barcodes

trimmed_forward= "./temp/minlen298_1P.fastq"

trimmed_reverse="./temp/minlen298_2P.fastq"

cutadapt_call="cutadapt -u 36 -U 36 -j 0 -e 0 --discard-untrimmed -g file:../../../Data/demultiplexed_sequencing_reads/validation/RC_pcr_for.fa -G file:../../../Data/demultiplexed_sequencing_reads/validation/RC_pcr_rev.fa -o ./demultiplex/{name1}-{name2}.for.fastq.gz -p ./demultiplex/{name1}-{name2}.rev.fastq.gz "+ trimmed_forward+ " " +trimmed_reverse

print(cutadapt_call)

subprocess.check_output(cutadapt_call, shell=True)

#%% Import library information datasheet
sample_dataframe=pd.read_csv(empty_datasheet, index_col=0)


#%% Update sample datasheet 

for sample_ID in list(sample_dataframe.index):
    sample_dataframe.loc[sample_ID, 'filepath']="./demultiplex/" + str(sample_dataframe.loc[sample_ID, 'RC_for_index'])+"-"+str(sample_dataframe.loc[sample_ID, 'RC_rev_index'])
# adds the demultiplexed filepaths to the sample DataFrame

#%% Unzip demultiplexed files

unzip_call="gunzip -f ./demultiplex/*.fastq.gz"



subprocess.check_output(unzip_call, shell =True)


#%%Merge reads with pandaseq

if "merged_reads" not in os.listdir("./"):
    os.mkdir("./merged_reads")
    



before_merge_depth_dict = {}
# empty container that wil hold the initial number of reads for each sample

def merge_reads(Sample_ID):
    # this function generates and runs a Pandaseq call to merge R1 and R2 for each demultiplexed sample.
    # it returns the starting number of reads in the file, which is also stored in a dict
    
    filepath = sample_dataframe.loc[Sample_ID]['filepath']
    filepath_for = filepath + ".for.fastq"
    filepath_rev = filepath + ".rev.fastq"
    # define the R1 and R2 file paths
    
    sample_read_count = 0
    # counter, set to 0
    
    with open(filepath_for, 'r') as source:
        for line in source:
            if line.startswith('@M'):
                sample_read_count += 1
                # opens file and loops thourgh it while counting the number of headers to count the initial number of reads
                
    before_merge_depth_dict[Sample_ID] = sample_read_count
    # adds inital read count to dict
                
    if sample_read_count >= 1:
        # for samples  with reads
                
        filepath_out = './merged_reads/'+str(Sample_ID)+'.fasta'
        # define the merging output filepath
        
        panda_seq_call = 'pandaseq -f '+filepath_for+' -r '+filepath_rev+ ' -L 340 -O 300 -g "./merged_reads/log.txt" -k 4 -B -N -t 0.5 -T 6 -w '+ filepath_out
        subprocess.check_output(panda_seq_call, shell=True)
        # generates and runs a pandseq call with the following parameters:
        #    -L 340             maximum length of assembled sequence: 340 pb
        #    -O 300             maximum overlap between reads: 300 pb
        #    -k 4               number of sequence locations per kmer (increased from default value of 2)
        #    -B                 allow input to lack a barcode
        #    -N                 eliminate all sequences with Ns
        #    -t 0.5             score threshold
        #    -T 6               use 6 threads
        #    -g  ./merged_reads/log.txt log in file instead of in standard error
        #    -w filepath_out    write output to filepath_out
        
    return sample_read_count
    # return initial number of reads
    
    

#%% Merge all libraries

for sample_ID in list(sample_dataframe.index):
    # for sample in DataFrame
    print(sample_ID)
    merge_reads(sample_ID)
    # merge reads

sample_dataframe["demult_reads"]=before_merge_depth_dict


#%% Count merged reads

after_merge_depth_dict = {}
# empty container that wil hold the post-merge number of reads for each sample

def count_merged(Sample_ID):
    
    filepath = './merged_reads/'
    filepath_merged = str(Sample_ID)+'.fasta'
    filepath+=filepath_merged
    # find merged reads fasta file path based on sample name
    
    after_merge_depth = 0
    # counter for read number
      
    with open(filepath, 'r') as source:
        for line in source:
            if line.startswith('>'):
                after_merge_depth += 1
                # open merged read file, loop through lines,
                # find the headers and increment by one for each header 

    after_merge_depth_dict[Sample_ID] = after_merge_depth
    # once file has been read through, add read number to the dictionary
    return after_merge_depth
    # return merged read counts
    
for sample_ID in list(sample_dataframe.index):
    # for each sample
    count_merged(sample_ID)
    # count the number of reads after merging
    
sample_dataframe['merged_read_count'] = pd.Series(after_merge_depth_dict)
sample_dataframe['merged_read_count']


plt.bar([x for x in range(0,78)], sample_dataframe['merged_read_count'])
plt.xlabel('Sample', fontsize=14)
plt.ylabel('Number of merged reads', fontsize=14)

#%% Export sample_datasheet

sample_dataframe.to_csv(demultiplexed_datasheet, sep=',')
