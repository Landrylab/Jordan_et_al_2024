# -*- coding: utf-8 -*-
"""
Created on Thu Apr 11 14:09:09 2024

@author: DAJOR4
"""

#From the raw reads, demultiplex into replicates of the DHFR-PCA screen based on barcodes. The results of this script are available in the Dryad repository

#%% Chunk 1 : Importing packages
import os

import subprocess

import pandas as pd
print(pd.__name__, pd.__version__)

import numpy as np
print(np.__name__, np.__version__)

import matplotlib.pyplot as plt
import matplotlib
print(matplotlib.__name__, matplotlib.__version__)

from collections import Counter
#import regex as regex

import seaborn as sns


#%% Chunk 2 : 

#Set this to the folder containing the script.  The common_sequences folder should be in the folder above this one.
os.chdir("/your/path/")
print(os.listdir())
universal_seqs_fasta = "../../../Data/demultiplexed_sequencing_reads/common_sequences/universal_amplicon.fa"
indexes_fasta = "../../../Data/demultiplexed_sequencing_reads/common_sequences/indexes.fa"
empty_datasheet="../../../Data/demultiplexed_sequencing_reads/surrounding_region/surrounding_region_sample_datasheet_empty.csv"
demultiplexed_datasheet='../../../Data/demultiplexed_sequencing_reads/surrounding_region/surrounding_region_sample_datasheet_demultiplexed_merged.csv'

#%% Chunk 3 : Basic functions and data import
def reverse_complement(dna):
    """ function that reverse complements DNA
    dna: input dna sequence
    """
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join([complement[base] for base in dna[::-1]])

def get_dict_of_seq(fasta_file):
    """ function that converts a fasta file to a dictionary of sequences
    fasta_file: the input fasta file
    """
    
    file_fasta_dict = {}
    # output dict of imported seqs
    
    with open(fasta_file, 'r') as fasta:    
        for line in fasta:
            # loops through the file

            if line.startswith('>') == True:
                seq_info = line.strip('>').strip('\n').split('\t')[0]
                file_fasta_dict[seq_info] = ''
                # checks if seq header

            else:
                file_fasta_dict[seq_info] += line.strip('\n')
                # If not, append nt to current seq
                
    return file_fasta_dict



universal_seqs = get_dict_of_seq(universal_seqs_fasta)
print(universal_seqs)


indexes = get_dict_of_seq(indexes_fasta)
print(indexes)

degen_seq = 'NNNNN'

del(universal_seqs_fasta, indexes_fasta)

#%% Name stuff
experiment_name = 'Pbs2_surrounding_region'

#Change to the path to the raw sequencing files
path_to_R1_file = './path/to/R1/file.fastq.gz'
path_to_R2_file = './path/to/R2/file.fastq.gz'


#%%Run FastQC on reads
subprocess.check_output('fastqc '+path_to_R1_file, shell=True)#  path to the R1 (forward) compressed sequencing file
subprocess.check_output('fastqc '+path_to_R2_file, shell=True)# path to the R2 (reverse) compressed sequencing file



#%% Trim reads with TRIMMOMATIC


#Change the path to your installation of trimmomatic
trim_call = 'java -jar /home/adminlandry/anaconda3/share/trimmomatic-0.39-2/trimmomatic.jar PE -threads 6 -trimlog log.txt '
trim_call += path_to_R1_file+' '+path_to_R2_file+' '+'-baseout ./temp/minlen250.fastq MINLEN:250'
#Only keep reads above 250 bp
subprocess.check_output(trim_call, shell=True)


#%% Check quality with FastQC after trimming
subprocess.check_output('fastqc '+'./temp/minlen250_1P.fastq', shell=True)
subprocess.check_output('fastqc '+'./temp/minlen250_2P.fastq', shell=True)



#%%Make folders for next steps
if 'amplicon_sequences' not in os.listdir('./'):
    os.mkdir('./amplicon_sequences')
    
if 'bowtie_indexes' not in os.listdir('./'):
    os.mkdir('.bowtie_indexes')
    
#%%Check abundance of libraries

plate_pcr_amplicon_for = './amplicon_sequences/plate_pcr_for.fa'
plate_pcr_amplicon_rev = './amplicon_sequences/plate_pcr_rev.fa'

for_plate_indexes = [21,1]

rev_plate_indexes = [21,1]

print(for_plate_indexes)
print(rev_plate_indexes)


with open(plate_pcr_amplicon_for, 'w') as fasta_out:
    
   for index in for_plate_indexes:

        index_seq = indexes[str(index)]
        
        header = '>'+str(index)+'\n'
        fasta_out.write(header)
        
        amplicon_frag = degen_seq+index_seq+universal_seqs['plate_fivep_sticky']+'\n'
        fasta_out.write(amplicon_frag)


for_plate_index = './bowtie_indexes/'+experiment_name+'_for_plate'
for_plate_bowtie_index_call = 'bowtie-build -f -r -o 4 '+plate_pcr_amplicon_for+' '+for_plate_index
for_plate_indexing_log = subprocess.check_output(for_plate_bowtie_index_call, shell=True)
        
        
with open(plate_pcr_amplicon_rev, 'w') as fasta_out:
    
   for index in rev_plate_indexes:

        index_seq = indexes[str(index)]
        
        header = '>'+str(index)+'\n'
        fasta_out.write(header)
        
        amplicon_frag = degen_seq + index_seq+ reverse_complement(universal_seqs['plate_threep_sticky'])+'\n'
        fasta_out.write(amplicon_frag)    
    
rev_plate_index = './bowtie_indexes/'+experiment_name+'_rev_plate'
rev_plate_bowtie_index_call = 'bowtie-build -f -r -o 4 '+plate_pcr_amplicon_rev+' '+rev_plate_index
rev_plate_indexing_log = subprocess.check_output(rev_plate_bowtie_index_call, shell=True)


#%%Align for Barcode

#Align Forward plate barcode
test_align_call_1P = 'bowtie -t -v 3 -p 6 -k 1 --trim3 230 --trim5 5 --norc --chunkmbs 256 '
# adjusted for 50pb trimming
test_align_call_1P += for_plate_index+' '
test_align_call_1P += './temp/minlen250_1P.fastq '
test_align_call_1P += './temp/plate_for_align.txt'
print (test_align_call_1P)
subprocess.check_output(test_align_call_1P, shell = True)

#Align Reverse plate barcode
test_align_call_2P = 'bowtie -t -v 3 -p 6 -k 1 --trim3 230 --trim5 5 --norc --chunkmbs 256 '
# adjusted for 50pb trimming
test_align_call_2P += rev_plate_index+' '
test_align_call_2P += './temp/minlen250_2P.fastq '
test_align_call_2P += './temp/plate_rev_align.txt'
print (test_align_call_2P)
subprocess.check_output(test_align_call_2P, shell = True)



##% Calculations with the aligned plate barcodes
plate_align_for_output = './temp/plate_for_align.txt'
plate_align_rev_output = './temp/plate_rev_align.txt'


plate_index_found_for = {}
plate_index_found_rev = {}
# empty containers that will hold    read_name: index    data pairs for the forward and reverse reads that
# had valid alignments


with open(plate_align_for_output, 'r') as for_index_positives:
    for line in for_index_positives:
        # opens P1 index alignment pack file and loops through lines

        line = line.split('\t')
        # split line to extract data
        read_name = line[0].split(' ')[0]
        # Only keep the cluster coordinates and run associated info, not the illumina annotations in the read 
        # name. See http://support.illumina.com/content/dam/illumina-support/help/BaseSpaceHelp_v2/Content/Vault/
        #           Informatics/Sequencing_Analysis/BS/swSEQ_mBS_FASTQFiles.htm
        # For info on illumina read name convention
        plate_index_found_for[read_name] = line[2]
        # add read P1 index to dict

#log_entry['for_index_pos'] = len(index_found_for.keys())
# log entry: number of reads with valid P1 index alignments

with open(plate_align_rev_output, 'r') as rev_index_positives:
    for line in rev_index_positives:
        # opens P2 index alignment pack file and loops through lines

        line = line.split('\t')
        # split line to extract data
        read_name = line[0].split(' ')[0]
        # Only keep the cluster coordinates and run associated info
        plate_index_found_rev[read_name] = line[2]
        # add read P2 index to dict

#log_entry['rev_index_pos'] = len(index_found_rev.keys())
# log entry: number of reads with valid P2 index alignments


plate_both = list(set(plate_index_found_for.keys()) & set(plate_index_found_rev.keys()))
#log_entry['both_index'] = len(both_indexes)
# Use the set data format to quickly id reads with valid alignments for both P1 and P2 indexes, and stores the
# taht number for the read pack in the log
print(len(plate_both))

#%% Graph of reads

def plot_reads_per_plate(shared_index_list, index_for_dict, index_rev_dict):
    
    plate_read_counter = Counter()
    
    for read in shared_index_list:
        
        for_index = index_for_dict[read]
        rev_index = index_rev_dict[read]
        
        plate = for_index+'-'+rev_index
        
        plate_read_counter[plate] +=1
    
    plt.figure(figsize=(14,6))
    plt.subplot(121)
    plt.bar(range(len(plate_read_counter)), list(plate_read_counter.values()), align='center', color = 'blue', alpha=0.5)
    plt.xticks(range(len(plate_read_counter)), list(plate_read_counter.keys()))
    plt.xlabel('Plate barcode combinations')
    plt.ylabel('Number of reads')
    
    plt.subplot(122)
    plt.bar(range(len(plate_read_counter)), [np.log10(x) for x in plate_read_counter.values()], align='center', color = 'blue', alpha=0.5)
    plt.xticks(range(len(plate_read_counter)), list(plate_read_counter.keys()))
    plt.xlabel('Plate barcode combinations')
    plt.ylabel('Number of reads (log10)')
    
    
plot_reads_per_plate(plate_both, plate_index_found_for, plate_index_found_rev)

#%%Make the first sample dataframe, the barcodes for each pool only
sample_dataframe = pd.read_csv(empty_datasheet, index_col=0)

sample_dataframe

#%%Align Row/Column Barcodes with Bowtie

RC_pcr_amplicon_for = './amplicon_sequences/RC_pcr_for.fa'
RC_pcr_amplicon_rev = './amplicon_sequences/RC_pcr_rev.fa'

RC_for_indexes = list(sample_dataframe.RC_for_index.unique())

RC_rev_indexes = list(sample_dataframe.RC_rev_index.unique())

print(RC_for_indexes)
print(RC_rev_indexes)


with open(RC_pcr_amplicon_for, 'w') as fasta_out:
    
    for index in RC_for_indexes:

        index_seq = indexes[str(index)]
        
        header = '>'+str(index)+'\n'
        fasta_out.write(header)
        
        amplicon_frag = universal_seqs['plate_fivep_sticky']+index_seq+universal_seqs['RC_fivep_sticky']+'\n'
        fasta_out.write(amplicon_frag)
        
        
for_RC_index = './bowtie_indexes/'+experiment_name+'_for_RC'
for_RC_bowtie_index_call = 'bowtie-build -f -r -o 4 '+RC_pcr_amplicon_for+' '+for_RC_index
for_RC_indexing_log = subprocess.check_output(for_RC_bowtie_index_call, shell=True)
       
        
with open(RC_pcr_amplicon_rev, 'w') as fasta_out:
    
    for index in RC_rev_indexes:

        index_seq = indexes[str(index)]
        
        header = '>'+str(index)+'\n'
        fasta_out.write(header)
        
        amplicon_frag = reverse_complement(universal_seqs['plate_threep_sticky']) + index_seq
        amplicon_frag += reverse_complement(universal_seqs['RC_threep_sticky'])+'\n'
        fasta_out.write(amplicon_frag)    
    
rev_RC_index = './bowtie_indexes/'+experiment_name+'_rev_RC'
rev_RC_bowtie_index_call = 'bowtie-build -f -r -o 4 '+RC_pcr_amplicon_rev+' '+rev_RC_index
rev_RC_indexing_log = subprocess.check_output(rev_RC_bowtie_index_call, shell=True)

#%% Align RC barcodes with Bowtie2
test_align_call_for = 'bowtie -t -v 3 -p 6 -k 1 --trim3 219 --trim5 28 --norc --chunkmbs 256 '


test_align_call_for += for_RC_index+' '
test_align_call_for += './temp/minlen250_1P.fastq '
test_align_call_for += './temp/RC_for_align.txt'

subprocess.check_output(test_align_call_for, shell = True)



test_align_call_rev = 'bowtie -t -v 3 -p 6 -k 1 --trim3 221 --trim5 26 --norc --chunkmbs 256 '

test_align_call_rev += rev_RC_index+' '
test_align_call_rev += './temp/minlen250_2P.fastq '
test_align_call_rev += './temp/RC_rev_align.txt'

subprocess.check_output(test_align_call_rev, shell = True)

#%%Calculations with Row/Column barcode alignements

RC_align_for_output = './temp/RC_for_align.txt'
RC_align_rev_output = './temp/RC_rev_align.txt'


index_found_for_RC = {}
index_found_rev_RC = {}
# empty containers that will hold    read_name: index    data pairs for the forward and reverse reads that
# had valid alignments


with open(RC_align_for_output, 'r') as for_index_positives:
    for line in for_index_positives:
        # opens P1 index alignment pack file and loops through lines

        line = line.split('\t')
        # split line to extract data
        read_name = line[0].split(' ')[0]
        # Only keep the cluster coordinates and run associated info, not the illumina annotations in the read 
        # name. See http://support.illumina.com/content/dam/illumina-support/help/BaseSpaceHelp_v2/Content/Vault/
        #           Informatics/Sequencing_Analysis/BS/swSEQ_mBS_FASTQFiles.htm
        # For info on illumina read name convention
        index_found_for_RC[read_name] = line[2]
        # add read P1 index to dict

#log_entry['for_index_pos'] = len(index_found_for.keys())
# log entry: number of reads with valid P1 index alignments

with open(RC_align_rev_output, 'r') as rev_index_positives:
    for line in rev_index_positives:
        # opens P2 index alignment pack file and loops through lines

        line = line.split('\t')
        # split line to extract data
        read_name = line[0].split(' ')[0]
        # Only keep the cluster coordinates and run associated info
        index_found_rev_RC[read_name] = line[2]
        # add read P2 index to dict

#log_entry['rev_index_pos'] = len(index_found_rev.keys())
# log entry: number of reads with valid P2 index alignments


both_indexes_RC = list(set(index_found_for_RC.keys()) & set(index_found_rev_RC.keys()))
#log_entry['both_index'] = len(both_indexes)
# Use the set data format to quickly id reads with valid alignments for both P1 and P2 indexes, and stores the
# that number for the read pack in the log
print(len(both_indexes_RC))


#%% Assemble the indexes of plate PCR and Row/Column PCR into a dataframe

valid_plates = [tuple((21,1))]

def plot_reads_in_plate(plate_both, RC_both, plate_index_found_for, plate_index_found_rev, index_found_for_RC, index_found_rev_RC, plate,silent=True):
    
    plate_rows = sorted(set(index_found_for_RC.values()))
    plate_columns = sorted(set(index_found_rev_RC.values()))
    
    read_plate_RC_dict = {}
    grid_dict = {}
    
    for row in plate_rows:
        
        grid_dict[int(row)] = {}
        
        #print(row)
        
        for column in plate_columns:
            
            #print(type(column))
            
            grid_dict[int(row)][int(column)] = 0
            
    if silent == False:
        
        print(grid_dict)
    
    plate_and_rc_both = list((set(plate_both)&set(both_indexes_RC)))
    
    for read in plate_and_rc_both:
        
        for_index_plate = int(plate_index_found_for[read])
        rev_index_plate = int(plate_index_found_rev[read])
        
        #print(for_index_plate, rev_index_plate)
        
        #print(plate)
        #print(tuple((for_index_plate,rev_index_plate)))
        
        if tuple((int(for_index_plate), int(rev_index_plate))) == plate:
            
            #print('ok')
            
            for_index_RC = int(index_found_for_RC[read])
            rev_index_RC = int(index_found_rev_RC[read])
            
            #print(for_index_RC, rev_index_RC)
            #print(grid_dict[for_index_RC][rev_index_RC])
            
            grid_dict[for_index_RC][rev_index_RC] += 1
            
            read_plate_RC_dict[read] = [int(for_index_plate), int(rev_index_plate), for_index_RC, rev_index_RC]
            
            
            
    plate_df = pd.DataFrame.from_dict(grid_dict, orient='index')
    return plate_df, read_plate_RC_dict
            
###I
RC_index = plot_reads_in_plate(plate_both, both_indexes_RC, plate_index_found_for, plate_index_found_rev, index_found_for_RC, index_found_rev_RC, valid_plates[0], silent = False)



#%% Make a figure of reads per Row/Column barcode

col_order = [1,2,3]

plt.figure(figsize=(14/1.5, 10/1.5))


sns.heatmap(RC_index[0][col_order], cmap='viridis')
plt.xlabel('Column Barcode', fontsize = 16)
plt.ylabel('Row Barcode', fontsize = 16)
plt.xticks(fontsize = 14)
plt.yticks(fontsize = 14)

plt.axhline(2, color='white')
plt.axhline(4, color='white')
plt.axhline(6, color='white')
plt.axhline(8, color='white')
plt.title('Reads per pool')

#
#%% Preparation for demultiplexing reads
all_plate_index_dicts = []

for plate in valid_plates:
    
    mapping = plot_reads_in_plate(plate_both, both_indexes_RC, plate_index_found_for, plate_index_found_rev, index_found_for_RC, 
                                index_found_rev_RC, plate, silent=True)[1]
    
    all_plate_index_dicts.append(mapping)
    

del(valid_plates)
#%% Make a folder for the demultiplexed reads

os.mkdir("./demultiplexed_reads")

#%% Create files to sort the demultiplexed reads into

plate_RC_combi_to_dict = {}

out_fastq_path_dict_for_df = {}

out_fastq_path_dict_demult = {}

dir_path = './demultiplexed_reads/'


for sample in sample_dataframe.index:
    
    plate_for = sample_dataframe.loc[sample]['PE1_index']
    RC_for = sample_dataframe.loc[sample]['RC_for_index']
    RC_rev = sample_dataframe.loc[sample]['RC_rev_index']
    plate_rev = sample_dataframe.loc[sample]['PE2_index']
    
    plate_RC_combi_to_dict[sample] = [plate_for, plate_rev, RC_for, RC_rev]
    
    out_path = dir_path+str(sample)+'_indexes_'+str(plate_for)+'_'+str(plate_rev)+'_'+str(RC_for)+'_'+str(RC_rev)
    
    out_path_1 = dir_path+str(sample)+'_indexes_'+str(plate_for)+'_'+str(plate_rev)+'_'+str(RC_for)+'_'+str(RC_rev)+'_for'+'.fastq'
    out_path_2 = dir_path+str(sample)+'_indexes_'+str(plate_for)+'_'+str(plate_rev)+'_'+str(RC_for)+'_'+str(RC_rev)+'_rev'+'.fastq'
    
    out_fastq_path_dict_for_df[sample] = out_path
    
    out_fastq_path_dict_demult[tuple([plate_for, plate_rev, RC_for, RC_rev])] = out_path
    
    create_file_for = 'touch '+ out_path_1
    create_file_rev = 'touch '+ out_path_2
    
    subprocess.check_output(create_file_for, shell=True)
    subprocess.check_output(create_file_rev, shell=True)
  


for plate_dict in all_plate_index_dicts:
    
    print(len(plate_dict.keys()))

    for read in list(plate_dict.keys()):

        index_combi = plate_dict[read]

        if index_combi not in plate_RC_combi_to_dict.values():

            plate_dict.pop(read, None)
        
    print(len(plate_dict.keys()))

#%% Demultiplexing

def split_reads_from_fastq(file, file_type, plate_mapping_list):
    
    to_write = 'none'
    read_name = 'none'
    
    
    with open(file, 'r') as source_fastq:
        
        for line in source_fastq:
            
            #print(read_name)
            #print(line)

            if line.startswith('@M') == True:
                
                if to_write == 'none':
                    
                    read_name = line.split(' ')[0].strip('@')
                    to_write = line

                else:
                    
                    for plate_read_dict in plate_mapping_list:
                        
                                                
                        if read_name in plate_read_dict.keys():

                            indexes = tuple(plate_read_dict[read_name])
                            
                            if file_type == 'forward':
                                
                                file_path = out_fastq_path_dict_demult[indexes]+'_for.fastq'
                                
                            elif file_type == 'reverse':
                                
                                file_path = out_fastq_path_dict_demult[indexes]+'_rev.fastq'

                            
                            
                            with open(file_path, 'a') as dest:
                                
                                   dest.write(to_write)

                            break





                    read_name = line.split(' ')[0].strip('@')

                    to_write = line
                

            else:

                to_write += line
                
split_reads_from_fastq('./temp/minlen250_1P.fastq', 'forward', all_plate_index_dicts)
split_reads_from_fastq('./temp/minlen250_2P.fastq', 'reverse', all_plate_index_dicts)

#%% Associate each entry in dataframe with a demultiplexed file

for sample_ID in sample_dataframe.index:
    print(sample_ID)
        
    sample_dataframe.loc[sample_ID, 'filepath'] = f"./demultiplexed/sample_{sample_ID}"


del(sample_ID)
#%% Merge the reads

if './merged_reads' not in os.listdir('./'):
    os.mkdir('./merged_reads')
    
    
before_merge_depth_dict = {}


def merge_reads(Sample_ID):
    
    filepath = sample_dataframe.loc[Sample_ID]['filepath']
    
    filepath_for = filepath+'_for.fastq'
    filepath_rev = filepath+'_rev.fastq'
    
    sample_read_count = 0
    
    with open(filepath, 'r') as source:
        
        for line in source:
            
            if line.startswith('>'):
                
                sample_read_count += 1
                
    before_merge_depth_dict[Sample_ID] = sample_read_count
                
    filepath_out = './merged_reads/'+str(Sample_ID)+'.fasta'

    panda_seq_call = 'pandaseq -f '+filepath_for+' -r '+filepath_rev+ ' -L 550 -O 400 -k 4 -B -N -t 0.5 -T 6 -g log.txt -w '+ filepath_out

    subprocess.check_output(panda_seq_call, shell=True)
    
        
    return sample_read_count
    
#%% Run the merge_reads function
for sample_ID in list(sample_dataframe.index):
    
    merge_reads(sample_ID)
    
#%% Add read count info to the dataframe
sample_dataframe['demult_read_count'] = pd.Series(before_merge_depth_dict)


#%% Count merged reads per library

after_merge_depth_dict = {}

def count_merged(Sample_ID):
    
    filepath = './demultiplexed_reads/'
    
    filepath_merged = 'sample_' +str(Sample_ID)+'_for.fastq'
    
    filepath+=filepath_merged
    
    after_merge_depth = 0
    
    
    with open(filepath, 'r') as source:

        for line in source:

            if line.startswith('@'):

                after_merge_depth += 1

    after_merge_depth_dict[Sample_ID] = after_merge_depth
    
    return after_merge_depth
    
    
for sample_ID in list(sample_dataframe.index):
    
    count_merged(sample_ID)
    sample_dataframe.loc[sample_ID, "merged_file"]=f"./merged/sample_{sample_ID}.fasta"
    
sample_dataframe['merged_read_count'] = pd.Series(after_merge_depth_dict)
# add read count info to the dataframe



#%% Save this version of the dataframe to csv



sample_dataframe.to_csv(demultiplexed_datasheet, sep=',')