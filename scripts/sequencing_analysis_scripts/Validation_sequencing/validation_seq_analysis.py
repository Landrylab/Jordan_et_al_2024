#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 15:37:39 2024

@author: adminlandry
"""

#Using the demultiplexed reads obtained using the validation_seq_demultiplexing.py script and uploaded to the SRA,
#this script counts how many times each Pbs2 mutation appears in each library

#%% Importing packages
import os

import subprocess

import pandas as pd
print(pd.__name__, pd.__version__)

import numpy as np
print(np.__name__, np.__version__)

import scipy
print(scipy.__name__, scipy.__version__)

#from collections import Counter
#import regex as regex

import seaborn as sns
print(sns.__name__, sns.__version__)

#%% Set folder

###Set this to the folder containing the script. It should also contain the sample_datasheet_demultiplexed_merged.csv file, which is obtained from the demultiplexing script. it should also contain the demultiplexed reads from the SRA, in a folder named merged_reads. The common_sequences folder should be in the folder above this one.
 
os.chdir("/your/path")

print(os.listdir())


#%% Import library information datasheet

#obtained
sample_dataframe=pd.read_csv('./validation_sample_datasheet_demultiplexed_merged.csv', index_col=0)


#%% Trim and aggregate reads

if "trimmed_reads" not in os.listdir("./"):
    os.mkdir("./trimmed_reads")
    
if "aggregated_reads" not in os.listdir("./"):
    os.mkdir("./aggregated_reads")

def trim_aggregate(Sample_ID):
    # This function generates and runs two vsearch calls to trim and then aggregate the merged reads.
    # Trimming removes the variable degenerate sequences, which then  allows us to aggregate the sequences that have the same sequence with minimal influence from amplicon regions used for multiplexing.
       
    filepath_in = './merged_reads/'    
    filepath_merged = filepath_in+str(Sample_ID)+'.fasta'
    # get the merged read filepath based on the sample
    trim_folder='./trimmed_reads/'
    agg_folder='./aggregated_reads/'
    
    trim_path = trim_folder + str(Sample_ID)+ '_trim.fasta'
    aggregate_path = agg_folder +str(Sample_ID)+ '_agg.fasta'
    # generate output filepaths for the trimming and aggregating steps
    
    vsearch_trim_call = 'vsearch --fastx_filter '+filepath_merged+' --fastq_stripleft 18 --fastq_stripright 18 --fastaout ' + trim_path
    subprocess.check_output(vsearch_trim_call, shell=True)
    # generates and runs a vsearch call with the following parameters:
    #    --fastx_filter filepath_merged    vsearch program to be used and input: fastxfilter can perform different
    #                                      operations on fastq files
    #    --fastq_stripleft 18              removes 18 bases from the 5p end
    #    --fastq_stripright 18             removes 18 bases from the 3p end
    #    --fastaout trim_path              ouput trimmed sequences the trim_path file                 
    
    vsearch_aggregate_call = 'vsearch --fastx_uniques '+ trim_path +' --relabel seq --fastaout '+aggregate_path+' --sizeout'
    subprocess.check_output(vsearch_aggregate_call, shell=True)
    # generates and runs a vsearch call with the following parameters:
    #    --fastx_uniques trim_path    vsearch program to be used and input. fastx_uniques aggregates perfectly identical 
    #                                    sequences in the input fasta, dramatically reducing the number of sequences that then 
    #                                    need to be aligned
    #    --relabel seq                   changes header so that all sequence names begin by seq
    #    --output                        path where output will be written
    #    --sizeout                       append the size of the sequence cluster to the fasta header
    
    
for sample_ID in list(sample_dataframe.index):
    
    trim_aggregate(sample_ID)
    
del(sample_ID)
    
#%% Align to Pbs2 WT sequence

if 'amplicon_align' not in os.listdir("./"):
    os.mkdir('./amplicon_align')
    
    
for Sample_ID in list(sample_dataframe.index):
    
    print(Sample_ID)
    #Define file paths
    file_in="./aggregated_reads/"+ str(Sample_ID) + "_agg.fasta"
    file_ref="../common_sequences/Pbs2_sequence.fasta"
    file_out= f'./amplicon_align/lib{Sample_ID}_aligned.needle'
    
    
    #Write command line to run the needle alignment tool
    cmdline_needle=f'needle -auto -gapopen 50 -asequence {file_ref} -bsequence {file_in} -aformat markx10 -outfile {file_out}'
    #print(cmdline_needle)
    #Run command 
    subprocess.check_output(cmdline_needle, shell = True)
    
del(file_in, file_out, Sample_ID, file_ref, cmdline_needle)


#%% define function to parse the alignement file
def parse_needle_output(Sample_ID):
    # this function parses the needle alignments and extract the aligned sequences reference and query seuqences. It takes as input the 
    # path to the needle output
    
    n_aligns = 0
    # counter for the number of alignments
    
    align_seqs_dict = {}
    # empty container that will hold the aligned sequences
                
    needle_align_path = "./amplicon_align/lib"+ str(Sample_ID) + "_aligned.needle"
    # path to the alignments
        
    with open(needle_align_path, 'r') as source:
        # open the alignment
        
        current_align = ''
        current_qseq = ''
        current_sseq = ''
        # empty container objects for data processing
        
        qseq_done = 0
        #counter for the number of alignments processed
        
        for line in source:
            # loop through the file
            
            if line.startswith('>>>') == True:
                # detect headers
                
                n_aligns +=1
                # increment alignment counter by one
                
                align_name = line.strip('>>>')
                   
                if n_aligns != 1:
                    # if this is not the first alignment
                    
                    align_seqs_dict[current_align] = [current_qseq, current_sseq]
                    # add the information on the previous alignment to the dict
                    
                    current_align = align_name
                    # update the name for the new entry
                    
                    current_qseq = ''
                    current_sseq = ''
                    # reset temporary variables
                    
                    qseq_done = 0
                    # reset indicator for query sequence extraction
                    
                else:
                    current_align = align_name
                    # for the fisrt sequence, just need to store the align name
                    
                    
            elif line.startswith(';') == False and line.startswith('>') == False and line.startswith('\n') == False and line.startswith('#') == False:
                # skip all the useless lines to process only the aligned sequences
                
                if qseq_done == 1:
                    current_sseq += line.strip('\n')
                    # if the query seq is done (qseq = 1), add sequence to the subject
                    
                else:
                    current_qseq += line.strip('\n')
                    # if the query seq is not done, continue to update it
                    
            elif line.startswith('#--') == True:
                align_seqs_dict[align_name] = [current_qseq, current_sseq]
                # update dict with info from the last entry in the alignment sequence
                
            else:
                if qseq_done == 0 and current_qseq != '':
                    qseq_done =1
                    # if the qseq is recorded, update value
                
                    
    return align_seqs_dict, n_aligns

#%% Define function to identify mutations from alignement file

def find_mutations(Sample_ID):
    # For a given path and reference sequence, runs the Needle alignment parser and
    # then goes through the extracted alignments to find the differences between the aligned query and the reference.
    #
    
    dict_of_mut_dicts[Sample_ID]={}
    #allele_dict = {}
    # empty dict that will hold the sequence name :mutations information
   
    
    align_dict, align_count = parse_needle_output(Sample_ID)
    # get the parsed alignments for the sample
    
    for entry in list(align_dict.keys()):
        # loop through aligned sequences
        
        read_var_list = []
        # temporary holder for mutations found in sequence
        
        query_seq = align_dict[entry][1]
        # aligned cds sequence of the strain
        
        align_ref = align_dict[entry][0]
        # aligned cds sequence of the reference
        
        gap_adjust = 0
        # value used to adjust the cds sequence index for the presence of insertions in the strain sequence vs the 
        # reference cds
        
        backtrack_adjust = 0
        # value used to adjust the cds sequence index for the presence of deletions in the strain sequence vs the 
        # reference cds
        
        temp_var = None
        # temporary variable to hold the sequence of an insertion or deletion as a string. When the gap ends, annotation 
        # will be added to read_var_list
        
        indel_start = 0
        # position start of the indel annotation in the reference sequence, with adjustment for gap presence
        
        ref_seq_no_gaps = align_ref.replace('-','')
        # reference sequence with gaps removed
        
        align_start = (amplicon_dict[ref_orf].index(ref_seq_no_gaps))+1
        # position in the reference sequence where the alignment starts
        
        query_seq_no_gaps = len(query_seq.replace('-',''))
        # length of the query when gaps are removed
        
        for nt in range(0, len(align_ref)):
            # iterates through the entire alignment of the strain prot sequence
            
            if query_seq[nt] == '-':
                # detect a deletion variant
                
                # logic for indel detection/annotation:
                #
                # suppose we have this alignment  
                #
                # 1 2 3 4 5 6 7 8 9
                # A T - - A A A T G    strain variant: del gaps are indexed because the aa index is based on reference
                # A T K P A - - T G
                # 1 2 3 4 5     6 7    reference: insert gaps not indexed because aa positions do exist in reference
                #
                # following this logic, every time an insertion is detected and annotated, the gap_adjust value is 
                # incremented by the length of the gap and used to adjust the variant mapping to make it match the 
                # reference index values. The indel aa postion is the first residue detected as part of the indel
                
                
                if indel_start == 0:
                    # checks if the character is the start or the continuation of a gap in the alignment
                    
                    temp_var = 'del'+ align_ref[nt]
                    indel_start = (nt+1-gap_adjust)
                    # if it is, starts a new annotation entry with a start position compensated for previous insertions
                    # (if any)
                    
                    backtrack_adjust += 1
                    
                else:
                  
                    temp_var += align_ref[nt]
                    # if it is not, adds the following aa to the deletion annotation
                    
                    backtrack_adjust += 1
                    
                    
            elif align_ref[nt] == '-':
                # detects an insertion variant
                
                if indel_start == 0:
                    # checks if the character is the start or the continuation of a gap in the alignment
                    
                    temp_var = 'ins'+ query_seq[nt]
                    
                    indel_start = (nt+1-gap_adjust)
                    # if it is, starts a new annotation entry with a start position compensated for previous insertions
                    # (if any)
                    
                    gap_adjust += 1
                    # increments the gap adjust for the this added aa in the strain sequence                   
                    
                    
                else:
                  
                    temp_var += query_seq[nt]
                    # if it is not, adds the following aa to the insertion annotation
                    
                    gap_adjust += 1
                    # increments the gap adjust for the this added aa in the strain sequence
                    
                    
            elif query_seq[nt] != align_ref[nt]:
                # detects a mismatch between the strain sequence and the reference
                
                
                variant = align_ref[nt]+'|'+str((nt+1-gap_adjust))+'|'+query_seq[nt]
                read_var_list.append(variant)
                # creates an annotation for the strain-reference aa mismatch and appends it to the list of 
                # annotations
                
            else:
              
                 if indel_start != 0:
                    # detects if there is currently an open gap entry. If there is, then the detected mismatch means 
                    # that it has now concluded
                    
                    read_var_list.append(str((indel_start))+temp_var)
                    temp_var = None
                    indel_start = 0
                    # adds the indel annotation to the strain variant list and resets temporary variables for the next 
                    # indel entry
                    
                    
        if query_seq_no_gaps >=  len(ref_seq_no_gaps)*0.8 and len(read_var_list)<25:            
            dict_of_mut_dicts[Sample_ID][entry] = read_var_list, align_start
            # apply a filter for alignment quality: the alignment must cover at least 80% of the reference sequence and
            # there must be less than 25 differences between the query and the reference (insertions, deletions or SNVs)
                           
    #return dict_of_mut_dicts[Sample_ID]

#%% Identify mutations in all aligned sequences

dict_of_mut_dicts={}
amplicon_dict={"Pbs2": "AGTAATCAAAGCGAGCAAGACAAAGGCAGTTCACAATCACCTAAACATATTCAGCAGATTGTTAATAAGCCATTGCCGCCTCTTCCCGTAGCAGGAAGTTCTAAGGTTTCACAAAGAATGAGTAGCCAAGTCGTGCAAGCGTCCTCCAAGAGCACTCTTAAGAACGTT"}
ref_orf="Pbs2"

for Sample_ID in list(sample_dataframe.index):
    print(Sample_ID)
    find_mutations(Sample_ID)

del(Sample_ID)
#%% Get variant count from mutation information

def get_variant_count(mutation_set, ref_seq, codon_start, n_aa):
    
    
    for Sample_ID in list(sample_dataframe.index):
        var_dict_of_dicts[Sample_ID] = {}
        
        variants = list(dict_of_mut_dicts[Sample_ID].keys())
        codon_groups = {}
        
        codon = 0
       
        
        wt_count =0
        valid_seq=0
        
        for nt in range(0, n_aa*3):
             pos = nt+1
            #so the indexing does not start at 0, and mess up the numbering with regards to the codons
             if nt % 3 == 0:
                 codon += 1
                
             codon_groups[pos] = codon
            
             var_dict_of_dicts[Sample_ID][codon] = {}
            
        wt_codons = {}
            
        ref = amplicon_dict[ref_seq]
        
        for aa in range(0, n_aa):
            offset = 0#33
            
            start = offset+(aa*3)
            
            wt_codon=ref[start:(start+3)]
            
            wt_codons[(aa+codon_start)] = wt_codon
            
            #This is now a dict of dict of a dict. For each codon, there is a dict of all variants 
            var_dict_of_dicts[Sample_ID][aa+codon_start][wt_codon]=np.nan
            
            
            
        for variant in variants:
            
            var_info = variant.split(',')
            var_count =int(var_info[1].split(';')[1].strip('size='))
                   
            mut_list = mutation_set[Sample_ID][variant][0]
            #mutation_set is pool_pos_muts
            
            
            if var_count>=20:
                #print(mut_list)
                
                for mutation in mut_list:
                  
                    if 'del' in mutation:
                        mut_info = mutation.split('del')
                        mut_pos = int(mut_info[0])
                        
                        if mut_pos == 1:
                            mut_list.remove(mutation)
                            
                if len(mut_list) <=3 and 'ins' not in str(mut_list) and 'del' not in str(mut_list):
                    
                    if len(mut_list) ==0:
                        wt_count += var_count
                        
                    else:
                        #print(mut_list)
                        mut_nt_list = []
                        
                        out_list = []
                        
                        for mutation in mut_list:
                            
                            mut_pos = int(mutation.split('|')[1])
                            
                            if mut_pos >= codon_start and mut_pos <=(n_aa*3):
                            
                                mut_nt_list.append(codon_groups[mut_pos])
                                
                            else:
                                out_list.append(mut_pos)
                            
                        if len(set(mut_nt_list)) == 1:
                            valid_seq+=var_count
                            
                            codon = int(list(set(mut_nt_list))[0])
                            # info
                            wt_seq = wt_codons[codon]
                                                                          
                            new_seq = [x for x in wt_seq]
                            
                            for mutation in mut_list:
                                mut_pos = int(mutation.split('|')[1])
                                mutation = mutation.split('|')[2]
                                
                                codon_pos = (mut_pos-1)%3
                                
                                new_seq[codon_pos] = mutation
                                
                            new_codon = ''.join(new_seq)
                                
                            #print(wt_seq, mut_list, new_codon)
                            
                            if new_codon in list(var_dict_of_dicts[Sample_ID][codon].keys()):
                                
                                
                                var_dict_of_dicts[Sample_ID][codon][new_codon]+=var_count
                                
                            else:
                                var_dict_of_dicts[Sample_ID][codon][new_codon]=var_count
                                
                                
                            
                        elif len(set(mut_nt_list)) == 0 and len(out_list)>=1:
                            #print(mut_nt_list, out_list)
                            wt_count+=var_count
                        
                        
            elif var_count < 20:
              
                for mutation in mut_list:
                  
                    if 'del' in mutation:
                        mut_info = mutation.split('del')
                        mut_pos = int(mut_info[0])
                        
                        if mut_pos == 1:
                            mut_list.remove(mutation)
                        
                if len(mut_list) <=3 and 'ins' not in str(mut_list) and 'del' not in str(mut_list):
                  
                    if len(mut_list) ==0:
                        wt_count += var_count
                        
                    else:
                        #print(mut_list)
                        mut_nt_list = []
                        
                        out_list = []
                        
                        for mutation in mut_list:
                          
                            mut_pos = int(mutation.split('|')[1])
                            
                            if mut_pos >= (codon_start*3) and mut_pos <= ((codon_start*3)+(n_aa*3)):
                            #if mut_pos >=1 and mut_pos <= 48:
                            
                                mut_nt_list.append(codon_groups[mut_pos])
                                
                            else:
                                out_list.append(mut_pos)
                                
                        if len(set(mut_nt_list)) == 1:
                            #print(mut_list, var_count)
                            valid_seq+=var_count
                            
                            codon = int(list(set(mut_nt_list))[0])
                            
                            wt_seq = wt_codons[codon]
                            
                            new_seq = [x for x in wt_seq]
                            
                            for mutation in mut_list:
                                mut_pos = int(mutation.split('|')[1])
                                mutation = mutation.split('|')[2]
                                
                                codon_pos = (mut_pos-1)%3
                                
                                new_seq[codon_pos] = mutation
                                
                            new_codon = ''.join(new_seq)
                            
                            #print(wt_seq, mut_list, new_codon)
                            
                            if new_codon in list(var_dict_of_dicts[Sample_ID][codon].keys()):
                                
                                
                                var_dict_of_dicts[Sample_ID][codon][new_codon]+=var_count
                                
                                
                            else:
                                var_dict_of_dicts[Sample_ID][codon][new_codon]=var_count
                                
                            
              
              
        wt_count_dict[Sample_ID] = wt_count    
        print(Sample_ID)
    return var_dict_of_dicts[Sample_ID], wt_count

#%% Get counts per variant
var_dict_of_dicts = {}
wt_count_dict = {}


#mutation_set = dict of dict per library
#ref_seq= Sequence that is mutated (key in amplicon dictionary)
#codon_start
#n_aa = Number of amino acids in sequence which is coded for 
#def get_variant_count_1(mutation_set, ref_seq, codon_start, n_aa):
get_variant_count(dict_of_mut_dicts, 'Pbs2', 1, 56)

del(amplicon_dict, ref_orf)
#%% Convert variant count dict to df

if not os.path.exists("./variant_counts/"):
    os.makedirs("./variant_counts/")

#Transform dicts of dicts into dict of dataframes, then save dataframes as csv files
df_dict = {}
for Sample_ID in list(sample_dataframe.index):
  
    variant_file = f'./variant_counts/library{Sample_ID}.csv'
    sample_dataframe.loc[Sample_ID, "variant_file"]=variant_file
    df_dict[Sample_ID] = pd.DataFrame(var_dict_of_dicts[Sample_ID]) 
    df_dict[Sample_ID].to_csv(variant_file, sep = ',')
 
#Transform WT count dictionary to dataframe, then save as csv file     
wt_count_df = pd.DataFrame.from_dict(wt_count_dict, columns=["wt_count"], orient = 'index')
wt_count_df.to_csv('./variant_counts/wt_counts.csv', sep=',')


#%% Copy and save sample dataframe

for i in list(sample_dataframe.index):
    sample_dataframe.loc[i, "wt_reads"] = wt_count_df.loc[i, "wt_count"]
    sample_dataframe.loc[i, "variant_reads"]=df_dict[i].fillna(0).sum().sum()
del(i)    
sample_dataframe["final_reads"]= sample_dataframe["wt_reads"] + sample_dataframe["variant_reads"]

sample_dataframe.to_csv("./sample_sheet_final.csv", index=True)