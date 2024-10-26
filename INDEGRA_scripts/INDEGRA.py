#!/usr/bin/env python3


import argparse
import gffutils
import logging
import os
import shutil
import pandas as pd
import numpy as np 
import pysam
import tempfile
import time
from tqdm import tqdm
import math
from math import *
import scipy.stats as stats
from pathlib import Path


TestSize = 30
Censor=150
NBins = 10

a=1
b=8
s=1800
def DTI_Function(x):    
    if x==0:
        return 10    
    if 1/x<s:
        return a+(b-a)*(1/x-1)**(1/2)/(s-1)**(1/2)
    else:
        return 10+(b-10)*(1/x-1)**(-(b-1)/2/(10-b))/(s-1)**(-(b-1)/2/(10-b))




# set logging parameters 
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S')

# function to parse bam :: for each primary alignment, returns read_id, sequence_length, and then length of all cigar operations 
def parse_bam(bam_file):
    start = time.time()

    # use pysam to read the bam
    with pysam.AlignmentFile(bam_file, "rb") as bam_file:
        
        # create a dict for bam data
        bam_data = {}

        # loop over records in the bam file
        for read in bam_file.fetch():
            # discard secondary/supplementary
            if read.is_secondary or read.is_supplementary:
                continue

            # basecalled sequence length
            seq_len = len(read.query_sequence)

            # get cigar operations
            cigar_ops = read.cigartuples

            # calculate cumulative length for each cigar operator - note, 0 = match, 1 = ins, 2 = del, 3 = splice, 4 = softclip, 5 = hardclip 
            cigar_lengths = {}
            cigar_maxs = {0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0}
            for i, cigartuple in enumerate(cigar_ops):
                cigar_op, cigar_len = cigartuple
                if cigar_op in cigar_lengths:
                    cigar_lengths[cigar_op] += cigar_len
                else:
                    cigar_lengths[cigar_op] = cigar_len
                if cigar_op in cigar_maxs:
                    if cigar_len > cigar_maxs[cigar_op]:
                        cigar_maxs[cigar_op] = cigar_len
                else:
                    cigar_maxs[cigar_op] = cigar_len

            # filter if max insertion >80
            if cigar_maxs[1]>args.ins:
                continue
            # filter if max deletion >160
            if cigar_maxs[2]>args.del:
                continue
            # filter if soft clip >200
            if cigar_maxs[4]>args.sc:
                continue

            # strip version number from transcript id
            #reference_transcript = read.reference_name.split('.')[0]
            reference_transcript = read.reference_name

            # write this data into our bam dict
            bam_data[read.query_name] = {"sequence_length": seq_len, "cigar_lengths": cigar_lengths, "reference_transcript": reference_transcript, "start_position": read.reference_start, "end_position": read.reference_end}

        elapsed_time = round(time.time() - start, 2)
        logging.info(f'Parsed {len(bam_data)} lines from BAM file in {elapsed_time} seconds')

        # print a preview of bam_data if verbosity == 2
        if args.verbosity == 2:
            header = "\n".join([f"{key}: {value}" for key, value in list(bam_data.items())[:5]])
            logging.info(f"BAM Data Header:\n{header}")

    return bam_data 



# function to merge data from bam, gtf, sequencing_summary, using transcript_id as key 
# includes the MDI_read_length, the sum of cigar match, deletion, insertion 
def merge_data(gtf_data, bam_data, summary_data, sample_name):
    start = time.time()

    # header for merged data 
    merged_data = "sample_name\tread_id\tref_transcript\tsequence_length\tcigar_match\tcigar_ins\tcigar_del\tcigar_splice\tcigar_softclip_total\tcigar_hardclip_total\tstart_position\tend_position\tend_reason\ttranscript_annotated_length\ttranscript_biotype\n"

    # make list to store any transcripts which have merging error with gtf data 
    gtf_error_transcripts = []
    
    # make list for transcripts which have a merging error with sequencing_summary 
    ss_error_transcripts = []

    # create list to store read names which have discrepancy between cigar lengths and read length
    cigar_discrepancy_reads = []
    cigar_discrepancy_dict = {}

    # loop over read_id in bam_data 
    for read_id, bam_info in bam_data.items():

        # select reference transcript 
        ref_transcript = bam_info["reference_transcript"]

        # as a sanity check, get total cigar length and compare to read length. They should be the same. 
        read_length = bam_info.get("sequence_length", 0) 
        
        # get cigar lengths - note, 0 = match, 1 = ins, 2 = del, 3 = splice, 4 = softclip, 5 = hardclip 
        cigar_lengths = bam_info.get("cigar_lengths", {})
        cigar_values = "\t".join(str(cigar_lengths.get(i, 0)) for i in range(6))
        
        # get start and end positions
        start_position = bam_info.get("start_position", {})
        end_position = bam_info.get("end_position", {})
        
        # check if aligned length matches read length 
        #aligned_length = sum(cigar_lengths[0] + cigar_lengths[1] + cigar_lengths[3] + cigar_lengths[4] + cigar_lengths[5] + 1)
        aligned_length = cigar_lengths.get(0, 0) + cigar_lengths.get(1, 0) + cigar_lengths.get(3, 0) + cigar_lengths.get(4, 0) + cigar_lengths.get(5, 0) 

        # compare the CIGAR length and sequenced read length
        if aligned_length != read_length:
            cigar_discrepancy_reads.append(read_id) # if there's an error, write the read_id to our discrepency list 
            cigar_discrepancy_dict[read_id] = aligned_length - read_length # note the size of the error as well 

        # now, check if the reference transcript in the bam is represented in gtf_data, else print to gtf_error_transcript 
        if ref_transcript in gtf_data:
            gtf_info = gtf_data[ref_transcript] # fetch the annotated transcript length and the annotated biotype 
            
            # now, check if the read id is in the sequencing_summary file 
            if read_id in summary_data:
                end_reason = summary_data[read_id] # fetch the read_end_reason 
            
            else:
                ss_error_transcripts.append(ref_transcript)
        else:
            gtf_error_transcripts.append(ref_transcript)

        # finally, combine everything into merged_data
        merged_data += f'{sample_name}\t{read_id}\t{ref_transcript}\t{read_length}\t{cigar_values}\t{start_position}\t{end_position}\t{end_reason}\t{gtf_info["length"]}\t{gtf_info["biotype"]}\n'


    # log read/transcripts pairs that have errors merging with the gtf 
    # note, we get a few in the test data, because the gtf is chromosome only whereas the FASTA that was used has scaffolds, alt contigs etc. 
    if gtf_error_transcripts:
        unique_gtf_error_transcripts = set(gtf_error_transcripts)
        with open(os.path.join(os.path.dirname(args.output_file), 'gtf_merge_errors.txt'), 'w') as f:
            f.write("\n".join(unique_gtf_error_transcripts))
        logging.error(f'Error merging data for {len(gtf_error_transcripts)} reads with gtf data. Error transcripts written to gtf_merge_errors.txt in the output directory.')
    
    # similarly, log read/transcript pairs with errors merging against the sequencing summary 
    if ss_error_transcripts:
        unique_ss_error_transcripts = set(ss_error_transcripts)
        with open(os.path.join(os.path.dirname(args.output_file), 'ss_merge_errors.txt'), 'w') as f:
            f.write("\n".join(unique_ss_error_transcripts))
        logging.error(f'Error merging data for {len(ss_error_transcripts)} reads with sequencing summary. Error transcripts written to ss_merge_errors.txt in the output directory.')

    # record cases where we have any discrepency in the cigar
    if cigar_discrepancy_reads:
        with open(os.path.join(os.path.dirname(args.output_file), '_cigar_discrepancy_errors.txt'), 'w') as f:
            for read_id, discrepancy in cigar_discrepancy_dict.items():
                f.write(f'{read_id}\t{discrepancy}\n')
        logging.error(f'Error in cigar lengths for {len(cigar_discrepancy_reads)} reads. Discrepancies written to cigar_discrepancy_errors.txt in the output directory.')
    
    # relay stats about the merging -- currently output by default, but could be moved to verbosity == 1
    logging.info(f'Total Reads Merged: {len(bam_data)}')
    logging.info(f'Total Transcripts Merged: {len(gtf_data)}')
    logging.info(f'Reads in BAM which clash with with GTF: {len(gtf_error_transcripts)}, transcripts these reads map to: {len(set(gtf_error_transcripts))}')
    logging.info(f'Reads in BAM which do not match sequencing summary: {len(ss_error_transcripts)}, transcripts these reads map to: {len(set(ss_error_transcripts))}')

    # print 10 lines of merged data if verbosity == 2
    if args.verbosity == 2:
        preview_data = merged_data.splitlines()[:10]
        preview_data_str = '\n'.join(preview_data)
        logging.info("First 10 lines of the merged data:\n%s", preview_data_str)
    
    # write merged data to tmpfile
    temp_file_path = None
    temp_file_path = tmp_output_file + sample_name + "_merged_data.tsv"
    with open(temp_file_path, 'w') as file:
        file.write(merged_data)
    logging.info(f'Merged data written to temporary file: {temp_file_path}')

    # update user 
    elapsed_time = round(time.time() - start, 2)
    logging.info(f'Merged data in {elapsed_time} seconds')

    return merged_data, temp_file_path


# function to merge data from bam, gtf, sequencing_summary, using transcript_id as key 
# includes the MDI_read_length, the sum of cigar match, deletion, insertion 
def merge_data_justbam(bam_data,sample_name):
    start = time.time()

    # header for merged data 
    merged_data = "sample_name\tread_id\tref_transcript\tsequence_length\tcigar_match\tcigar_ins\tcigar_del\tcigar_splice\tcigar_softclip_total\tcigar_hardclip_total\tstart_position\tend_position\n"
  
    # create list to store read names which have discrepancy between cigar lengths and read length
    cigar_discrepancy_reads = []
    cigar_discrepancy_dict = {}

    # loop over read_id in bam_data 
    for read_id, bam_info in bam_data.items():

        # select reference transcript 
        ref_transcript = bam_info["reference_transcript"]

        # as a sanity check, get total cigar length and compare to read length. They should be the same. 
        read_length = bam_info.get("sequence_length", 0) 
        
        # get cigar lengths - note, 0 = match, 1 = ins, 2 = del, 3 = splice, 4 = softclip, 5 = hardclip 
        cigar_lengths = bam_info.get("cigar_lengths", {})
        cigar_values = "\t".join(str(cigar_lengths.get(i, 0)) for i in range(6))

        # get start and end positions
        start_position = bam_info.get("start_position", {})
        end_position = bam_info.get("end_position", {})

        # check if aligned length matches read length 
        #aligned_length = sum(cigar_lengths[0] + cigar_lengths[1] + cigar_lengths[3] + cigar_lengths[4] + cigar_lengths[5] + 1)
        aligned_length = cigar_lengths.get(0, 0) + cigar_lengths.get(1, 0) + cigar_lengths.get(3, 0) + cigar_lengths.get(4, 0) + cigar_lengths.get(5, 0) 

        # compare the CIGAR length and sequenced read length
        if aligned_length != read_length:
            cigar_discrepancy_reads.append(read_id) # if there's an error, write the read_id to our discrepency list 
            cigar_discrepancy_dict[read_id] = aligned_length - read_length # note the size of the error as well 

        # finally, combine everything into merged_data
        merged_data += f'{sample_name}\t{read_id}\t{ref_transcript}\t{read_length}\t{cigar_values}\t{start_position}\t{end_position}\n'


    # similarly, log read/transcript pairs with errors merging against the sequencing summary 

    # record cases where we have any discrepency in the cigar
    if cigar_discrepancy_reads:
        with open(os.path.join(os.path.dirname(args.output_file), '_cigar_discrepancy_errors.txt'), 'w') as f:
            for read_id, discrepancy in cigar_discrepancy_dict.items():
                f.write(f'{read_id}\t{discrepancy}\n')
        logging.error(f'Error in cigar lengths for {len(cigar_discrepancy_reads)} reads. Discrepancies written to cigar_discrepancy_errors.txt in the output directory.')
    
    # relay stats about the merging -- currently output by default, but could be moved to verbosity == 1
    logging.info(f'Total Reads Merged: {len(bam_data)}')
    logging.info(f'Total Transcripts Merged: {len(merged_data)}')

    # print 10 lines of merged data if verbosity == 2
    if args.verbosity == 2:
        preview_data = merged_data.splitlines()[:10]
        preview_data_str = '\n'.join(preview_data)
        logging.info("First 10 lines of the merged data:\n%s", preview_data_str)
    
    # write merged data to tmpfile
    temp_file_path = None
    temp_file_path = tmp_output_file + sample_name + "_merged_data.tsv"
    with open(temp_file_path, 'w') as file:
        file.write(merged_data)
    logging.info(f'Merged data written to temporary file: {temp_file_path}')

    # update user 
    elapsed_time = round(time.time() - start, 2)
    logging.info(f'Merged data in {elapsed_time} seconds')

    return merged_data, temp_file_path


BINSIZE = 10
def saturation(y):
    bins_len = max(y) +1 - BINSIZE
    freq1= np.zeros(max(y))
    for value in y:
         freq1[value-1] +=1
    if max(freq1)>len(y)/5:
         b = freq1[::-1]
         return len(b) - np.argmax(b)    
    freq_arr = np.zeros(bins_len)
    for value in y:
         for index in range(max(value-BINSIZE,1), min(value+1,bins_len)):
             freq_arr[index] +=1
    b = freq_arr[::-1]
    if np.argmax(b)>BINSIZE-1:
         return (len(b) - np.argmax(b)  +1)  
    else:     
         return len(b) - np.argmax(b) +BINSIZE - 1  


def saturation5p(y):
    bins_len = max(y) +1 - BINSIZE
    freq1= np.zeros(max(y)+1)
    for value in y:
         freq1[value] +=1
    if max(freq1)>max(len(y)/10,max(3,n)):
         return np.argmax(freq1)   
    return 0


# filter reads with inappropriate 3'end 
def MergeAll(merged_data_files):
    start = time.time()
    FinalData= pd.read_csv(merged_data_files[0], sep='\t')
    if n>1:
        for i in range(1,n):
            temp_data = pd.read_csv(merged_data_files[i], sep='\t')
            FinalData=pd.concat([FinalData,temp_data],ignore_index=True)
    return(FinalData)            



# filter reads with inappropriate 3'end 
def filter_reads_3pend_noadditional(merged_data, sample, output_path):
    start = time.time()
    combined_data = "sample_name\tread_id\tref_transcript\tstart_position\tend_position\tsaturation\tRealTranscriptLength\tnew_read_length\tfull_length\n"
    # read in the data 
    # group by ref_transcript
    grouped_by_transcript = merged_data.groupby("ref_transcript", group_keys=False)
    # Find 3'end saturation
    Start_sat = grouped_by_transcript["start_position"].apply(lambda x:saturation5p(x)).reset_index(name="Saturation_5p")  
    Start_sat = Start_sat.set_index('ref_transcript')  
    Saturation_df = grouped_by_transcript["end_position"].apply(lambda x:saturation(x)).reset_index(name="Saturation_3p")
    Saturation_df = Saturation_df.set_index('ref_transcript') 
    for index, read_info in merged_data.iterrows():
        # select reference transcript 
        ref_transcript = read_info["ref_transcript"]
        sample_name = read_info["sample_name"]
        sat=Saturation_df.loc[ref_transcript]
        init=Start_sat.loc[ref_transcript]
        sat = sat.min()
        init = init.min()
        if abs(read_info["end_position"]-sat)>50:
            continue 
        start = read_info["start_position"] 
        end = read_info["end_position"] 
        new_length = sat - max(start, init) 
        if new_length < Censor:
            continue
        full_length = abs(init - max(start, init)) <15
        rd_id = read_info["read_id"]
        TL = sat - init
        #combined_data += f'{read_info["read_id"]}\t{read_info["ref_transcript"]}\t{read_info["sequence_length"]}\t{read_info["cigar_match"]}\t{read_info["cigar_ins"]}\t{read_info["cigar_del"]}\t{read_info["cigar_splice"]}\t{read_info["cigar_softclip_total"]}\t{read_info["cigar_hardclip_total"]}\t{read_info["start_position"]}\t{read_info["end_position"]}\t{read_info["end_reason"]}\t{read_info["transcript_annotated_length"]}\t{read_info["transcript_biotype"]}\t{sat}\t{new_length}\t{full_length}\n'
        combined_data += f'{sample_name}\t{rd_id}\t{ref_transcript}\t{start}\t{end}\t{sat}\t{TL}\t{new_length}\t{full_length}\n'
    logging.info(f"Writing the processed data to {output_path + sample}_filtered_data.tsv...")
    file_path = output_path + sample + "_filtered_data.tsv"
    with open(file_path, 'w') as file:
        file.write(combined_data)
    # update user 
    elapsed_time = round(time.time() - start, 2)
    logging.info(f'Filtered data in {elapsed_time} seconds')
    return combined_data, file_path





# estimate fragmentation rate, convert to DTI, and test for random fragmentation
def Estimate_fragmentation_nogtf(filtered_data, output_path):
    start = time.time()
    grouped_by_sample = filtered_data.groupby("sample_name", group_keys=False) 
    Sample_DTI=[]
    for sampleName, info in grouped_by_sample:       
        logging.info(f"Computing fragmentation for sample {sampleName}")
        final_data = "transcript\tread_count\tNon_full_length_reads\ttotal_read_length\tfragmentation_rate\tDTI\tpvalue\tsaturation\tRealTranscriptLength\n"
        
        # group by ref_transcript
        grouped_by_transcript = info.groupby("ref_transcript", group_keys=False)
    
        # Find 3'end saturation
        Mean_df = grouped_by_transcript["Weighted_read_length"].sum().reset_index(name="Sum_length")
        FL_df = grouped_by_transcript["FullLength_weight"].sum().reset_index(name="Weight_number_Full_length")
        count_reads_df = grouped_by_transcript.size().reset_index(name="number_of_reads")
        
        # merge the three dataframes, add useful information
        temp_df = pd.merge(Mean_df, FL_df, on="ref_transcript")    
        result_df = pd.merge(temp_df, count_reads_df, on="ref_transcript")
        result_df["Saturation"] = grouped_by_transcript["saturation"].mean().reset_index(name="Saturation")["Saturation"]   
        result_df["RealTranscriptLength"] = grouped_by_transcript["RealTranscriptLength"].mean().reset_index(name="RealTranscriptLength")["RealTranscriptLength"]        
            
        # estimate fragmentation
        frag_rate = (result_df["number_of_reads"]-result_df["Weight_number_Full_length"])/(result_df["Sum_length"]-result_df["Weight_number_Full_length"])
        result_df['Fragmentation_rate'] = frag_rate
        result_df['DTI'] = list(map(DTI_Function, frag_rate)) 
        result_df = result_df.set_index('ref_transcript')
        Sample_DTI.append(round(result_df['DTI'].median(),2))
        
        # test for random fragmentation    
        for transcript, read_info in grouped_by_transcript:
            Result_info = result_df.loc[transcript]
            NbReads = int(Result_info["number_of_reads"])
            Sat =  int(Result_info["Saturation"])
            TL =  int(Result_info["RealTranscriptLength"])
            NbNonFullLength = int(NbReads-Result_info["Weight_number_Full_length"])
            SumLength = int(Result_info["Sum_length"]- NbReads)
            if NbReads<5:
                continue
            Kappa = Result_info["Fragmentation_rate"]
            DTI = Result_info["DTI"]
            y= read_info["Weighted_read_length"]
            y = np.array([TL if math.isnan(x) else x for x in y]).astype(int)
            Nbb=max(10,len(y)//10)
            pval = RandomTest(y,TL,Kappa,Nbb)
            final_data += f'{transcript}\t{NbReads}\t{NbNonFullLength}\t{SumLength}\t{Kappa}\t{DTI}\t{pval}\t{Sat}\t{TL}\n'
        # update user 
        logging.info(f"Writing the processed data to {output_path + sampleName}_process.txt...")
        with open(output_path + sampleName + "_process.txt", 'w') as file:
            file.write(final_data)
    DTI_df=pd.DataFrame({'Sample': sample_names,'DTI': Sample_DTI})
    DTI_df.to_csv(output_path + 'DTI' + sample + '.csv', index=False)
    elapsed_time = round(time.time() - start, 2)
    logging.info(f'Computed fragmentation rates of all samples in {elapsed_time} seconds')



def RandomTest(y,Sa,kappa,Nb): # Sa = saturation of 3'end value, kappa = fragmentation rate
    shift = (Sa - Censor) % (Nb)
    bins_size = max(1,(Sa - Censor) // (Nb))
    Emp_arr = np.zeros(Nb)
    for value in y:
         index = (value - Censor - shift) // bins_size
         Emp_arr[max(0,min(index,Nb - 1))] +=1         
    Theo_arr = np.zeros(Nb) 
    Theo_arr[0]+= 1-(1-kappa)**(bins_size+shift)
    for i in range(1,Nb):
         Theo_arr[i] += (1-kappa)**(i*bins_size+shift)*(1-(1-kappa)**(bins_size))
    Theo_arr[Nb-1] += (1-kappa)**(Sa-Censor-1)
    Theo_arr=Theo_arr/sum(Theo_arr)
    chi_square_test_statistic, p_value = stats.chisquare(Emp_arr, Theo_arr*len(y))
    return(p_value)

    
def process_vector(binary_vector):
    n = len(binary_vector)
    result_vector = [0] * n
    i1 = -1  # Initialize the first index of the block of 0
    i=1
    result_vector[0] = 1
    while i<n:
        if binary_vector[i-1] == 1:
            result_vector[i] = 1
            if binary_vector[i] == 0:
                a = i  # Update the first index of the block of 0
                i=i+1
                while binary_vector[i] == 0:
                    result_vector[i] = 1
                    i=i+1
                b=i-1
                result_vector[i] = (n - a) / (n - b-1) 
                i=i+1  
            else:
                i=i+1                
        else: 
            result_vector[i] = 1
            i=i+1
    return result_vector

    

def mask(input_path, censor):
    start = time.time()
    logging.info(f"Reading data from {input_path}...")
    
    # read and sort by transcript and md_read_end
    merged_data = pd.read_csv(input_path, sep='\t').sort_values(['ref_transcript', 'new_read_length'])
    logging.info("Data sorted.")
       
    # initialize new column
    merged_data['censoring_weight'] = [1]*len(merged_data['new_read_length'])

    # now: we use a nested censoring function to perform censoring
    def censoring_function_censor(group):
        # produce binary vector: 1 if read not censored, 0 otherwise
        binary_censor = group['end_reason'] != 'unblock_mux_change'
        group['censoring_weight'] = process_vector(binary_censor)
        return group

    # apply censoring function to each transcript group
    if censor:
        merged_data = merged_data.groupby('ref_transcript').apply(censoring_function_censor)    
    merged_data['FullLength_weight'] = merged_data['censoring_weight']*merged_data['full_length']
    merged_data['Weighted_read_length'] = merged_data['censoring_weight']*merged_data['new_read_length']
    return merged_data

def subset_bam(sample_name, input_bam, read_names_file):
    read_names = set()
    with open(read_names_file, "r") as file:
        # Iterate through each line in the file
        for line in file:
            # Split the line into columns based on a delimiter (e.g., tab)
            columns = line.strip().split("\t")
            # Check if the first column matches the given string
            if columns[0] == sample_name and len(columns) > 1:
                # Print the second column
                read_names.add(columns[1])
                
    # Open the input BAM file
    bam_in = pysam.AlignmentFile(input_bam, "rb")
    output_bam = args.output_file + sample_name + "_cleaned.bam"
    
    # Create a new BAM file to write the subset of reads
    bam_out = pysam.AlignmentFile(output_bam, "wb", template=bam_in)

    # Fetch reads by name and write them to the new BAM file
    for read in bam_in:
        if read.query_name in read_names:
            bam_out.write(read)

    # Close the input and output BAM files
    bam_in.close()
    bam_out.close()


# main 
def main(args):
    
    # logging verbosity 
    if args.verbosity == 1:
        logging.basicConfig(level=logging.INFO)
    elif args.verbosity >= 2:
        logging.basicConfig(level=logging.DEBUG)

    mergeFiles=[]
    for i in range(n):
        logging.info(f'Processing BAM file of sample {sample_names[i]}')
        bam_data = parse_bam(bam_files[i])
        merged_data, merge_temp_file = merge_data_justbam(bam_data,sample_names[i])
        mergeFiles.append(merge_temp_file)
    Merged_data_All=MergeAll(mergeFiles)
    logging.info('Filtering on 3prime end')
    filtered_data, filtered_file = filter_reads_3pend_noadditional(Merged_data_All,sample,args.output_file)
    censored_data = mask(filtered_file, False)                
    logging.info('Computing fragmentation and testing for random fragmentation')
    Estimate_fragmentation_nogtf(censored_data,args.output_file)
    
    # create subsetted bam files, if clean_bam == true 
    if args.clean_bam:
        for i in range(n):
            logging.info(f'Creating BAM file of sample {sample_names[i]}')
            subset_bam(sample_names[i], bam_files[i], filtered_file)

    # delete temp files, if keep == false (default)
    if not args.keep_temp:
        for temp_file in [tmp_output_file + sample_name + "_merged_data.tsv" for sample_name in sample_names]:
            if temp_file is not None and os.path.exists(temp_file):
                os.remove(temp_file)

if __name__ == "__main__":

    # parse arguments
    parser = argparse.ArgumentParser(description='Process GTF, BAM, and sequencing summary files.')
    parser.add_argument('--bam_file', help='Input sorted + indexed bam file.')
    parser.add_argument('--gtf_file', help='Input gtf/gff2 annnotation.')
    parser.add_argument('--summary_file', help='Input sequencing summary file from guppy.')
    parser.add_argument('--output_file', help='Output table.')
    parser.add_argument('--Condition',type=str, default='', help='Name of the group of samples.')    
    parser.add_argument('--samples',type=str, default='', help='comma delimited names of the samples.')     
    parser.add_argument('--ins', type=int, default=80, help='Maximum insertion length filtering') 
    parser.add_argument('--del', type=int, default=160, help='Maximum deletion length filtering')  
    parser.add_argument('--sc', type=int, default=200, help='Maximum soft-clip length filtering')     
    parser.add_argument('-t', '--threads', type=int, default=1, help='Number of threads to use (not implemented).') ## not implemented 
    parser.add_argument('-v', '--verbosity', type=int, default=0, help='Verbosity: 0 = Minimum, 1 = Information, 2 = Debugging.')
    parser.add_argument('-k', '--keep_temp', action='store_true', default=False, help='Keep temporary files in output directory.')
    parser.add_argument('-c', '--clean_bam', action='store_true', default=False, help='Produce a bam file containing only reads that were not discarded.')
        
    # flags for mask mode 
    parser.add_argument('--input_file', help='Input text file for mask mode.')
    parser.add_argument('--output_file_mask', help='Output text file for mask mode.')

    args = parser.parse_args()

    # validate args
    if not args.bam_file:
        parser.error("Degradation estimation requires at least a bam_file.")
    # adjust temp directory
    tmp_output_file = args.output_file+ "/tmp/"
    Path(tmp_output_file).mkdir(parents=True, exist_ok=True)
    sample = args.Condition
    sample_names =  args.samples.split(',')
    bam_files =  args.bam_file.split(',')
    if len(bam_files)!=len(sample_names):
        parser.error("Different sample numbers and bam files were provided.")
    n=len(bam_files)
    os.environ['TMPDIR'] = tmp_output_file

    main(args)
