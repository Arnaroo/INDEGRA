#!/usr/bin/env python3

"""\
Basic python3 implementation of Nanograd4
Written by AS on 2023-06-29
aditya.sethi@anu.edu.au
"""

import argparse
import gffutils
import logging
import os
import pandas as pd
import numpy as np 
import pysam
import tempfile
import time
from tqdm import tqdm
import math
from math import *
import scipy.stats as stats


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
            if cigar_maxs[1]>80:
                continue
            # filter if max deletion >160
            if cigar_maxs[2]>160:
                continue
            # filter if soft clip >200
            if cigar_maxs[4]>200:
                continue

            # strip version number from transcript id
            reference_transcript = read.reference_name.split('.')[0]

            # write this data into our bam dict
            bam_data[read.query_name] = {"sequence_length": seq_len, "cigar_lengths": cigar_lengths, "reference_transcript": reference_transcript, "start_position": read.reference_start, "end_position": read.reference_end}

        elapsed_time = round(time.time() - start, 2)
        logging.info(f'Parsed {len(bam_data)} lines from BAM file in {elapsed_time} seconds')

        # print a preview of bam_data if verbosity == 2
        if args.verbosity == 2:
            header = "\n".join([f"{key}: {value}" for key, value in list(bam_data.items())[:5]])
            logging.info(f"BAM Data Header:\n{header}")

        # write the bam dict to a temporary file, if the keep flag is enabled 
        temp_file_path = None
        if args.keep_temp:
            temp_file_path = tempfile.mktemp(suffix=sample + "_parsed_bam.tsv", dir=os.path.dirname(os.path.abspath(args.output_file)))
            df = pd.DataFrame.from_dict(bam_data, orient="index")
            df.to_csv(temp_file_path, sep="\t", index=True)

    return bam_data, temp_file_path 

# function to filter gtf for transcripts seen in the bam (this speeds up the gtf parsing later on)
def filter_gtf(bam_data, gtf_file_path):
    start = time.time()

    # get transcript_ids of interest from the BAM we just parsed 
    transcript_ids = set([info["reference_transcript"] for info in bam_data.values()])
    temp_file_path = tempfile.mktemp(dir=os.path.dirname(os.path.abspath(gtf_file_path)))

    lines = []
    with open(gtf_file_path, 'r') as gtf_file:

        # keep cases where the transcript_id in the gtf file coincides with a transcript_id in the bam file 
        lines = [line for line in gtf_file if not line.startswith('#') and line.split('\t')[8].split(';')[2].split('"')[1].split('.')[0] in transcript_ids]

    # record how many entries we keep 
    original_line_count = sum(1 for line in open(gtf_file_path, 'r') if not line.startswith('#'))
    filtered_line_count = len(lines)

    # write the filtered lines to the temp file
    with open(temp_file_path, 'w') as temp_file:
        temp_file.writelines(lines)

    # update user 
    elapsed_time = round(time.time() - start, 2)
    logging.info(f'Filtered {filtered_line_count} GTF features from {original_line_count} GTF features in {elapsed_time} seconds')

    # print header of parsed gtf if verbosity == 2 
    if args.verbosity == 2:
        with open(temp_file_path, 'r') as temp_file:
            lines = [next(temp_file) for x in range(10)]
            logging.info("First 10 lines of the filtered GTF file:\n" + "\n".join(lines))

    return temp_file_path

# function to parse gtf :: returns, transcript_id (w/o version), transcript_biotype, transcript_length 
def parse_gtf(gtf_file_path):
    start = time.time()

    # use gffutils to parse the gtf file 
    # note, there are small differences between ensembl, gencode, other gtf files 
    # this will need to be tested and refined depending on what annotation is being used 
    gtf_db = gffutils.create_db(gtf_file_path, dbfn='gtf.db', force=True, keep_order=True, merge_strategy='merge', sort_attribute_values=True, disable_infer_genes = True, disable_infer_transcripts = True)
    
    # create a dict for transcripts 
    gtf_data = {}
    
    for feature in gtf_db.all_features(featuretype='transcript'):
        
        # remove transcript_version from transcript_id (otherwise, there can be discrecpency with the bam)
        transcript_id = feature.id.split(".")[0]  # remove version

        # fetch biotype or return NA
        # important: GENCODE uses "transcript_type" whereas other annotations use different field names to store the transcript biotype information 
        # I have seen before, transcript_biotype, gene_type etc. 
        biotype = feature.attributes['transcript_type'][0] if 'transcript_type' in feature.attributes else 'NA'
        
        # get transcript length 
        # note, we calculate trancsript length from the sum of the cumulative exons 
        # otherwise, we get the unspliced length
        length = sum(len(i) for i in gtf_db.children(feature, featuretype='exon'))  

        # store the data in a dict 
        gtf_data[transcript_id] = {"biotype": biotype, "length": length}

    elapsed_time = round(time.time() - start, 2)
    logging.info(f'Parsed {len(gtf_data)} lines from GTF file in {elapsed_time} seconds')
    
    # print 10 lines if verbosity == 2
    if args.verbosity == 2:
        preview_data = dict(list(gtf_data.items())[:10])
        logging.info(f"First 10 lines of the GTF data:\n{preview_data}")

    # write the bam dict to a temporary file, if the keep flag is enabled 
    temp_file_path = None
    if args.keep_temp:
        temp_file_path = tempfile.mktemp(suffix=sample + "_parsed_gtf.tsv", dir=os.path.dirname(os.path.abspath(args.output_file)))
        df = pd.DataFrame.from_dict(gtf_data, orient="index")
        df.to_csv(temp_file_path, sep="\t", index=True)
 
    return gtf_data, temp_file_path

# function to parse sequencing_summary file :: returns, read_id, read_end_reason 
def parse_summary(summary_file_path):
    start = time.time()

    # read in the data and select the read_id and read_end_reason 
    summary_df = pd.read_csv(summary_file_path, sep='\t')
    summary_data = {row["read_id"]: row["end_reason"] for _, row in summary_df.iterrows()}
    
    elapsed_time = round(time.time() - start, 2)
    logging.info(f'Parsed {len(summary_data)} lines from summary file in {elapsed_time} seconds')
    
    # print 10 lines if verbosity == 2
    if args.verbosity == 2:
        preview_data = dict(list(summary_data.items())[:10])
        logging.info(f"First 10 lines of the summary data:\n{preview_data}")

    # write the bam dict to a temporary file, if the keep flag is enabled 
    temp_file_path = None
    if args.keep_temp:
        temp_file_path = tempfile.mktemp(suffix=sample + "_parsed_summary.tsv", dir=os.path.dirname(os.path.abspath(args.output_file)))
        df = pd.DataFrame.from_dict(summary_data, orient="index")
        df.to_csv(temp_file_path, sep="\t", index=True)

    return summary_data, temp_file_path

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
    if args.keep_temp:
        temp_file_path = tempfile.mktemp(suffix=sample + "_merged_data.tsv", dir=os.path.dirname(os.path.abspath(args.output_file)))
        with open(temp_file_path, 'w') as file:
            file.write(merged_data)
    logging.info(f'Merged data written to temporary file: {temp_file_path}')

    # update user 
    elapsed_time = round(time.time() - start, 2)
    logging.info(f'Merged data in {elapsed_time} seconds')

    return merged_data, temp_file_path

# function to merge data from bam, gtf, sequencing_summary, using transcript_id as key 
# includes the MDI_read_length, the sum of cigar match, deletion, insertion 
def merge_data_gtf_noseq(gtf_data, bam_data,sample_name):
    start = time.time()

    # header for merged data 
    merged_data = "sample_name\tread_id\tref_transcript\tsequence_length\tcigar_match\tcigar_ins\tcigar_del\tcigar_splice\tcigar_softclip_total\tcigar_hardclip_total\tstart_position\tend_position\ttranscript_annotated_length\ttranscript_biotype\n"

    # make list to store any transcripts which have merging error with gtf data 
    gtf_error_transcripts = []
    
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
            
        else:
            gtf_error_transcripts.append(ref_transcript)

        # finally, combine everything into merged_data
        merged_data += f'{sample_name}\t{read_id}\t{ref_transcript}\t{read_length}\t{cigar_values}\t{start_position}\t{end_position}\t{gtf_info["length"]}\t{gtf_info["biotype"]}\n'


    # log read/transcripts pairs that have errors merging with the gtf 
    # note, we get a few in the test data, because the gtf is chromosome only whereas the FASTA that was used has scaffolds, alt contigs etc. 
    if gtf_error_transcripts:
        unique_gtf_error_transcripts = set(gtf_error_transcripts)
        with open(os.path.join(os.path.dirname(args.output_file), 'gtf_merge_errors.txt'), 'w') as f:
            f.write("\n".join(unique_gtf_error_transcripts))
        logging.error(f'Error merging data for {len(gtf_error_transcripts)} reads with gtf data. Error transcripts written to gtf_merge_errors.txt in the output directory.')
    
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

    # print 10 lines of merged data if verbosity == 2
    if args.verbosity == 2:
        preview_data = merged_data.splitlines()[:10]
        preview_data_str = '\n'.join(preview_data)
        logging.info("First 10 lines of the merged data:\n%s", preview_data_str)
    
    # write merged data to tmpfile
    temp_file_path = None
    if args.keep_temp:
        temp_file_path = tempfile.mktemp(suffix=sample + "_merged_data.tsv", dir=os.path.dirname(os.path.abspath(args.output_file)))
        with open(temp_file_path, 'w') as file:
            file.write(merged_data)
    logging.info(f'Merged data written to temporary file: {temp_file_path}')

    # update user 
    elapsed_time = round(time.time() - start, 2)
    logging.info(f'Merged data in {elapsed_time} seconds')

    return merged_data, temp_file_path

# function to merge data from bam, gtf, sequencing_summary, using transcript_id as key 
# includes the MDI_read_length, the sum of cigar match, deletion, insertion 
def merge_data_nogtf(bam_data, summary_data,sample_name):
    start = time.time()

    # header for merged data 
    merged_data = "sample_name\tread_id\tref_transcript\tsequence_length\tcigar_match\tcigar_ins\tcigar_del\tcigar_splice\tcigar_softclip_total\tcigar_hardclip_total\tstart_position\tend_position\tend_reason\n"
  
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


            # now, check if the read id is in the sequencing_summary file 
        if read_id in summary_data:
            end_reason = summary_data[read_id] # fetch the read_end_reason            
        else:
            ss_error_transcripts.append(ref_transcript)

        # finally, combine everything into merged_data
        merged_data += f'{sample_name}\t{read_id}\t{ref_transcript}\t{read_length}\t{cigar_values}\t{start_position}\t{end_position}\t{end_reason}\n'


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
    logging.info(f'Total Transcripts Merged: {len(merged_data)}')
    logging.info(f'Reads in BAM which do not match sequencing summary: {len(ss_error_transcripts)}, transcripts these reads map to: {len(set(ss_error_transcripts))}')

    # print 10 lines of merged data if verbosity == 2
    if args.verbosity == 2:
        preview_data = merged_data.splitlines()[:10]
        preview_data_str = '\n'.join(preview_data)
        logging.info("First 10 lines of the merged data:\n%s", preview_data_str)
    
    # write merged data to tmpfile
    temp_file_path = None
    if args.keep_temp:
        temp_file_path = tempfile.mktemp(suffix=sample + "_merged_data.tsv", dir=os.path.dirname(os.path.abspath(args.output_file)))
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
    if args.keep_temp:
        temp_file_path = tempfile.mktemp(suffix=sample + "_merged_data.tsv", dir=os.path.dirname(os.path.abspath(args.output_file)))
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
def filter_reads_3pend(merged_data):
    start = time.time()
    combined_data = "sample_name\tread_id\tref_transcript\ttranscript_annotated_length\ttranscript_biotype\tstart_position\tend_position\tend_reason\tsaturation\tRealTranscriptLength\tnew_read_length\tfull_length\n"
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
        end_reason = read_info["end_reason"]
        ref_len =  read_info["transcript_annotated_length"]
        biotype = read_info["transcript_biotype"]
        TL = sat - init
        #combined_data += f'{read_info["read_id"]}\t{read_info["ref_transcript"]}\t{read_info["sequence_length"]}\t{read_info["cigar_match"]}\t{read_info["cigar_ins"]}\t{read_info["cigar_del"]}\t{read_info["cigar_splice"]}\t{read_info["cigar_softclip_total"]}\t{read_info["cigar_hardclip_total"]}\t{read_info["start_position"]}\t{read_info["end_position"]}\t{read_info["end_reason"]}\t{read_info["transcript_annotated_length"]}\t{read_info["transcript_biotype"]}\t{sat}\t{new_length}\t{full_length}\n'
        combined_data += f'{sample_name}\t{rd_id}\t{ref_transcript}\t{ref_len}\t{biotype}\t{start}\t{end}\t{end_reason}\t{sat}\t{TL}\t{new_length}\t{full_length}\n'
    # write combined data to tmpfile
    temp_file_path = None
    if args.keep_temp:
        temp_file_path = tempfile.mktemp(suffix=sample + "_filtered_data.tsv", dir=os.path.dirname(os.path.abspath(args.output_file)))
        with open(temp_file_path, 'w') as file:
            file.write(combined_data)
    logging.info(f'Filtered data written to temporary file: {temp_file_path}')
    # update user 
    elapsed_time = round(time.time() - start, 2)
    logging.info(f'Filtered data in {elapsed_time} seconds')
    return combined_data, temp_file_path

# filter reads with inappropriate 3'end 
def filter_reads_3pend_nogtf(merged_data):
    start = time.time()
    combined_data = "sample_name\tread_id\tref_transcript\tstart_position\tend_position\tend_reason\tsaturation\tRealTranscriptLength\tnew_read_length\tfull_length\n"
    # read in the data 
    # group by ref_transcript
    grouped_by_transcript = merged_data.groupby("ref_transcript", group_keys=False)
    # Find 3'end saturation
    Saturation_df = grouped_by_transcript["end_position"].apply(lambda x:saturation(x)).reset_index(name="Saturation_3p")
    Saturation_df = Saturation_df.set_index('ref_transcript') 
    Start_sat = grouped_by_transcript["start_position"].apply(lambda x:saturation5p(x)).reset_index(name="Saturation_5p")  
    Start_sat = Start_sat.set_index('ref_transcript')  
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
        end_reason = read_info["end_reason"]
        TL = sat - init
        #combined_data += f'{read_info["read_id"]}\t{read_info["ref_transcript"]}\t{read_info["sequence_length"]}\t{read_info["cigar_match"]}\t{read_info["cigar_ins"]}\t{read_info["cigar_del"]}\t{read_info["cigar_splice"]}\t{read_info["cigar_softclip_total"]}\t{read_info["cigar_hardclip_total"]}\t{read_info["start_position"]}\t{read_info["end_position"]}\t{read_info["end_reason"]}\t{read_info["transcript_annotated_length"]}\t{read_info["transcript_biotype"]}\t{sat}\t{new_length}\t{full_length}\n'
        combined_data += f'{sample_name}\t{rd_id}\t{ref_transcript}\t{start}\t{end}\t{end_reason}\t{sat}\t{TL}\t{new_length}\t{full_length}\n'
    # write combined data to tmpfile
    temp_file_path = None
    if args.keep_temp:
        temp_file_path = tempfile.mktemp(suffix=sample + "_filtered_data.tsv", dir=os.path.dirname(os.path.abspath(args.output_file)))
        with open(temp_file_path, 'w') as file:
            file.write(combined_data)
    logging.info(f'Filtered data written to temporary file: {temp_file_path}')
    # update user 
    elapsed_time = round(time.time() - start, 2)
    logging.info(f'Filtered data in {elapsed_time} seconds')
    return combined_data, temp_file_path


# filter reads with inappropriate 3'end 
def filter_reads_3pend_nosummary(merged_data):
    start = time.time()
    combined_data = "sample_name\tread_id\tref_transcript\ttranscript_annotated_length\ttranscript_biotype\tstart_position\tend_position\tsaturation\tRealTranscriptLength\tnew_read_length\tfull_length\n"
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
        ref_len =  read_info["transcript_annotated_length"]
        biotype = read_info["transcript_biotype"]
        TL = sat - init
        #combined_data += f'{read_info["read_id"]}\t{read_info["ref_transcript"]}\t{read_info["sequence_length"]}\t{read_info["cigar_match"]}\t{read_info["cigar_ins"]}\t{read_info["cigar_del"]}\t{read_info["cigar_splice"]}\t{read_info["cigar_softclip_total"]}\t{read_info["cigar_hardclip_total"]}\t{read_info["start_position"]}\t{read_info["end_position"]}\t{read_info["end_reason"]}\t{read_info["transcript_annotated_length"]}\t{read_info["transcript_biotype"]}\t{sat}\t{new_length}\t{full_length}\n'
        combined_data += f'{sample_name}\t{rd_id}\t{ref_transcript}\t{ref_len}\t{biotype}\t{start}\t{end}\t{sat}\t{TL}\t{new_length}\t{full_length}\n'
    # write combined data to tmpfile
    temp_file_path = None
    if args.keep_temp:
        temp_file_path = tempfile.mktemp(suffix=sample + "_filtered_data.tsv", dir=os.path.dirname(os.path.abspath(args.output_file)))
        with open(temp_file_path, 'w') as file:
            file.write(combined_data)
    logging.info(f'Filtered data written to temporary file: {temp_file_path}')
    # update user 
    elapsed_time = round(time.time() - start, 2)
    logging.info(f'Filtered data in {elapsed_time} seconds')
    return combined_data, temp_file_path

# filter reads with inappropriate 3'end 
def filter_reads_3pend_noadditional(merged_data):
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
    # write combined data to tmpfile
    temp_file_path = None
    if args.keep_temp:
        temp_file_path = tempfile.mktemp(suffix=sample + "_filtered_data.tsv", dir=os.path.dirname(os.path.abspath(args.output_file)))
        with open(temp_file_path, 'w') as file:
            file.write(combined_data)
    logging.info(f'Filtered data written to temporary file: {temp_file_path}')
    # update user 
    elapsed_time = round(time.time() - start, 2)
    logging.info(f'Filtered data in {elapsed_time} seconds')
    return combined_data, temp_file_path



# estimate fragmentation rate, convert to DTI, and test for random fragmentation
def Estimate_fragmentation(filtered_data, output_path):
    start = time.time()
    grouped_by_sample = filtered_data.groupby("sample_name", group_keys=False) 
    for sampleName, info in grouped_by_sample:
        final_data = "transcript\tread_count\tfragmentation_rate\tDTI\tpvalue\tsaturation\tRealTranscriptLength\ttranscript_annotated_length\ttranscript_biotype\n"
        
        # group by ref_transcript
        grouped_by_transcript = filtered_data.groupby("ref_transcript", group_keys=False)
    
        # Find 3'end saturation
        Mean_df = grouped_by_transcript["Weighted_read_length"].sum().reset_index(name="Sum_length")
        FL_df = grouped_by_transcript["FullLength_weight"].sum().reset_index(name="Weight_number_Full_length")
        count_reads_df = grouped_by_transcript.size().reset_index(name="number_of_reads")
    
        # merge the three dataframes, add useful information
        temp_df = pd.merge(Mean_df, FL_df, on="ref_transcript")    
        result_df = pd.merge(temp_df, count_reads_df, on="ref_transcript")
        result_df["Annotated_length"] = grouped_by_transcript["transcript_annotated_length"].mean().reset_index(name="Annotated_Length")["Annotated_Length"] 
        result_df["transcript_biotype"] = grouped_by_transcript["transcript_biotype"].unique().reset_index(name="transcript_biotype")["transcript_biotype"]  
        result_df["Saturation"] = grouped_by_transcript["saturation"].mean().reset_index(name="Saturation")["Saturation"]  
        result_df["RealTranscriptLength"] = grouped_by_transcript["RealTranscriptLength"].mean().reset_index(name="RealTranscriptLength")["RealTranscriptLength"]        
        
        # estimate fragmentation
        frag_rate = (result_df["number_of_reads"]-result_df["Weight_number_Full_length"])/(result_df["Sum_length"]-result_df["Weight_number_Full_length"])
        result_df['Fragmentation_rate'] = frag_rate
        result_df['DTI'] = list(map(DTI_Function, frag_rate)) 
        result_df = result_df.set_index('ref_transcript')
        
        # test for random fragmentation
        for transcript, read_info in grouped_by_transcript:    
            Result_info = result_df.loc[transcript]  	
            Sat =  int(Result_info["Saturation"])
            TL =  int(Result_info["RealTranscriptLength"])
            NbReads = int(Result_info["number_of_reads"])
            Biotype = Result_info["transcript_biotype"]
            AnnoLen = Result_info["Annotated_length"]
            if NbReads<5:
                continue
            Kappa = Result_info["Fragmentation_rate"]
            DTI = Result_info["DTI"]
            y= read_info["Weighted_read_length"]
            y = np.array([TL if math.isnan(x) else x for x in y]).astype(int)
            Nbb=max(10,len(y)//5)
            pval = RandomTest(y,TL,Kappa,Nbb)
            final_data += f'{transcript}\t{NbReads}\t{Kappa}\t{DTI}\t{pval}\t{Sat}\t{TL}\t{AnnoLen}\t{Biotype}\n'
        # update user 
        logging.info(f"Writing the processed data to {output_path + sampleName}_process.txt...")
        with open(output_path + sampleName + "_process.txt", 'w') as file:
            file.write(final_data)
    elapsed_time = round(time.time() - start, 2)
    logging.info(f'Computed fragmentation rates of all samples in {elapsed_time} seconds')



# estimate fragmentation rate, convert to DTI, and test for random fragmentation
def Estimate_fragmentation_nogtf(filtered_data, output_path):
    start = time.time()
    grouped_by_sample = filtered_data.groupby("sample_name", group_keys=False) 
    for sampleName, info in grouped_by_sample:       
        logging.info(f"Computing fragmentation for sample {sampleName}")
        final_data = "transcript\tread_count\tfragmentation_rate\tDTI\tpvalue\tsaturation\tRealTranscriptLength\n"
        
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
        
        # test for random fragmentation    
        for transcript, read_info in grouped_by_transcript:
            Result_info = result_df.loc[transcript]
            Sat =  int(Result_info["Saturation"])
            TL =  int(Result_info["RealTranscriptLength"])
            NbReads = int(Result_info["number_of_reads"])
            if NbReads<5:
                continue
            Kappa = Result_info["Fragmentation_rate"]
            DTI = Result_info["DTI"]
            y= read_info["Weighted_read_length"]
            y = np.array([TL if math.isnan(x) else x for x in y]).astype(int)
            Nbb=max(10,len(y)//10)
            pval = RandomTest(y,TL,Kappa,Nbb)
            final_data += f'{transcript}\t{NbReads}\t{Kappa}\t{DTI}\t{pval}\t{Sat}\t{TL}\n'
        # update user 
        logging.info(f"Writing the processed data to {output_path + sampleName}_process.txt...")
        with open(output_path + sampleName + "_process.txt", 'w') as file:
            file.write(final_data)
    elapsed_time = round(time.time() - start, 2)
    logging.info(f'Computed fragmentation rates of all samples in {elapsed_time} seconds')

def RandomTest(y,Sa,kappa,Nb): # Sa = saturation of 3'end value, kappa = fragmentation rate
    shift = (Sa - Censor) % (Nb)
    bins_size = (Sa - Censor) // (Nb)
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
        bam_data, bam_temp_file = parse_bam(bam_files[i])
        merged_data, merge_temp_file = merge_data_justbam(bam_data,sample_names[i])
        mergeFiles.append(merge_temp_file)
    Merged_data_All=MergeAll(mergeFiles)
    logging.info('Filtering on 3prime end')
    filtered_data, filtered_temp_file = filter_reads_3pend_noadditional(Merged_data_All)
    censored_data = mask(filtered_temp_file, False)                
    logging.info('Computing fragmentation and testing for random fragmentation')
    Estimate_fragmentation_nogtf(censored_data,args.output_file)
    
    # delete temp files, if keep == false (default)
    if not args.keep_temp:
        for temp_file in [bam_temp_file, gtf_temp_file, summary_temp_file, merge_temp_file]:
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
    parser.add_argument('-t', '--threads', type=int, default=1, help='Number of threads to use (not implemented).') ## not implemented 
    parser.add_argument('-v', '--verbosity', type=int, default=0, help='Verbosity: 0 = Minimum, 1 = Information, 2 = Debugging.')
    parser.add_argument('-k', '--keep_temp', action='store_true', default=False, help='Keep temporary files in output directory.')
    
    # flags for mask mode 
    parser.add_argument('--input_file', help='Input text file for mask mode.')
    parser.add_argument('--output_file_mask', help='Output text file for mask mode.')

    args = parser.parse_args()

    # validate args
    if not args.bam_file:
        parser.error("Degradation estimation requires at least a bam_file.")
    # adjust temp directory
    tmp_output_file = args.output_file
    sample = args.Condition
    sample_names =  args.samples.split(',')
    bam_files =  args.bam_file.split(',')
    if len(bam_files)!=len(sample_names):
        parser.error("Different sample numbers and bam files were provided.")
    n=len(bam_files)
    os.environ['TMPDIR'] = os.path.dirname(os.path.abspath(tmp_output_file))

    main(args)
