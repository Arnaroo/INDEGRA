# INDEGRA: INtegrity and DEGradation of RNA Analysis 
Evaluation of transcriptome-wide RNA degradation from long-read sequencing 


## About


------------------------------------------
# Table of Contents
------------------------------------------

   * [Dependencies](#dependencies)
   * [Installation](#installation)
   * [Preprocessing](#preprocess-drs-signals)
     *   [Basecalling](#basecalling)   
     *   [Alignment](#alignment)
   * [Estimate Degradation](#estimate-degradation)
     * [Transcriptome-wide RNA degradation evaluation](#transcriptome-wide-rna-degradation-evaluation)

------------------------------------------
# Dependencies
------------------------------------------
```
python=3.9
numpy==1.19.2
pysam==0.22.0
scipy==1.11.4
```

------------------------------------------
# Installation
------------------------------------------

Simply clone from github:
```
git clone https://github.com/Arnaroo/INDEGRA/ && cd INDEGRA
```


------------------------------------------
# Preprocess DRS signals
------------------------------------------

INDEGRA was tested on data basecalled with guppy 6.4.6 and aligned with minimap 2.24

## Basecalling

Recommended parameters:
```
guppy_basecaller -i $INPUTDIR --recursive -s $output_path -c guppy/ont-guppy/data/rna_r9.4.1_70bps_hac.cfg --device cuda:all:100% --compress_fastq --gpu_runners_per_device 2
```


## Alignment
minimap 2.24 for alignment and samtools 1.12 for quality checks

Recommended parameters:
-k 14 for human transcriptomes 

```
minimap2 -ax map-ont -k 14 ${fasta} ${input_path}/guppy.fastq | samtools sort -o ${output_path}.bam
samtools index ${output_path}.bam

samtools view -b -F 2324  ${output_path}.bam > ${output_path}_filtered.bam
samtools index ${output_path}_filtered.bam
```


------------------------------------------
# Estimate Degradation
------------------------------------------

## Transcriptome-wide RNA degradation evaluation


Example bash code to run INDEGRA Direct Transcript Integrity estimation.

```
export BamFiles="bam_filtered1.bam,bam_filtered2.bam" # comma-separated list of path to bam files to co-evaluate
export Condition="PatientSamples" # name for the group of bam files to be co-evaluated
export Samples="Patient1,Patient2" # comma-separated list of sample names corresponding to the bam files
export out="/home" # output directory 
export OUT=Hek293_mRNA_pU
export GTF="Homo_sapiens.GRCh38.109.gtf"     #for optional transcript annotation
export SUMMARY="sequencing_summary_run1,sequencing_summary_run2"	#for optional use of censoring information from sequencer 

python3 ./INDEGRA_scripts/INDEGRA.py --bam_file ${BamFiles} --Condition ${Condition} --samples ${Samples} --output_file "${out}/" [--gtf_file ${GTF}] [--summary_file ${SUMMARY}]
```


### Example of output files

INDEGRA produces several temporary and final output files.

The "sample_process.txt" files contain 7 to 9 columns indicating transcript name; read count of the transcript, estimated fragmentation rate, DTI estimation, p-value of random fragmentation test, 3' end estimation of the transcript, transcript length estimation, transcript annotated length and transcript biotype. The last two columns are only provided if a gtf file is used:

```
transcript	read_count	fragmentation_rate	DTI	pvalue	saturation	RealTranscriptLength
ENST00000000233	19	0.0008488410055501142	6.662195770179479	0.889744117123946	1032	1032
ENST00000001008	8	0.0005374899220639613	8.112481522335202	0.5978979204727686	2234	2234
ENST00000002165	10	0.0008426966292134832	6.682818211723543	0.09573402350945918	1710	1710
ENST00000003100	6	0.0004384426517011575	8.67865805911778	0.49690385624058464	3155	3155
ENST00000007390	12	0.0008776128929312271	6.568526696075884	0.5528488484255603	1209	1209
```

The "condition_DTI.txt" file contains 2 columns indicating sample name and sample DTI estimation:
```
sample	DTI
sample1	8.1
sample2	8.4
```
