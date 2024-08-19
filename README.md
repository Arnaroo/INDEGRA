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
     *   [Transcriptome-wide RNA degradation evaluation](#transcriptome-wide-rna-degradation-evaluation)
     *   [Example of output files](#example-of-output-files)
   * [Differential Biological Degradation](#differential-biological-degradation)
   * [Differential Transcript Expression](#differential-transcript-expression)
     *   [Computing normalisation offsets](#computing-normalisation-offsets)
     *   [Example usage with DESeq2](#example-usage-with-deseq2)
     *   [Example usage with edgeR](#example-usage-with-edgeR)

------------------------------------------
# Dependencies
------------------------------------------
Python library dependencies for DTI metric computations
```
python=3.9
numpy==1.19.2
pysam==0.22.0
scipy==1.11.4
```

R package dependencies for differential degradation and expression
```
R=4.2
library(dplyr)
library(ggpubr)
library(aroma.light)
library(DESeq2)
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

INDEGRA was tested on data basecalled with guppy 6.4.6 and aligned with minimap 2.24. 
Its imput is a bam file, necessarily aligned to a transcriptome reference.

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

To allow more robust estimation and comparison of degradation rates in different samples, we recommend running INDEGRA on all samples simultaneously. This will still return one file per sample containing transcript- and sample-specific degradation rates, as well as a summary file containing DTI metrics for each sample. Running all samples simultaneously allows to obtain a better estimate of the transcripts 3' and 5' ends, but does not limit analysis to transcripts shared by all samples.

Example bash code to run INDEGRA Direct Transcript Integrity estimation.

```
export BamFiles="bam_filtered1.bam,bam_filtered2.bam" # comma-separated list of path to bam files to co-evaluate
export Condition="PatientSamples" # name for the group of bam files to be co-evaluated
export Samples="Patient1,Patient2" # comma-separated list of sample names corresponding to the bam files
export out="/home" # output directory 

python3 ./INDEGRA_scripts/INDEGRA.py --bam_file ${BamFiles} --Condition ${Condition} --samples ${Samples} --output_file "${out}/" 
```

Two additional options can be used, -k (equivalent to --keep_temp) to keep temporary files in the /tmp folder, and -c (equivalent to --clean_bam) to produce a bam files containing only reads that were not discarded during any of the INDEGRA steps. Note that this option is not recommended by default as it significantly slows down the runtime of the pipeline.

### Example of output files

INDEGRA produces several temporary and final output files.

The "sample_process.txt" files contain 9 columns indicating transcript name; read count of the transcript, number of reads that are not full-length, sum of all read length, estimated fragmentation rate, DTI estimation, p-value of random fragmentation test, 3' end estimation of the transcript, and transcript length estimation.

```
transcript	read_count	Non_full_length_reads	total_read_length	fragmentation_rate	DTI	pvalue	saturation	RealTranscriptLength
ENSMUST00000000090	7	4	4409	0.0009064128710627691	6.4792676098960555	0.9418699497384452	645	645
ENSMUST00000000756	90	36	60074	0.0005989020129762103	7.741778551478916	0.11054252722032465	780	704
ENSMUST00000001452	8	8	14315	0.0005585422048453537	7.98124920533038	0.49786479239934633	1977	1977
ENSMUST00000001460	6	6	7555	0.0007935458272715249	6.856310131199388	0.18660788759013464	3474	3474
ENSMUST00000001513	5	5	6293	0.0007939028262940616	6.854992216232237	0.20950722331285973	1758	1758

```

The "DTI.csv" file contains 2 columns indicating sample name and sample DTI estimation:
```
Sample,DTI
Patient1,5.85
Patient2,6.10
```

Additionnaly, if option -c (equivalent to --clean_bam) is used, INDEGRA will produce a bam file per sample containing the reads retained through the pipeline.

------------------------------------------
# Differential Biological Degradation
------------------------------------------

Running the pipeline simply requires providing the "_process.txt" files and choosing a prior probability p of two transcripts having different biological degradation rates. This  parameter p controls the false positive rate. For instance p can be chosen as 0.5 in the case of no prior expectation on the different degradation rates in two samples, or as 0.05 when very few transcripts are expected to differ.
```
library(dplyr)
library(ggpubr)
source("INDEGRA_scripts/Functions_Degradation.R")


process1="Patient1_process.txt"
process2="Patient2_process.txt"

p=0.5
Result<-Test_Degradation(process1,process2,p)
PlotResults(Result)
```

It is possible to change the prior probability without re-running the whole pipeline with the following command:

```
Result2=Change_Prior(Result,0.05)
PlotResults(Result2, samplenames=c("Patient 1","Patient 2"), labels=TRUE)
```

A summary of the significant hits can be obtained with the Get_Significant function, providing read counts in each sample, biological degradation rate estimates, log-fold change and posterior probability of difference in rates:
```
Sig=Get_Significant(Result2)
```



------------------------------------------
# Differential Transcript Expression
------------------------------------------

### Computing normalisation offsets


Correction of degradation bias relies on Lowess regression of log-counts to DTI values. Offset matrices can then be provided to standard differential expression tools such as DESeq2 or edgeR. 
Running the pipeline simply requires providing the sample names and their "_process.txt" files. At this stage information on groups is not required
```
library(ggpubr)
library(aroma.light)
source("INDEGRA_scripts/Functions_DTE.R")

samples=c("Sample1","Sample2","Sample3","Sample4")
samplefiles=c("PathTo/Sample1_process.txt", "PathTo/Sample2_process.txt", "PathTo/Sample3_process.txt", "PathTo/Sample4_process.txt") 

A=Normalize_DTI(samples,samplefiles)
Plot_Normalisation(A)

```

### Example usage with DESeq2

RawCounts together with normalization offset should be provided. Those are saved in the Normalize_DTI output:

```
library(DESeq2)

Group=c('Group1','Group1','Group2','Group2')
coldata=data.frame(Condition=factor(Group))
rownames(coldata)=colnames(A$RawCounts)

normFactors <- exp(-1 * A$MatrixOffset)
normFactors <- normFactors / exp(rowMeans(log(normFactors)))



ddsOffset <-DESeqDataSetFromMatrix(countData = A$RawCounts,
                              colData = coldata,
                              design = ~ Condition)
normalizationFactors(ddsOffset) <- normFactors
ddsOffset <- DESeq(ddsOffset)
resOffset <- results(ddsOffset)
resOffset

```

### Example usage with edgeR

RawCounts together with normalization offset should be provided. Those are saved in the Normalize_DTI output:

```
library(edgeR)

Group=c('Group1','Group1','Group2','Group2')
coldata=data.frame(Condition=factor(Group))
rownames(coldata)=colnames(A$RawCounts)

design <- model.matrix(~Condition, data=coldata)

yOff <- DGEList(counts=A$RawCounts,
             group=coldata$Condition)
yOff$offset <- -A$MatrixOffset
yOff <- estimateDisp(yOff, design)
fitoff <- glmFit(yOff, design)
lrtoff <- glmLRT(fitoff, coef=2)
resOff=topTags(lrtoff,n=nrow(A$RawCounts))
```

