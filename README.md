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
   * [Differential Transcript Abundance](#differential-transcript-abundance)
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

Example bash code to run INDEGRA Direct Transcript Integrity estimation on the provided test datasets.

```
export BamFiles="Test_Data/Undegraded_Rep1_chr1_reads.bam,Test_Data/Undegraded_Rep2_chr1_reads.bam,Test_Data/Deg400s_Rep1_chr1_reads.bam,Test_Data/Deg400s_Rep2_chr1_reads.bam" # comma-separated list of path to bam files to co-evaluate
export Condition="Test_Data" # name for the group of bam files to be co-evaluated
export Samples="Undegraded_Rep1,Undegraded_Rep2,Deg400s_Rep1,Deg400s_Rep2" # comma-separated list of sample names corresponding to the bam files
export out="./Output" # output directory 

mkdir -p ${out} 2>/dev/null
python3 ./INDEGRA_scripts/INDEGRA.py --bam_file ${BamFiles} --Condition ${Condition} --samples ${Samples} --output_file "${out}/" 
```

Several additional options can be used:
    -k (equivalent to --keep_temp) to keep temporary files in the /tmp folder, 
    -c (equivalent to --clean_bam) to produce a bam files containing only reads that were not discarded during any of the INDEGRA steps,
    --ins (default 80) to change the maximum insertion length filtering,
    --deletion (default 160) to change the maximum deletion length filtering,
    --sc (default 200) to change the maximum softclip length filtering, and
    -s (equivalent to --readlength, default 150) to change the threshold to correct on the shortest reads.

### Example of output files

INDEGRA produces several temporary and final output files.

The "*_process.txt" files contain 9 columns indicating transcript name; read count of the transcript, number of reads that are not full-length, sum of all read length, estimated fragmentation rate, DTI estimation, p-value of random fragmentation test, 3' end estimation of the transcript, and transcript length estimation.

```
transcript	read_count	Non_full_length_reads	total_read_length	fragmentation_rate	DTI	pvalue	saturation	RealTranscriptLength
ENST00000003912.7	7	7	3253	0.002147239263803681	4.557750349397514	0.9988469099341366	5481	5481
ENST00000040877.2	10	10	7646	0.0013061650992685476	5.563516290521175	0.9992865862137844	5184	5184
ENST00000054666.11	39	38	22515	0.0016849199663016006	5.017229430306395	0.2724024593898182	2178	2178
ENST00000060969.6	11	11	7353	0.0014937533948940792	5.266958301258209	0.9566233947303409	5494	5494
ENST00000194214.10	29	29	13482	0.0021463992302568277	4.558447975239521	0.002348012622075157	685	685

```

The "DTI.csv" file contains 2 columns indicating sample name and sample DTI estimation:
```
Sample,DTI
Undegraded_Rep1,4.78
Undegraded_Rep2,4.82
Deg400s_Rep1,7.03
Deg400s_Rep2,5.77
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


process1="Undegraded_Rep1_process.txt"
process2="Deg400s_Rep1_process.txt"

p=0.1
Result<-Test_Degradation(process1,process2,p)
PlotResults(Result)
```

By default, the code will subsample the transcripts with more than 2000 reads. This default can be modified with Test_Degradation(process1,process2,p,thresh=5000).
It is possible to change the prior probability without re-running the whole pipeline with the following command:

```
Result2=Change_Prior(Result,0.05)
PlotResults(Result2, samplenames=c("Undegraded","Deg400s"), labels=FALSE, GeneLabel=c("# INDEGRA: INtegrity and DEGradation of RNA Analysis 
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
   * [Differential Transcript Abundance](#differential-transcript-abundance)
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

Example bash code to run INDEGRA Direct Transcript Integrity estimation on the provided test datasets.

```
export BamFiles="Test_Data/Undegraded_Rep1_chr1_reads.bam,Test_Data/Undegraded_Rep2_chr1_reads.bam,Test_Data/Deg400s_Rep1_chr1_reads.bam,Test_Data/Deg400s_Rep2_chr1_reads.bam" # comma-separated list of path to bam files to co-evaluate
export Condition="Test_Data" # name for the group of bam files to be co-evaluated
export Samples="Undegraded_Rep1,Undegraded_Rep2,Deg400s_Rep1,Deg400s_Rep2" # comma-separated list of sample names corresponding to the bam files
export out="./Output" # output directory 

mkdir -p ${out} 2>/dev/null
python3 ./INDEGRA_scripts/INDEGRA.py --bam_file ${BamFiles} --Condition ${Condition} --samples ${Samples} --output_file "${out}/" 
```

Several additional options can be used:
    -k (equivalent to --keep_temp) to keep temporary files in the /tmp folder, 
    -c (equivalent to --clean_bam) to produce a bam files containing only reads that were not discarded during any of the INDEGRA steps,
    --ins (default 80) to change the maximum insertion length filtering,
    --deletion (default 160) to change the maximum deletion length filtering,
    --sc (default 200) to change the maximum softclip length filtering, and
    -s (equivalent to --readlength, default 150) to change the threshold to correct on the shortest reads.

### Example of output files

INDEGRA produces several temporary and final output files.

The "*_process.txt" files contain 9 columns indicating transcript name; read count of the transcript, number of reads that are not full-length, sum of all read length, estimated fragmentation rate, DTI estimation, p-value of random fragmentation test, 3' end estimation of the transcript, and transcript length estimation.

```
transcript	read_count	Non_full_length_reads	total_read_length	fragmentation_rate	DTI	pvalue	saturation	RealTranscriptLength
ENST00000003912.7	7	7	3253	0.002147239263803681	4.557750349397514	0.9988469099341366	5481	5481
ENST00000040877.2	10	10	7646	0.0013061650992685476	5.563516290521175	0.9992865862137844	5184	5184
ENST00000054666.11	39	38	22515	0.0016849199663016006	5.017229430306395	0.2724024593898182	2178	2178
ENST00000060969.6	11	11	7353	0.0014937533948940792	5.266958301258209	0.9566233947303409	5494	5494
ENST00000194214.10	29	29	13482	0.0021463992302568277	4.558447975239521	0.002348012622075157	685	685

```

The "DTI.csv" file contains 2 columns indicating sample name and sample DTI estimation:
```
Sample,DTI
Undegraded_Rep1,4.78
Undegraded_Rep2,4.82
Deg400s_Rep1,7.03
Deg400s_Rep2,5.77
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


process1="Undegraded_Rep1_process.txt"
process2="Deg400s_Rep1_process.txt"

p=0.1
Result<-Test_Degradation(process1,process2,p)
PlotResults(Result)
```

By default, the code will subsample the transcripts with more than 2000 reads. This default can be modified with Test_Degradation(process1,process2,p,thresh=5000).
It is possible to change the prior probability without re-running the whole pipeline with the following command:

```
Result2=Change_Prior(Result,0.05)
PlotResults(Result2, samplenames=c("Undegraded","Deg400s"), labels=FALSE, GeneLabel=c("# INDEGRA: INtegrity and DEGradation of RNA Analysis 
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
   * [Differential Transcript Abundance](#differential-transcript-abundance)
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

Example bash code to run INDEGRA Direct Transcript Integrity estimation on the provided test datasets.

```
export BamFiles="Test_Data/Undegraded_Rep1_chr1_reads.bam,Test_Data/Undegraded_Rep2_chr1_reads.bam,Test_Data/Deg400s_Rep1_chr1_reads.bam,Test_Data/Deg400s_Rep2_chr1_reads.bam" # comma-separated list of path to bam files to co-evaluate
export Condition="Test_Data" # name for the group of bam files to be co-evaluated
export Samples="Undegraded_Rep1,Undegraded_Rep2,Deg400s_Rep1,Deg400s_Rep2" # comma-separated list of sample names corresponding to the bam files
export out="./Output" # output directory 

mkdir -p ${out} 2>/dev/null
python3 ./INDEGRA_scripts/INDEGRA.py --bam_file ${BamFiles} --Condition ${Condition} --samples ${Samples} --output_file "${out}/" 
```

Several additional options can be used:
    -k (equivalent to --keep_temp) to keep temporary files in the /tmp folder, 
    -c (equivalent to --clean_bam) to produce a bam files containing only reads that were not discarded during any of the INDEGRA steps,
    --ins (default 80) to change the maximum insertion length filtering,
    --deletion (default 160) to change the maximum deletion length filtering,
    --sc (default 200) to change the maximum softclip length filtering, and
    -s (equivalent to --readlength, default 150) to change the threshold to correct on the shortest reads.

### Example of output files

INDEGRA produces several temporary and final output files.

The "*_process.txt" files contain 9 columns indicating transcript name; read count of the transcript, number of reads that are not full-length, sum of all read length, estimated fragmentation rate, DTI estimation, p-value of random fragmentation test, 3' end estimation of the transcript, and transcript length estimation.

```
transcript	read_count	Non_full_length_reads	total_read_length	fragmentation_rate	DTI	pvalue	saturation	RealTranscriptLength
ENST00000003912.7	7	7	3253	0.002147239263803681	4.557750349397514	0.9988469099341366	5481	5481
ENST00000040877.2	10	10	7646	0.0013061650992685476	5.563516290521175	0.9992865862137844	5184	5184
ENST00000054666.11	39	38	22515	0.0016849199663016006	5.017229430306395	0.2724024593898182	2178	2178
ENST00000060969.6	11	11	7353	0.0014937533948940792	5.266958301258209	0.9566233947303409	5494	5494
ENST00000194214.10	29	29	13482	0.0021463992302568277	4.558447975239521	0.002348012622075157	685	685

```

The "DTI.csv" file contains 2 columns indicating sample name and sample DTI estimation:
```
Sample,DTI
Undegraded_Rep1,4.78
Undegraded_Rep2,4.82
Deg400s_Rep1,7.03
Deg400s_Rep2,5.77
```

Additionnaly, if option -c (equivalent to --clean_bam) is used, INDEGRA will produce a bam file per sample containing the reads retained through the pipeline.

------------------------------------------
# Differential Biological Degradation
------------------------------------------

Running the pipeline simply requires providing the "_process.txt" files and choosing a prior probability p of two transcripts having different biological degradation rates. This  parameter p controls the false positive rate. For instance p can be chosen as 0.5 in the case of no prior expectation on the different degradation rates in two samples, or as 0.05 when very few transcripts are expected to differ.
```
library(dplyr)
library(ggpubr)
library(ggrepel)
source("INDEGRA_scripts/Functions_Degradation.R")


process1="Output/Undegraded_Rep1_process.txt"
process2="Output/Deg400s_Rep1_process.txt"

p=0.1
Result<-Test_Degradation(process1,process2,p)
PlotResults(Result)
```

By default, the code will subsample the transcripts with more than 2000 reads. This default can be modified with Test_Degradation(process1,process2,p,thresh=5000).
It is possible to change the prior probability without re-running the whole pipeline with the following command:

```
Result2=Change_Prior(Result,0.05)
PlotResults(Result2, samplenames=c("Undegraded","Deg400s"), labels=FALSE, GeneLabel=c("ENST00000003912.7","ENST00000054666.11"))
```

The GeneLabel option allows to visualize transcripts of particular interest.
A summary of the significant hits can be obtained with the Get_Significant function, providing read counts in each sample, biological degradation rate estimates, log-fold change and posterior probability of difference in rates:
```
Sig=Get_Significant(Result2)
```



------------------------------------------
# Differential Transcript Abundance
------------------------------------------

### Computing normalisation offsets


Correction of degradation bias relies on Lowess regression of log-counts to DTI values. Offset matrices can then be provided to standard differential expression tools such as DESeq2 or edgeR. 
Running the pipeline simply requires providing the sample names and their "_process.txt" files. At this stage information on groups is not required
```
library(ggpubr)
library(aroma.light)
source("INDEGRA_scripts/Functions_DTE.R")

samples=c("Undegraded_Rep1","Undegraded_Rep2","Deg400s_Rep1","Deg400s_Rep2")
samplefiles=c("Output/Undegraded_Rep1_process.txt", "Output/Undegraded_Rep2_process.txt", "Output/Deg400s_Rep1_process.txt", "Output/Deg400s_Rep2_process.txt") 

A=Normalize_DTI(samples,samplefiles)
Plot_Normalisation(A)

```

### Example usage with DESeq2

RawCounts together with normalization offset should be provided. Those are saved in the Normalize_DTI output:

```
library(DESeq2)

Group=c('Undegraded','Undegraded','Degraded','Degraded')
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

Group=c('Undegraded','Undegraded','Degraded','Degraded')
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




