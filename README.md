# Re-analysis of data from 'Dental cell type atlas reveals stem and differentiated cell types in mouse and human teeth'

## Introduction

Human dental tissue can't regenarate and this statement has forced human to use artificial implants. However, those implants are not always tolerated by human body and there are numbers of associated pathologies. One different approach to this issue is the study of rodent dental epithelium. Due to they way of living, rodents' incisive don't stop growing during there life. This implicates that some stem cell niches remain active during their whole life. Following the framework of the paper from Krikanez et al. we reanalyzed a dataset of single cell RNA seq extracted from mouse dental gyrus. The aim of this work was to identificate the different cellular population using differential expression and to associate each population to the corresponding cellular type by studying genetic markers. Then the final objective was to approach the temporal dynamic of the different cell type with the method of RNA velocity.

## Getting the data
The data were collected on NCBI (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE146123). The dataset used in this work was obtained by selecting an homogeneous subset of the original study i.e. : Mus musculus, healthy, incisor, smart seq2. The corresponding subset included 2555 cells.

The accession numbers of the cells that fullfil all the criteria were sum up in an SSR accession number list and the script `fastqdump.sh` permitted to get the files with **fastq-dump**.

## Evaluation of the quality

In order to evaluate the quality of sequencing the **fastqc** tool was run on a subset of our dataset with the script `fastqc.sh`. Then a **multiqc** overview of the data permitted to identify that the quality of the fragment was not always statisfying on the ends of the reads. No significant contamination nor sequencing error were reported.

## Purification

In order to purify the dataset from the detected bias, **trimmomatic** was run with `trimmomatic.sh`. 

The targeted bias were:

* Remove the adaptators
* Remove poor quality bases
* Remove poor quality reads
* Remove too short reads

A second **fastqc** analysis on the trimmed data confirmed the efficiency of the purification as illustrated in the following figures.

![Mean per base quality](Pictures/fastqc_per_base_sequence_quality_plot.png)

The poor quality extremities have been removed.

## First alignment using Salmon

The first idea for identifying cellular types based on RNAseq data was to align the reads on a reference mouse transcriptome in order to count the number of reads per genes and to be able to go further in comparisons. This alignement on transcriptome was performed with **salmon** in a two steps framework.
First a salmon index was obtained with the script `salmon_index.sh` using as reference `Mus_musculus.GRCm38.cdna.all.fa` and `mouse_index_unzipped`.

Then, the Java package **salmon** xas run on all the cells with `alignment.sh`

## Analysis of the data with Seurat

The R package Seurat permits to conduct statistical analysis on single cell data. 

### Purification on gene counts

With the precise number of counts per genes a more precise purification can be run on the dataset. The three following values were used as reference of sequencing quality and a threshold was set on those values to ensure the consistence of the data.

![Data before purification](Pictures/before_purif.png) 

* Total number of reads: remove the lowest and highest 5%
* Number of different genes: remove the lowest 5%
* Percent of motochondrial RNA: only cells with less than 15% of mitochondrial RNA were kept.

1962 cells corresponded to these criteria and the Fig testifies that extreme values have been removed.

![Data after purification](Pictures/Features_after_purif.png)

### Estimating variance of the Dataset and rescaling the data

As the following analysis aim to detect different cell type it is important to select genes that exhibit variability in their expression. Based on the Fig, the 5% more expressed genes were selected. This step is important because it the first dimensional reduction of our data and it will prevent some bias for the next dimensional reductions. A lot of poorly variable genes will have a deleterious impact on the PCA because it will increase sparsity and introduice correlated variables that are not consistent with the goal of this work.

### Dimensional reduction
After the previous purification steps the dataset was still high dimension (5548 genes on 1962 cells) and then non usable for clustering algorithms. For this reason a PCA method was run on the data using Seurat. 
As it will be discussed in the part about clustering, PCA presents some weakness in terms of display (the clusters tend to overlap and for the machine it is not an issue but it make the figure looks strange)
### Clustering
### Annotation of the clusters

## Toward RNA velocity 1: Alignement on the whole genome

## Toward RNA velocity 2: Sorting the counts and creating the loom file

## Toward RNA velocity 3: Vizualizing velocity and identifying genes prone to explain it

## Conclusion


