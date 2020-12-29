# Re-analysis of data from 'Dental cell type atlas reveals stem and differentiated cell types in mouse and human teeth'

## Introduction

Human dental tissue can't regenarate and this statement has forced human to use artificial implants. However, those implants are not always tolerated by human body and there are numbers of associated pathologies. One different approach to this issue is the study of rodent dental epithelium. Due to they way of living, rodents' incisive don't stop growing during there life. This implicates that some stem cell niches remain active during their whole life. Following the framework of the paper from Krikanez et al. we reanalyzed a dataset of single cell RNA seq extracted from mouse dental gyrus. The aim of this work was to identificate the different cellular population using differential expression and to associate each population to the corresponding cellular type by studying genetic markers. Then the final objective was to approach the temporal dynamic of the different cell type with the method of RNA velocity.

## Getting the data
The data were collected on NCBI (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE146123). The dataset used in this work was obtained by selecting an homogeneous subset of the original study i.e. : Mus musculus, healthy, incisor, smart seq2.
The accession numbers of the cells that fullfil all the criteria were sum up in an SSR accession number list and the script `fastqdump.sh` permitted to get the files with **fastq-dump**.

## Evaluation of the quality

In order to evaluate the quality of sequencing the sc

## Purification

## First alignment using Salmon

## Analisis of the data with Seurat

### Estimating variance of the Dataset
### Dimensional reduction
### Clustering
### Annotation of the clusters

## Toward RNA velocity 1: Alignement on the whole genome

## Toward RNA velocity 2: Sorting the counts and creating the loom file

## Toward RNA velocity 3: Vizualizing velocity and identifying genes prone to explain it

## Conclusion


