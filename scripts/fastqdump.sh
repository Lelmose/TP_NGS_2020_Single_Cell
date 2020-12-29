#! /bin/bash

data="/ifb/data/mydatalocal"

# Create data and sra_data directories if not existing
mkdir -p $data
cd $data

mkdir -p sra_data
cd sra_data

# Create each line of SRR corresponds to an accession number to get on NCBI
SRR=`cat /ifb/data/mydatalocal/SRR_Acc_List.txt`

# For each element in our SRR file use the accession number to download the corresponding data
for srr in $SRR
do

# echo $srr
fastq-dump $srr --gzip
done