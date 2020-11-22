#! /bin/bash

data="/ifb/data/mydatalocal"
mkdir -p $data
cd $data
mkdir -p sra_data
cd sra_data

SRR=`cat /ifb/data/mydatalocal/SRR_Acc_List.txt`

# For each element in our SRR file use the accession number to download the corresponding data
for srr in $SRR
do

# echo $srr
fastq-dump $srr --gzip
done