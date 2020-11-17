#! /bin/bash

data="/ifb/data/mydatalocal"
mkdir -p $data
cd $data

# Create a new directory for fastqc_data
mkdir -p fastqc_data


# Create a new directory for fastqc_data
FASTQC=`find /ifb/data/mydatalocal/sra_data/. -type f | shuf -n 100`

# For each element in our SRR file use the accession number to download the corresponding data
for fastqc in $FASTQC
do
fastqc $fastqc -o ./fastqc_data -t 7
echo $fastqc
# fastq-dump $srr --gzip
done