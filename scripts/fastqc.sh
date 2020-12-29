#! /bin/bash

data="/ifb/data/mydatalocal"
mkdir -p $data
cd $data

# Create a new directory for fastqc_data if not existing
mkdir -p fastqc_data_trimed

# Selection randomly 10 cells in the dataset to run fastqc on them
FASTQC=`find /ifb/data/mydatalocal/trimmomatic_data/. -type f | shuf -n 10`

# For each cell produce a fastqc report of quality
for fastqc in $FASTQC
do

fastqc $fastqc -o ./fastqc_data_trimed 

done