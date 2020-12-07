#! /bin/bash

data="/ifb/data/mydatalocal"
mkdir -p $data
cd $data
# Quick script to clean extensions 
mkdir -p STAR_data_sorted
cd STAR_data_sorted

for file in *.fastq.gzAligned.out.bam.sorted; do
    mv "$file" "$(basename "$file" .fastq.gzAligned.out.bam.sorted).bam"
done
