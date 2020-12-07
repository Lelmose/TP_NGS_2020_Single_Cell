#! /bin/bash

data="/ifb/data/mydatalocal"
mkdir -p $data
cd $data

# Create a new directory for loom data
mkdir -p loom_files/MyTooth


startdir=/ifb/data/mydatalocal/STAR_data_sorted/

targetdir=/ifb/data/mydatalocal/loom_files/
refdir=/ifb/data/mydatalocal/STAR_index/Mus_musculus.GRCm38.101.gtf
/home/rstudio/.local/bin/velocyto run-smartseq2 -o $targetdir\
 /ifb/data/mydatalocal/STAR_data_sorted/*.bam $refdir

