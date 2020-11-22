#!/bin/bash

data="/home/rstudio/data/mydatalocal"
mkdir -p $data
cd $data
mkdir -p STAR_index
cd STAR_index
ls
# The runMode 'genomeGenerate' generate a STAR compatible genome file that we will 
# need to perform alignments
# --sjbdOverhang indicate what is the maximum possible stretch of sequence that can
# be found on one side of a spicing site. As our fragments are at most 43 we set 
# this value to 42.
STAR --runMode genomeGenerate --genomeDir Indexes \
--runThreadN 7 \
--genomeFastaFiles Mus_musculus.GRCm38.dna.primary_assembly.fa \
--sjdbGTFfile Mus_musculus.GRCm38.101.gtf \
--sjdbOverhang 42