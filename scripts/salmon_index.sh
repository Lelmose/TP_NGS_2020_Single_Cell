#!/bin/bash

data="/home/rstudio/data/mydatalocal"
mkdir -p $data
cd $data

mkdir -p Ref_seq
cd Ref_seq

# Salmon index is build from the reference sequence bank: Mus_musculus.GRCm38.cdna.all.fa 
# and the annotations: mouse_index_unzipped
salmon index -t Mus_musculus.GRCm38.cdna.all.fa -i mouse_index_unzipped  -k 31 