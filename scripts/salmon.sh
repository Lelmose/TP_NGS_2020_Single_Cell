#!/bin/bash

data="/home/rstudio/data/mydatalocal"
mkdir -p $data

cd $data

mkdir -p Ref_seq
cd Ref_seq
salmon index -t Mus_musculus.GRCm38.cdna.all.fa -i mouse_index_unzipped  -k 31 