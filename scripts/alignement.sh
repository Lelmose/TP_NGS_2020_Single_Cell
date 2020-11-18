#! /bin/bash

data="/ifb/data/mydatalocal"
mkdir -p $data
cd $data

# Create a new directory for aligned data
mkdir -p alignement_data

# Create a list of all the downloaded cells
ALIGN=`ls /home/rstudio/data/mydatalocal/trimmomatic_data`

# For each element in our liqt run the Java package Trimmomatic
for align in $ALIGN
do
startdir=/ifb/data/mydatalocal/sra_data/$align
targetdir=/ifb/data/mydatalocal/alignement_data/$align
salmon quant -i /ifb/data/mydatalocal/Ref_seq/mouse_index_unzipped -l SR -r $startdir --validateMappings -o $targetdir --gcBias
done