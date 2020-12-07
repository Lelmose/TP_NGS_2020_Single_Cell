#! /bin/bash

data="/ifb/data/mydatalocal"
mkdir -p $data
cd $data

# Create a new directory for STARed data
mkdir -p STAR_data_sorted

# Create a list of all the Aligned cells
SORT=`ls /home/rstudio/data/mydatalocal/STAR_data_bam`


# For each element in our list, run samtools sort
for sort in $SORT
do
startdir=/ifb/data/mydatalocal/STAR_data_bam/$sort
targetdir=/ifb/data/mydatalocal/STAR_data_sorted/$sort 
samtools sort $startdir -o $targetdir.sorted
done
