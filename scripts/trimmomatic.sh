#! /bin/bash

data="/ifb/data/mydatalocal"
mkdir -p $data
cd $data

# Create a new directory for purified data
mkdir -p trimmomatic_data

# Create a list of all the downloaded cells
TRIM=`ls /home/rstudio/data/mydatalocal/sra_data`

# For each element in our list run the Java package Trimmomatic
for trim in $TRIM
do

# Define the different directory
startdir=/ifb/data/mydatalocal/sra_data/$trim
targetdir=/ifb/data/mydatalocal/trimmomatic_data/$trim

java -jar  /softwares/Trimmomatic-0.39/trimmomatic-0.39.jar SE $startdir $targetdir \
ILLUMINACLIP:/softwares/Trimmomatic-0.39/adapters/TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

done
