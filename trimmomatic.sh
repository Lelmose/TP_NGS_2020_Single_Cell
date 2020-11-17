#! /bin/bash

data="/ifb/data/mydatalocal"
mkdir -p $data
cd $data

# Create a new directory for prufied data
mkdir -p trimmomatic_data

# Create a list of all the downloaded cells
TRIM=`ls /home/rstudio/data/mydatalocal/sra_data`

# For each element in our liqt run the Java package Trimmomatic
for trim in $TRIM
do
startdir=/ifb/data/mydatalocal/sra_data/$trim
targetdir=/ifb/data/mydatalocal/trimmomatic_data/$trim
java -jar  /softwares/Trimmomatic-0.39/trimmomatic-0.39.jar SE $startdir $targetdir \
ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
done