#! /bin/bash

data="/ifb/data/mydatalocal"
mkdir -p $data
cd $data

# Create a new directory for STARed data
mkdir -p STAR_data

# Create a list of all the downloaded/purified cells
STAR=`ls /home/rstudio/data/mydatalocal/trimmomatic_data`

# As we treat a large amount of files (unusual) we will keep the index in memory to save time.
# We do so through STAR --genomeLoad
refdir=/ifb/data/mydatalocal/STAR_index/Indexes
STAR --genomeLoad LoadAndExit  --genomeDir $refdir

# For each element in our liqt run STAR GeneCounts
for star in $STAR
do
startdir=/ifb/data/mydatalocal/trimmomatic_data/$star
targetdir=/ifb/data/mydatalocal/STAR_data/$star 

# --readFilesCommand correspond to wether the input files are comressed or not. 
# zcat correspond to compressed files.
# --outSAMtype can be used to sort the counts in the same time they're generated. 
# In our cases there is too much files to do so. (We wiil do it later using samtoolsort.sh)

STAR --genomeDir $genomdir \
      --readFilesIn $startdir \
      --readFilesCommand zcat \
      --outSAMtype BAM unsorted \
      --quantMode GeneCounts \
      --outFileNamePrefix $targetdir
done

STAR --genomeLoad Remove  --genomeDir $refdir