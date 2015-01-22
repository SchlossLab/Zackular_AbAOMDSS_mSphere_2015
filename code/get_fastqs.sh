#!/bin/bash

################################################################################
#
# get_fastqs.sh
#
# Currently this pulls the fastq files from Joe's directory on Axiom. In the
# future, we'll have to replace the cp command below with a pull out of the SRA,
# convert to fastq, and decompression.
# 
#
# Dependencies...
# * The *.files file out of the data/process/ directory
#
# Output...
# * Puts the listed fastq files into the data/raw folder
#
################################################################################

# the name of the files file that will be used in make.contigs
FILES_FILE=$1   


# where are the raw fastqs?
ORIG_LOCATION=/mnt/EXT/Schloss-data/zackular/WT-Nod1-analyses


# where are we putting the fastqs?
RAW_FOLDER=data/raw


# let's clean the fastq files out of the raw folder
rm $RAW_FOLDER/*.fastq


# the second column of the files file contains the file names for read 1 and the
# third column contains the file names for read 2
FASTQ="$(cut -f 2 $FILES_FILE) $(cut -f 3 $FILES_FILE)" 


# cycle through the FASTQ file names and copy them from one place to the next
for FILE in $FASTQ
do
    echo $FILE
    cp "$ORIG_LOCATION/$FILE" "$RAW_FOLDER"
done

