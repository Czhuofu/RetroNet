#!/bin/bash
### calculate sequencing depth around supporting reads ###
### reads with abnormal depth will be removed later ###
sub=$1
workfolder=$2
retro=$3

### sort the bed file by coordiantes  ###
awk '{print $1"\t"$2"\t"$3}' \
   $workfolder/$retro/$sub.bed | \
   sort -k1,1 -k2,3n | sort -u \
   > $workfolder/$retro/$sub.uniq.bed \

### sequencing depth at each window ###
samtools bedcov \
   $workfolder/$retro/$sub.uniq.bed \
   $workfolder/align/$sub.final.bam \
   > $workfolder/$retro/$sub.depth
