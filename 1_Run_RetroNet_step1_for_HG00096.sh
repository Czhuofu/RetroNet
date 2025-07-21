#!/bin/bash

# 0) Path to Workfolder
Workfolder=/home/zhinanlin/scratch

# 1) Run Step 1 for HG00096
sub=HG00096
bash ${Workfolder}/RetroNet/Singularity_Slurm_RetroNet_step1.sh \
   -o ${Workfolder} \
   -j ${sub} \
   -m ${Workfolder}/RetroNet \
   -g hg38 \
   -p tiny \
   -i 1 \
   -b ${Workfolder}/HG00096/HG00096.clean.bam \
   -n 100