#!/bin/bash

# 0) Path to Workfolder
Workfolder=$yourdownloadpath

# 1) Run Step 2 for HG00096 by using NA12878 as control
sub=HG00096
bash ${Workfolder}/RetroNet/Singularity_Slurm_RetroNet_step2.sh \
   -o ${Workfolder} \
   -j ${sub} \
   -m ${Workfolder}/RetroNet \
   -g hg38 \
   -p tiny \
   -c ${Workfolder}/RetroNet/NA12878_100X_hg38_filter2 \
   -x 0.95 \
   -l 1