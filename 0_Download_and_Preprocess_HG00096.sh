#!/bin/bash

# 0) Path to Workfolder (Replace with your downlaod path)
Workfolder=$yourdownloadpath

# 1) Downlad HG00096 WGS from HGSVC
cd ${Workfolder}
wget ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR324/ERR3240114/HG00096.final.cram
# md5sum d3354f61a055adfcfc988470bc507b2d

# 2) Convert CRAM to BAM file
mkdir HG00096
samtools view -bh -F 0x0100 -F 0x400 -F 0x800 \
    -T ./GRCh38DH/GRCh38_full_analysis_set_plus_decoy_hla.fa \
    HG00096.final.cram \
    -o ${Workfolder}/HG00096/HG00096.clean.bam
samtools index ${Workfolder}/HG00096/HG00096.clean.bam