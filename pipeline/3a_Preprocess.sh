#!/bin/bash

usage="$(basename "$0") [-h] [-o j m v g] -- Preprocessing before calling the putative MEIs. Requires 1 core, and 20gb memory; Requires samtools, bedtools, bcftools and exonerate.

where:
    -h  show this help text
    -o  output folder path
    -j  subject ID
    -m  masterpath
    -v  version control for RetroSom (default 3)
    -g  reference genome (default hg38, supporting hg38, hg19 and b37)"

while getopts ":ho:j:m:v:c:g:" opt; do
  case $opt in
    h) echo "$usage"
       exit
       ;;
    o) outpath=$(readlink -f $OPTARG)
       ;;
    j) sub="$OPTARG"
       ;;
    m) masterpath=$(readlink -f $OPTARG)
       ;;
    v) ver="$OPTARG"
       ;;
    g) hg="$OPTARG"
       ;;
    :) printf "missing argument for -%s\n" "$OPTARG" >&2
       echo "$usage" >&2
       exit 1
       ;;
    \?) printf "illegal option: -%s\n" "$OPTARG" >&2
       echo "$usage" >&2
       exit 1
       ;;
  esac
done

# Activate RetroSom virtual env #
source /opt/miniconda3/bin/activate RetroSom

#
retro=retro_v$ver
workfolder=$outpath/$sub
echo "$sub ran on $(hostname -s)"

date '+%m/%d/%y %H:%M:%S'
echo "TE calling phase ... begins"

### Group1 => 5, convert part of PE supporting reads to SR ###
$masterpath/utls/00_GroupItoV.pl $masterpath $sub $workfolder/$retro
### select a window around supporting reads ###
$masterpath/utls/01_correct.bed.pl $masterpath $sub 300 $workfolder $retro $hg
### check the sequencing depth at the selected windows, require bedtools ### 
$masterpath/utls/02_depth_zero.sh $sub $workfolder $retro
date '+%m/%d/%y %H:%M:%S'
echo "Depth analysis ... finishes"