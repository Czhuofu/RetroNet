#!/bin/bash

usage="$(basename "$0") [-h] [-o j m v] -- Predict the confidence of each supporting read to be real MEI using RF, NB and LR

where:
    -h  show this help text
    -o  output folder path
    -j  subject ID
    -m  masterpath
    -v  version control for RetroSom (default 3)"

while getopts ":ho:j:m:v:" opt; do
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

echo "$outpath ran on $(hostname -s)"

Rscript $masterpath/ALU/02_PE_level1/01.RF.LogR.NB.r $outpath $ver 1 $sub $masterpath
$masterpath/ALU/02_PE_level1/02_combine_RF.pl $outpath $ver 1 $sub
$masterpath/ALU/02_PE_level1/03_3probs.pl $outpath $ver 1 $sub

Rscript $masterpath/ALU/02_PE_level1/01.RF.LogR.NB.r $outpath $ver 0 $sub $masterpath
$masterpath/ALU/02_PE_level1/02_combine_RF.pl $outpath $ver 0 $sub
$masterpath/ALU/02_PE_level1/03_3probs.pl $outpath $ver 0 $sub