#!/bin/bash
usage="$(basename "$0") [-h] [-o j t f v s m g i l] -- Visualization tool for MEI supporting reads, e.g.,

where:
    -h  show this help text
    -o  output folder path
    -j  subject ID
    -t  TE class: LINE/ALU
    -f  TE family: L1HS/AluYa5
    -v  version control for RetroSom (default 1)
    -s  strandness of the MEI (+ strand: 1, - strand: 0)
    -m  masterpath (default ~/masterpath)
    -g  huamn reference genome (hg19/hg38)
    -i  input file
    -l  label of images"

while getopts ":ho:j:t:f:v:s:m:g:i:l:" opt; do
  case $opt in
    h) echo "$usage"
       exit
       ;;
    o) outpath=$(readlink -f $OPTARG)
       ;;
    j) sub="$OPTARG"
       ;;
    t) TEclass="$OPTARG"
       ;;
    f) TEfamily="$OPTARG"
       ;;
    v) ver="$OPTARG"
       ;;
    s) strand="$OPTARG"
       ;;
    m) masterpath=$(readlink -f $OPTARG)
       ;;
    g) hg="$OPTARG"
       ;;
    i) file="$OPTARG"
       ;;
    l) label="$OPTARG"
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

while read line; do
         var1=$(echo $line|awk '{print $1}')
         var2=$(echo $line|awk '{print $2}')
         var3=$(echo $line|awk '{print $3}')
         echo "$var1 $var2 $var3"
         ./DeepVis.sh -i $sub -t $TEclass -f $TEfamily -c $var1 -d $var2 -e $var3 -r $ver -s $strand -p $outpath -m $masterpath -g $hg -l $label -a $file
done < ./$file