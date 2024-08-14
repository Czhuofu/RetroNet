#!/bin/bash
usage="$(basename "$0") [-h] [-o j t v m x] -- Visualization tool for MEI supporting reads, e.g.,

where:
    -h  show this help text
    -o  output folder path
    -j  subject ID
    -t  TE class: LINE/ALU
    -v  version control for RetroSom (default 1)
    -m  masterpath (default ~/masterpath)
    -g  hg(default hg38)
    -x  cutoff or threshold (default 0.5)"

while getopts ":ho:j:t:f:v:s:m:g:i:x:" opt; do
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
    v) ver="$OPTARG"
       ;;
    m) masterpath=$(readlink -f $OPTARG)
       ;;
    g) hg="$OPTARG"
       ;;
    x) cutoff="$OPTARG"
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

mkdir $outpath/$sub/RetroNet
python3 RetroNet.py $outpath $sub $TEclass $ver $masterpath $cutoff $hg
python3 $masterpath/RetroNet/pipeline/generate_bed.py $outpath $sub $TEclass $ver $masterpath $cutoff $hg
##Compress the png picture##
tar -zcvf $outpath/$sub/visual_$ver/$TEclass.tar.gz $outpath/$sub/visual_$ver/$TEclass/  --remove-files
##Remove temp file ##
rm -rf $outpath/$sub/visual_$ver/temp_*$TEclass*
