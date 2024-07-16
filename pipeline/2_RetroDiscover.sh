#!/bin/bash

usage="$(basename "$0") [-h] [-o j m v] -- Discover candidate supporting reads with modified Retroseq algorithm.

where:
    -h  show this help text
    -o  output folder path
    -j  subject ID
    -m  masterpath
    -v  version control for RetroSom (default 3)"

while getopts ":ho:j:m:v:c:" opt; do
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

Retro=retro_v$ver

# Activate RetroSom virtual env #
source /opt/miniconda3/bin/activate RetroSom
export TMPDIR=$outpath/$sub/temp
date '+%m/%d/%y %H:%M:%S'
echo "Discovery phase ... begins"
$masterpath/RetroSeq/retroseq.prun.pl -discover \
  -align \
  -srmode \
  -minclip 20 \
  -len 26 \
  -srcands $outpath/$sub/$Retro/$sub.sr.discover \
  -bam $outpath/$sub/align/$sub.final.bam \
  -eref $masterpath/refTE/TE_ALHS.bed \
  -output $outpath/$sub/$Retro/$sub.discover

### seperate the supporting reads from either strand ###
$masterpath/utls/21_par_direction.pl $sub $outpath $Retro
date '+%m/%d/%y %H:%M:%S'
echo "Discovery phase ... finish"
