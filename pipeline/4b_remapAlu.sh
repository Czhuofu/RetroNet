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
echo "$sub ran on $(hostname -s)"

retro=retro_v$ver
Retro=fix_GGCG
rm -r $outpath/$sub/$Retro
mkdir $outpath/$sub/$Retro
ln -s $outpath/$sub/$retro/$sub.discover $outpath/$sub/$Retro/$sub.discover
ln -s $outpath/$sub/$retro/$sub.sr.discover $outpath/$sub/$Retro/$sub.sr.discover
 $masterpath/ALU/00_remap/01_GGCG.sh $sub $outpath/$sub $Retro $masterpath

Retro=fix_GAGC
rm -r $outpath/$sub/$Retro
mkdir $outpath/$sub/$Retro
ln -s $outpath/$sub/$retro/$sub.discover $outpath/$sub/$Retro/$sub.discover
ln -s $outpath/$sub/$retro/$sub.sr.discover $outpath/$sub/$Retro/$sub.sr.discover
 $masterpath/ALU/00_remap/02_GAGC.sh $sub $outpath/$sub $Retro $masterpath

Retro=fix_CC
rm -r $outpath/$sub/$Retro
mkdir $outpath/$sub/$Retro
ln -s $outpath/$sub/$retro/$sub.discover $outpath/$sub/$Retro/$sub.discover
ln -s $outpath/$sub/$retro/$sub.sr.discover $outpath/$sub/$Retro/$sub.sr.discover
 $masterpath/ALU/00_remap/03_CC.sh $sub $outpath/$sub $Retro $masterpath

Retro=fix_CG
rm -r $outpath/$sub/$Retro
mkdir $outpath/$sub/$Retro
ln -s $outpath/$sub/$retro/$sub.discover $outpath/$sub/$Retro/$sub.discover
ln -s $outpath/$sub/$retro/$sub.sr.discover $outpath/$sub/$Retro/$sub.sr.discover
 $masterpath/ALU/00_remap/03_CG.sh $sub $outpath/$sub $Retro $masterpath

Retro=fix_TG
rm -r $outpath/$sub/$Retro
mkdir $outpath/$sub/$Retro
ln -s $outpath/$sub/$retro/$sub.discover $outpath/$sub/$Retro/$sub.discover
ln -s $outpath/$sub/$retro/$sub.sr.discover $outpath/$sub/$Retro/$sub.sr.discover
 $masterpath/ALU/00_remap/03_TG.sh $sub $outpath/$sub $Retro $masterpath

Retro=fix_AluYb
rm -r $outpath/$sub/$Retro
mkdir $outpath/$sub/$Retro
ln -s $outpath/$sub/$retro/$sub.discover $outpath/$sub/$Retro/$sub.discover
ln -s $outpath/$sub/$retro/$sub.sr.discover $outpath/$sub/$Retro/$sub.sr.discover
 $masterpath/ALU/00_remap/04_AluYb.sh $sub $outpath/$sub $Retro $masterpath

Retro=fix_AT
rm -r $outpath/$sub/$Retro
mkdir $outpath/$sub/$Retro
ln -s $outpath/$sub/$retro/$sub.discover $outpath/$sub/$Retro/$sub.discover
ln -s $outpath/$sub/$retro/$sub.sr.discover $outpath/$sub/$Retro/$sub.sr.discover
 $masterpath/ALU/00_remap/05_AT.sh $sub $outpath/$sub $Retro $masterpath