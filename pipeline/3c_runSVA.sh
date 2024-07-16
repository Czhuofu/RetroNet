#!/bin/bash -l

usage="$(basename "$0") [-h] [-o i m r s g f] -- Calling putative MEI insertions. Requires one core, and 20gb memory

where:
    -h  show this help text
    -o  output folder path
    -i  subject ID
    -m  masterpath (default ~/masterpath)
    -r  RetroSom version control (default 1)
    -s  strand (0, -strand; 1, +strand)
    -g  reference genome version (hg38, hg19 etc., default: hg38)
    -f  filter (0, no filter; 1, SEG filter; 2, SEG and size filter"

ver=1
masterpath=~/masterpath
hg=hg38
while getopts ":ho:i:m:r:s:g:i:f:" opt; do
  case $opt in
    h) echo "$usage"
       exit
       ;;
    o) outpath="$OPTARG"
       ;;
    i) sub="$OPTARG"
       ;;
    r) ver="$OPTARG"
       ;;
    m) masterpath="$OPTARG"
       ;;
    s) strand="$OPTARG"
       ;;
    g) hg="$OPTARG"
       ;;
    f) filter="$OPTARG"
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

retro=retro_v$ver\_$strand
rm -r $outpath/$sub/$retro/SVA
mkdir $outpath/$sub/$retro/SVA
tmpfolder=$outpath/$sub

date '+%m/%d/%y %H:%M:%S'
echo "TE calling phase ... begins subject = $sub"

awk '{if ($4 ~ /SVA/) print}' $outpath/$sub/$retro/$sub.sr.dedup.discover | grep 'OK' | \
   awk '{print $1"\t"$2"\t"$3"\t"$5}' | sort -T $tmpfolder -u | \
   sort -T $tmpfolder -k1,1 -k2,3n | \
   mergeBed -d 40 -c 4 -o distinct -delim ";" \
   -i stdin | \
   awk '{split ($4, num, ";"); print $1"\t"$2"\t"$3"\tsr,"length(num)","$4}' \
   > $outpath/$sub/$retro/SVA/$sub.SVA.SR.calls

### PE calling ###
### 1. the candidate reads ###
$masterpath/utls/03_filter_PE_alignment.pl $sub $outpath/$sub $retro 1 $strand

#$masterpath/SEG/seg $outpath/$sub/$retro/seg.pe.fasta 350 -h | \
#    awk '{if (/>/) print $1"\t"$2}' | \
#    sed s/complexity=// > $outpath/$sub/$retro/seg.pe.scores

#$masterpath/utls/04_filter_PE_complexity.pl $sub $outpath/$sub $retro 1

echo "Calling PE phase ... getting the candidate reads"
awk '{if ($4 ~ /SVA/) print}' $outpath/$sub/$retro/$sub.filter.discover | sort -T $tmpfolder -k1 -n -k2 \
    > $outpath/$sub/$retro/SVA/$sub.SVA.novel.sites

### SVA PE calls ###
date '+%m/%d/%y %H:%M:%S'
echo "Calling PE phase ... SVA calls"
awk '{if ($6 == "-") print $1"\t"$2-600"\t"$2+25"\t"NR; else print $1"\t"$3-25"\t"$3+600"\t"NR}' $outpath/$sub/$retro/SVA/$sub.SVA.novel.sites | \
   awk '{if ($2<0) $2=0; print}' OFS='\t' | \
   sort -T $tmpfolder -k1,1 -k2,3n | \
   mergeBed -d 0 -c 4 -o count,distinct -delim $'\t' \
   -i stdin \
   > $outpath/$sub/$retro/SVA/$sub.SVA.PE.calls

$masterpath/utls/07_filter_dup.pl $sub $retro SVA $outpath/$sub 1
$masterpath/utls/08_refine_depth.pl $sub $outpath/$sub $retro SVA

### combine SR and PE ###
date '+%m/%d/%y %H:%M:%S'
echo "Calling PE phase ... Combine SR and PE"
$masterpath/utls/09_merge.SR.PE.support.sh $sub $outpath/$sub $retro SVA $masterpath

### filter out the ref calls and the calls near reference masks ###
### filter out the insertions that are close to segdup, gaps(including telomere) and centromeres ##
date '+%m/%d/%y %H:%M:%S'
echo "Calling PE phase ... ref TE filter"

windowBed \
  -w 100 -v \
  -a $outpath/$sub/$retro/SVA/$sub.SVA.SR.PE.calls \
  -b $masterpath/refTE/position/$hg.fa_SVA\_$strand.bed \
  | windowBed -w 10 -v \
      -a stdin \
      -b $masterpath/refTE/position/$hg.mask.bed \
      > $outpath/$sub/$retro/SVA/$sub.SVA.noref.calls

### filter r1r2 ###
### remove redundancy, keep only one support read when  both read1 and read2 are putative supporting reads ###
### choose SR read over PE read ###
$masterpath/utls/16_prefilter_r1r2.pl $sub $outpath/$sub $retro SVA

