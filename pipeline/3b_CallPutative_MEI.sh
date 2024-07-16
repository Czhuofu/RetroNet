#!/bin/bash
usage="$(basename "$0") [-h] [-o j m v g s f] -- Calling putative MEI insertions. Requires one core, and 20gb memory.

where:
    -h  show this help text
    -o  output folder path
    -j  subject ID
    -m  masterpath
    -v  version control for RetroSom (default 3)
    -g  reference genome (default hg38, supporting hg38, hg19 and b37)
    -s  strand (0, -strand; 1, +strand)
    -f  filter (0, no filter; 1, SEG filter; 2, SEG and size filter"

while getopts ":ho:j:m:v:g:s:f:" opt; do
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
    s) strand="$OPTARG"
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

# Activate RetroSom virtual env #
source /opt/miniconda3/bin/activate RetroSom

retro=retro_v$ver\_$strand
rm -r $outpath/$sub/$retro
mkdir $outpath/$sub/$retro
mkdir $outpath/$sub/$retro/ALU
mkdir $outpath/$sub/$retro/LINE
ln -s $outpath/$sub/retro_v$ver/$sub.discover $outpath/$sub/$retro/$sub.discover
ln -s $outpath/$sub/retro_v$ver/$sub.sr.discover $outpath/$sub/$retro/$sub.sr.discover
ln -s $outpath/$sub/retro_v$ver/$sub.depth $outpath/$sub/$retro/$sub.depth
ln -s $outpath/$sub/retro_v$ver/$sub.bed $outpath/$sub/$retro/$sub.bed
tmpfolder=$outpath/$sub

date '+%m/%d/%y %H:%M:%S'
echo "TE calling phase ... begins"
### get the SEG scores ###
### add cigar to read names ###
awk -v str=$strand '{if (xor(str,($14>$15))) print ">"$5">"$4">"$7">"$12">"$13">"$14">"$15"\n"$16}' $outpath/$sub/$retro/$sub.sr.discover \
   > $outpath/$sub/$retro/seg.sr.fasta

$masterpath/SEG/seg $outpath/$sub/$retro/seg.sr.fasta 350 -h | \
   awk '{if (/>/) print $1"\t"$2}' | \
   sed s/complexity=// > $outpath/$sub/$retro/seg.sr.scores

date '+%m/%d/%y %H:%M:%S'
echo "Modifying SR mode calls"
### filter SR mode ###
$masterpath/utls/05_retroseq_SR_filter_srlength.pl $sub $outpath/$sub $retro $filter
$masterpath/utls/15_filter_SR_dup.pl $sub $outpath/$sub $retro

awk '{if ($4 ~ /Alu/) print}' $outpath/$sub/$retro/$sub.sr.dedup.discover | grep 'OK' | \
   awk '{print $1"\t"$2"\t"$3"\t"$5}' | sort -T $tmpfolder -u | \
   sort -T $tmpfolder -k1,1 -k2,3n | \
   mergeBed -d 40 -c 4 -o distinct -delim ";" \
   -i stdin | \
   awk '{split ($4, num, ";"); print $1"\t"$2"\t"$3"\tsr,"length(num)","$4}' \
   > $outpath/$sub/$retro/ALU/$sub.ALU.SR.calls

awk '{if ($4 ~ /L1/) print}' $outpath/$sub/$retro/$sub.sr.dedup.discover | grep 'OK' | \
   awk '{print $1"\t"$2"\t"$3"\t"$5}' | sort -T $tmpfolder -u | \
   sort -T $tmpfolder -k1,1 -k2,3n | \
   mergeBed -d 40 -c 4 -o distinct -delim ";" \
   -i stdin | \
   awk '{split ($4, num, ";"); print $1"\t"$2"\t"$3"\tsr,"length(num)","$4}' \
   > $outpath/$sub/$retro/LINE/$sub.LINE.SR.calls

### PE calling ###
### 1. the candidate reads ###
date '+%m/%d/%y %H:%M:%S'
echo "Calling PE phase ... filter the candidate support"
$masterpath/utls/03_filter_PE_alignment.pl $sub $outpath/$sub $retro 1 $strand

$masterpath/SEG/seg $outpath/$sub/$retro/seg.pe.fasta 350 -h | \
    awk '{if (/>/) print $1"\t"$2}' | \
    sed s/complexity=// > $outpath/$sub/$retro/seg.pe.scores

$masterpath/utls/04_filter_PE_complexity.pl $sub $outpath/$sub $retro 1

date '+%m/%d/%y %H:%M:%S'
echo "Calling PE phase ... getting the candidate reads"
awk '{if ($4 ~ /Alu/) print}' $outpath/$sub/$retro/$sub.filter.discover | sort -T $tmpfolder -k1 -n -k2 \
    > $outpath/$sub/$retro/ALU/$sub.ALU.novel.sites

awk '{if ($4 ~ /L1/) print}' $outpath/$sub/$retro/$sub.filter.discover | sort -T $tmpfolder -k1 -n -k2 \
    > $outpath/$sub/$retro/LINE/$sub.LINE.novel.sites

### ALU PE calls ###
date '+%m/%d/%y %H:%M:%S'
echo "Calling PE phase ... Alu calls"
awk '{if ($6 == "-") print $1"\t"$2-600"\t"$2+25"\t"NR; else print $1"\t"$3-25"\t"$3+600"\t"NR}' $outpath/$sub/$retro/ALU/$sub.ALU.novel.sites | \
   awk '{if ($2<0) $2=0; print}' OFS='\t' | \
   sort -T $tmpfolder -k1,1 -k2,3n | \
   mergeBed -d 0 -c 4 -o count,distinct -delim $'\t' \
   -i stdin \
   > $outpath/$sub/$retro/ALU/$sub.ALU.PE.calls

$masterpath/utls/07_filter_dup.pl $sub $retro ALU $outpath/$sub 1
$masterpath/utls/08_refine_depth.pl $sub $outpath/$sub $retro ALU

### L1 PE calls ###
date '+%m/%d/%y %H:%M:%S'
echo "Calling PE phase ... L1 calls"
awk '{if ($6 == "-") print $1"\t"$2-600"\t"$2+25"\t"NR; else print $1"\t"$3-25"\t"$3+600"\t"NR}' $outpath/$sub/$retro/LINE/$sub.LINE.novel.sites | \
   awk '{if ($2<0) $2=0; print}' OFS='\t' | \
   sort -T $tmpfolder -k1,1 -k2,3n | \
   mergeBed -d 0 -c 4 -o count,distinct -delim $'\t' \
   -i stdin \
   > $outpath/$sub/$retro/LINE/$sub.LINE.PE.calls

$masterpath/utls/07_filter_dup.pl $sub $retro LINE $outpath/$sub 1
$masterpath/utls/08_refine_depth.pl $sub $outpath/$sub $retro LINE

### combine SR and PE ###
#### need a bash for ubuntu
date '+%m/%d/%y %H:%M:%S'
echo "Calling PE phase ... Combine SR and PE"
bash $masterpath/utls/09_merge.SR.PE.support.sh $sub $outpath/$sub $retro ALU $masterpath
bash $masterpath/utls/09_merge.SR.PE.support.sh $sub $outpath/$sub $retro LINE $masterpath

### filter out the ref calls and the calls near reference masks ###
### filter out the insertions that are close to segdup, gaps(including telomere) and centromeres ##
date '+%m/%d/%y %H:%M:%S'
echo "Calling PE phase ... ref TE filter"

windowBed \
  -w 100 -v \
  -a $outpath/$sub/$retro/ALU/$sub.ALU.SR.PE.calls \
  -b $masterpath/refTE/position/$hg.fa_ALU\_$strand.bed \
  | windowBed -w 10 -v \
      -a stdin \
      -b $masterpath/refTE/position/$hg.mask.bed \
      > $outpath/$sub/$retro/ALU/$sub.ALU.noref.calls

windowBed \
  -w 100 -v \
  -a $outpath/$sub/$retro/LINE/$sub.LINE.SR.PE.calls \
  -b $masterpath/refTE/position/$hg.fa_LINE1\_$strand.bed \
  | windowBed -w 10 -v \
      -a stdin \
      -b $masterpath/refTE/position/$hg.mask.bed \
      > $outpath/$sub/$retro/LINE/$sub.LINE.noref.calls

### filter r1r2 ###
### remove redundancy, keep only one support read when  both read1 and read2 are putative supporting reads ###
### choose SR read over PE read ###
$masterpath/utls/16_prefilter_r1r2.pl $sub $outpath/$sub $retro ALU
$masterpath/utls/16_prefilter_r1r2.pl $sub $outpath/$sub $retro LINE
