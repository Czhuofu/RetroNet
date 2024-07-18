#!/bin/bash

usage="$(basename "$0") [-h] [-o j m v g s l] -- Pairing the supporting reads for the same MEI

where:
    -h  show this help text
    -o  output folder path
    -j  subject ID
    -m  masterpath
    -v  version control for RetroSom (default 3)
    -g  reference genome (default hg38, supporting hg38, hg19 and b37)
    -s  strandness (0, -strand; 1, +strand)
    -l  number of the libraries"

while getopts ":ho:j:m:v:c:g:s:l:" opt; do
  case $opt in
    h) echo "$usage"
       exit
       ;;
    o) outpath=$(readlink -f $OPTARG)
       ;;
    j) subject="$OPTARG"
       ;;
    m) masterpath=$(readlink -f $OPTARG)
       ;;
    v) ver="$OPTARG"
       ;;
    g) hg="$OPTARG"
       ;;
    s) strand="$OPTARG"
       ;;
    l) readgroup="$OPTARG"
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

#echo $1  ### subject ID ###
#echo $2  ### path ###
#echo retro_v$ver\_$strand ### output folder ###
#echo $4  ### reference genome hg19 G37 hg38... ###
#echo $5 ##strandness ###
sub=$subject\_Combined
sub2=$subject\-
startlib=1

TE=LINE
tmpfolder=$outpath/$sub

mkdir $outpath/$sub
mkdir $outpath/$sub/retro_v$ver\_$strand
rm -r $outpath/$sub/retro_v$ver\_$strand/$TE
mkdir $outpath/$sub/retro_v$ver\_$strand/$TE

### combine 6 invidual libraries ###
### combine SR support reads ###
rm $outpath/$sub/retro_v$ver\_$strand/$sub.sr.tabe.discover
for (( c=startlib; c<=$readgroup; c++ ))
     do
       awk -v lib=$c '$5="lib"lib"_"$5' OFS='\t' $outpath/$sub2$c/retro_v$ver\_$strand/$sub2$c.sr.tabe.discover >> $outpath/$sub/retro_v$ver\_$strand/$sub.sr.tabe.discover
     done

### combine the PE reads ###
rm $outpath/$sub/retro_v$ver\_$strand/$TE/$sub.$TE.novel.sites
for (( c=startlib; c<=$readgroup; c++ ))
     do
       awk -v lib=$c '$5="lib"lib"_"$5' OFS='\t' $outpath/$sub2$c/retro_v$ver\_$strand/$TE/$sub2$c.$TE.nodup.sites >> $outpath/$sub/retro_v$ver\_$strand/$TE/$sub.$TE.novel.sites
     done

### combine the anchor ###
echo -e "read\tdirect\tinsert\tseg\talign" > $outpath/$sub/retro_v$ver\_$strand/$TE/split-reads-anchor.txt
for (( c=startlib; c<=$readgroup; c++ ))
     do
       awk '{if (!/^read/) print $0}' $outpath/$sub2$c/retro_v$ver\_$strand/$TE/split-reads-anchor.txt | awk -v lib=$c '{split ($1, names, ":"); print names[1]":lib"lib"_"names[2]":"names[3]":"names[4]":"names[5]":"names[6]":"names[7]":"names[8]"\t"$2"\t"$3"\t"$4"\t"$5}' >> $outpath/$sub/retro_v$ver\_$strand/$TE/split-reads-anchor.txt
     done

### make calls ###
awk '{if ($4 ~ /L1/) print}' $outpath/$sub/retro_v$ver\_$strand/$sub.sr.tabe.discover | grep 'OK' | \
   awk '{print $1"\t"$2"\t"$3"\t"$5}' | sort -T $tmpfolder -u | \
   sort -T $tmpfolder -k1,1 -k2,3n | \
   mergeBed -d 40 -c 4 -o distinct -delim ";" \
   -i stdin | \
   awk '{split ($4, num, ";"); print $1"\t"$2"\t"$3"\tsr,"length(num)","$4}' \
   > $outpath/$sub/retro_v$ver\_$strand/$TE/$sub.$TE.SR.calls

### PE calling ###
### L1 PE calls ###
$masterpath/utls/06_ana_depth.pl $sub $outpath/$sub retro_v$ver\_$strand $TE 1000
mv $outpath/$sub/retro_v$ver\_$strand/$TE/$sub.$TE.PE.calls $outpath/$sub/retro_v$ver\_$strand/$TE/$sub.$TE.PE.nodup.calls
$masterpath/utls/08_refine_depth.pl $sub $outpath/$sub retro_v$ver\_$strand $TE

### combine SR and PE ###
date '+%m/%d/%y %H:%M:%S'
echo "Calling PE phase ... Combine SR and PE"
$masterpath/utls/09_merge.SR.PE.support.sh $sub $outpath/$sub retro_v$ver\_$strand $TE $masterpath

windowBed \
  -w 100 -v \
  -a $outpath/$sub/retro_v$ver\_$strand/LINE/$sub.LINE.SR.PE.calls \
  -b $masterpath/refTE/position/$hg\.fa_LINE1\_$strand.bed \
  | windowBed -w 10 -v \
      -a stdin \
      -b $masterpath/refTE/position/$hg.mask.bed \
      > $outpath/$sub/retro_v$ver\_$strand/LINE/$sub.LINE.noref.calls

$masterpath/utls/16_prefilter_r1r2.pl $sub $outpath/$sub retro_v$ver\_$strand LINE


for (( c=startlib; c<=$readgroup; c++ ))
     do
       awk -v lib=$c '$4="lib"lib"_"$4' OFS='\t' $outpath/$sub2$c/retro_v$ver\_$strand/$sub2$c.pe.LINE.matrix | awk 'FNR > 1' >> $outpath/$sub/retro_v$ver\_$strand/$sub.pe.LINE.matrix
       awk -v lib=$c '$4="lib"lib"_"$4' OFS='\t' $outpath/$sub2$c/retro_v$ver\_$strand/$sub2$c.sr.LINE.matrix | awk 'FNR > 1' >> $outpath/$sub/retro_v$ver\_$strand/$sub.sr.LINE.matrix
     done

TE=ALU
mkdir $outpath/$sub/retro_v$ver\_$strand/$TE

### combine the PE reads ###
rm $outpath/$sub/retro_v$ver\_$strand/$TE/$sub.$TE.novel.sites
for (( c=startlib; c<=$readgroup; c++ ))
     do
       awk -v lib=$c '$5="lib"lib"_"$5' OFS='\t' $outpath/$sub2$c/retro_v$ver\_$strand/$TE/$sub2$c.$TE.nodup.sites >> $outpath/$sub/retro_v$ver\_$strand/$TE/$sub.$TE.novel.sites
     done

### combine the anchor ###
echo -e "read\tdirect\tinsert\tseg\talign" > $outpath/$sub/retro_v$ver\_$strand/$TE/split-reads-anchor.txt
for (( c=startlib; c<=$readgroup; c++ ))
     do
       awk '{if (!/^read/) print $0}' $outpath/$sub2$c/retro_v$ver\_$strand/$TE/split-reads-anchor.txt | awk -v lib=$c '{split ($1, names, ":"); print names[1]":lib"lib"_"names[2]":"names[3]":"names[4]":"names[5]":"names[6]":"names[7]":"names[8]"\t"$2"\t"$3"\t"$4"\t"$5}' >> $outpath/$sub/retro_v$ver\_$strand/$TE/split-reads-anchor.txt
     done

### make calls ###
awk '{if ($4 ~ /Alu/) print}' $outpath/$sub/retro_v$ver\_$strand/$sub.sr.tabe.discover | grep 'OK' | \
   awk '{print $1"\t"$2"\t"$3"\t"$5}' | sort -T $tmpfolder -u | \
   sort -T $tmpfolder -k1,1 -k2,3n | \
   mergeBed -d 40 -c 4 -o distinct -delim ";" \
   -i stdin | \
   awk '{split ($4, num, ";"); print $1"\t"$2"\t"$3"\tsr,"length(num)","$4}' \
   > $outpath/$sub/retro_v$ver\_$strand/$TE/$sub.$TE.SR.calls

### PE calling ###
### L1 PE calls ###
$masterpath/utls/06_ana_depth.pl $sub $outpath/$sub retro_v$ver\_$strand $TE 1000
mv $outpath/$sub/retro_v$ver\_$strand/$TE/$sub.$TE.PE.calls $outpath/$sub/retro_v$ver\_$strand/$TE/$sub.$TE.PE.nodup.calls
$masterpath/utls/08_refine_depth.pl $sub $outpath/$sub retro_v$ver\_$strand $TE

### combine SR and PE ###
date '+%m/%d/%y %H:%M:%S'
echo "Calling PE phase ... Combine SR and PE"
$masterpath/utls/09_merge.SR.PE.support.sh $sub $outpath/$sub retro_v$ver\_$strand $TE $masterpath

windowBed \
  -w 100 -v \
  -a $outpath/$sub/retro_v$ver\_$strand/ALU/$sub.ALU.SR.PE.calls \
  -b $masterpath/refTE/position/$hg\.fa_ALU\_$strand.bed \
  | windowBed -w 10 -v \
      -a stdin \
      -b $masterpath/refTE/position/$hg.mask.bed \
      > $outpath/$sub/retro_v$ver\_$strand/ALU/$sub.ALU.noref.calls

$masterpath/utls/16_prefilter_r1r2.pl $sub $outpath/$sub retro_v$ver\_$strand ALU

head -n 1 $outpath/$sub2\1/retro_v$ver\_$strand/$sub2\1.pe.ALU.matrix > $outpath/$sub/retro_v$ver\_$strand/$sub.pe.ALU.matrix
head -n 1 $outpath/$sub2\1/retro_v$ver\_$strand/$sub2\1.sr.ALU.matrix > $outpath/$sub/retro_v$ver\_$strand/$sub.sr.ALU.matrix

for (( c=startlib; c<=$readgroup; c++ ))
     do
       awk -v lib=$c '$4="lib"lib"_"$4' OFS='\t' $outpath/$sub2$c/retro_v$ver\_$strand/$sub2$c.pe.ALU.matrix | awk 'FNR > 1' >> $outpath/$sub/retro_v$ver\_$strand/$sub.pe.ALU.matrix
       awk -v lib=$c '$4="lib"lib"_"$4' OFS='\t' $outpath/$sub2$c/retro_v$ver\_$strand/$sub2$c.sr.ALU.matrix | awk 'FNR > 1' >> $outpath/$sub/retro_v$ver\_$strand/$sub.sr.ALU.matrix
     done

###sva combine
TE=SVA
mkdir $outpath/$sub/retro_v$ver\_$strand/$TE

### combine the PE reads ###
rm $outpath/$sub/retro_v$ver\_$strand/$TE/$sub.$TE.novel.sites
for (( c=startlib; c<=$readgroup; c++ ))
     do
       awk -v lib=$c '$5="lib"lib"_"$5' OFS='\t' $outpath/$sub2$c/retro_v$ver\_$strand/$TE/$sub2$c.$TE.nodup.sites >> $outpath/$sub/retro_v$ver\_$strand/$TE/$sub.$TE.novel.sites
     done

### combine the anchor ###
echo -e "read\tdirect\tinsert\tseg\talign" > $outpath/$sub/retro_v$ver\_$strand/$TE/split-reads-anchor.txt
for (( c=startlib; c<=$readgroup; c++ ))
     do
       awk '{if (!/^read/) print $0}' $outpath/$sub2$c/retro_v$ver\_$strand/$TE/split-reads-anchor.txt | awk -v lib=$c '{split ($1, names, ":"); print names[1]":lib"lib"_"names[2]":"names[3]":"names[4]":"names[5]":"names[6]":"names[7]":"names[8]"\t"$2"\t"$3"\t"$4"\t"$5}' >> $outpath/$sub/retro_v$ver\_$strand/$TE/split-reads-anchor.txt
     done

### make calls ###
awk '{if ($4 ~ /Alu/) print}' $outpath/$sub/retro_v$ver\_$strand/$sub.sr.tabe.discover | grep 'OK' | \
   awk '{print $1"\t"$2"\t"$3"\t"$5}' | sort -T $tmpfolder -u | \
   sort -T $tmpfolder -k1,1 -k2,3n | \
   mergeBed -d 40 -c 4 -o distinct -delim ";" \
   -i stdin | \
   awk '{split ($4, num, ";"); print $1"\t"$2"\t"$3"\tsr,"length(num)","$4}' \
   > $outpath/$sub/retro_v$ver\_$strand/$TE/$sub.$TE.SR.calls

### PE calling ###
### L1 PE calls ###
$masterpath/utls/06_ana_depth.pl $sub $outpath/$sub retro_v$ver\_$strand $TE 1000
mv $outpath/$sub/retro_v$ver\_$strand/$TE/$sub.$TE.PE.calls $outpath/$sub/retro_v$ver\_$strand/$TE/$sub.$TE.PE.nodup.calls
$masterpath/utls/08_refine_depth.pl $sub $outpath/$sub retro_v$ver\_$strand $TE

### combine SR and PE ###
date '+%m/%d/%y %H:%M:%S'
echo "Calling PE phase ... Combine SR and PE"
$masterpath/utls/09_merge.SR.PE.support.sh $sub $outpath/$sub retro_v$ver\_$strand $TE $masterpath

windowBed \
  -w 100 -v \
  -a $outpath/$sub/retro_v$ver\_$strand/SVA/$sub.SVA.SR.PE.calls \
  -b $masterpath/refTE/position/$hg\.fa_SVA\_$strand.bed \
  | windowBed -w 10 -v \
      -a stdin \
      -b $masterpath/refTE/position/$hg.mask.bed \
      > $outpath/$sub/retro_v$ver\_$strand/SVA/$sub.SVA.noref.calls

$masterpath/utls/16_prefilter_r1r2.pl $sub $outpath/$sub retro_v$ver\_$strand SVA
