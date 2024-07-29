#!/bin/bash
usage="$(basename "$0") [-h] [-i t f c d e r s p m g] -- Visualization tool for MEI supporting reads, e.g.,
./RetroVis.sh -i neg_NoModel -t LINE -f L1HS -c chr3 -d 153168189 -e 153168814 -r 61 -s 0 -p /home/xwzhu/transfer/BulkSeq/12004/ -m /home/xwzhu/masterpath -g hg38

where:
    -h  show this help text
    -i  subject ID
    -t  TE class: LINE/ALU
    -f  TE family: L1HS/AluYa5
    -c  MEI chromosome
    -d  MEI coordinate 1
    -e  MEI coordinate 2
    -r  version control for RetroSom (default 1)
    -s  strandness of the MEI (+ strand: 1, - strand: 0)
    -p  datapath to the subject
    -m  masterpath (default ~/masterpath)
    -g  huamn reference genome (hg19/hg38)
"
while getopts ":hi:t:f:c:d:e:r:s:p:m:g:l:a:" opt; do
  case $opt in
    h) echo "$usage"
       exit
       ;;
    i) subject="$OPTARG"
       ;;
    t) TEclass="$OPTARG"
       ;;
    f) TEfamily="$OPTARG"
       ;;
    c) chr="$OPTARG"
       ;;
    d) cord1="$OPTARG"
       ;;
    e) cord2="$OPTARG"
       ;;
    r) ver="$OPTARG"
       ;;
    s) strand="$OPTARG"
       ;;
    p) datapath="$OPTARG"
       ;;
    m) masterpath="$OPTARG"
       ;;
    g) hg="$OPTARG"
       ;;
    l) label="$OPTARG"
       ;;
    a) filename="$OPTARG"
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

# mkdir $datapath/$subject/visual
rm -r $datapath/$subject/visual_${ver}/temp_$filename
mkdir $datapath/$subject/visual_${ver}/temp_$filename
echo -e $chr"\t"$cord1"\t"$cord2 > $datapath/$subject/visual_${ver}/temp_$filename/00_cord.bed
retro=retro_v$ver\_$strand
ref=$masterpath/refTE/sequence/$TEfamily.fa
callfile=$datapath/$subject/retro_v$ver\_$strand/$TEclass/$subject.$TEclass.SR.PE.calls
PEfile=$datapath/$subject/retro_v$ver\_$strand/$TEclass/$subject.$TEclass.novel.sites
SRfile=$datapath/$subject/retro_v$ver\_$strand/$subject.sr.tabe.discover

### Step1: extract the MEI overlapping with the query insertion ###
awk '{printf $1"\t"$2"\t"$3"\t"; for (i=6; i<=NF; ++i) {printf "%s,", $i; if ($i == $NF) printf "\n"}}' \
  $callfile | \
  intersectBed -wo \
   -a stdin \
   -b $datapath/$subject/visual_${ver}/temp_$filename/00_cord.bed \
   > $datapath/$subject/visual_${ver}/temp_$filename/01_overlapping_insertions.txt

### Step2: extract the supporting reads ###
./MEI_support_reads.pl $PEfile $SRfile $masterpath $TEclass $datapath/$subject/visual_${ver}/temp_$filename $datapath/$subject/$retro/$TEclass/split-reads-anchor.txt

### Step3: Realign the reads to TEfamily consensus sequence ###
petemp=$datapath/$subject/visual_${ver}/temp_$filename/02_PEsupport.fa
if [ -s "$petemp" ] && [ $TEclass == 'LINE' ]
then
env TMPDIR=$datapath/$subject exonerate --model affine:local \
   --bestn 1 --ryo "INFO: %qi %qal %pi %tS %ti %qab %qae %tab %tae %ql %tas\n" \
   $petemp \
   $ref \
   > $datapath/$subject/visual_${ver}/temp_$filename/04_PE.alignment
fi

if [ -s "$petemp" ] && [ $TEclass == 'ALU' ]
then
ref=$masterpath/refTE/sequence/AluYa5.fa
env TMPDIR=$datapath/$subject exonerate --model affine:local \
   --bestn 1 --ryo "INFO: %qi %qal %pi %tS %ti %qab %qae %tab %tae %ql %tas\n" \
   $petemp \
   $ref \
   > $datapath/$subject/visual_${ver}/temp_$filename/04Ya_PE.alignment

ref=$masterpath/refTE/sequence/AluYb8.fa
env TMPDIR=$datapath/$subject exonerate --model affine:local \
   --bestn 1 --ryo "INFO: %qi %qal %pi %tS %ti %qab %qae %tab %tae %ql %tas\n" \
   $petemp \
   $ref \
   > $datapath/$subject/visual_${ver}/temp_$filename/04Yb_PE.alignment

ref=$masterpath/refTE/sequence/AluYc1.fa
env TMPDIR=$datapath/$subject exonerate --model affine:local \
   --bestn 1 --ryo "INFO: %qi %qal %pi %tS %ti %qab %qae %tab %tae %ql %tas\n" \
   $petemp \
   $ref \
   > $datapath/$subject/visual_${ver}/temp_$filename/04Yc_PE.alignment

ref=$masterpath/refTE/sequence/AluYk13.fa
env TMPDIR=$datapath/$subject exonerate --model affine:local \
   --bestn 1 --ryo "INFO: %qi %qal %pi %tS %ti %qab %qae %tab %tae %ql %tas\n" \
   $petemp \
   $ref \
   > $datapath/$subject/visual_${ver}/temp_$filename/04Yk_PE.alignment
fi

if [ -s "$petemp" ] && [ $TEclass == 'SVA' ]
then
ref=$masterpath/refTE/sequence/SVA_E.fa
env TMPDIR=$datapath/$subject exonerate --model affine:local \
   --bestn 1 --ryo "INFO: %qi %qal %pi %tS %ti %qab %qae %tab %tae %ql %tas\n" \
   $petemp \
   $ref \
   > $datapath/$subject/visual_${ver}/temp_$filename/04SE_PE.alignment

ref=$masterpath/refTE/sequence/SVA_F.fa
env TMPDIR=$datapath/$subject exonerate --model affine:local \
   --bestn 1 --ryo "INFO: %qi %qal %pi %tS %ti %qab %qae %tab %tae %ql %tas\n" \
   $petemp \
   $ref \
   > $datapath/$subject/visual_${ver}/temp_$filename/04SF_PE.alignment
fi

srtemp=$datapath/$subject/visual_${ver}/temp_$filename/03_SRsupport.fa
if [ -s "$srtemp" ] && [ $TEclass  == 'LINE' ]
then
env TMPDIR=$datapath/$subject exonerate --model affine:local \
   --bestn 1 --ryo "INFO: %qi %qal %pi %tS %ti %qab %qae %tab %tae %tas\n" \
   $srtemp \
   $ref \
   > $datapath/$subject/visual_${ver}/temp_$filename/05_SR.alignment
fi

if [ -s "$srtemp" ] && [ $TEclass == 'ALU' ]
then
ref=$masterpath/refTE/sequence/AluYa5.fa
env TMPDIR=$datapath/$subject exonerate --model affine:local \
   --bestn 1 --ryo "INFO: %qi %qal %pi %tS %ti %qab %qae %tab %tae %tas\n" \
   $srtemp \
   $ref \
   > $datapath/$subject/visual_${ver}/temp_$filename/05Ya_SR.alignment

ref=$masterpath/refTE/sequence/AluYb8.fa
env TMPDIR=$datapath/$subject exonerate --model affine:local \
   --bestn 1 --ryo "INFO: %qi %qal %pi %tS %ti %qab %qae %tab %tae %tas\n" \
   $srtemp \
   $ref \
   > $datapath/$subject/visual_${ver}/temp_$filename/05Yb_SR.alignment

ref=$masterpath/refTE/sequence/AluYc1.fa
env TMPDIR=$datapath/$subject exonerate --model affine:local \
   --bestn 1 --ryo "INFO: %qi %qal %pi %tS %ti %qab %qae %tab %tae %tas\n" \
   $srtemp \
   $ref \
   > $datapath/$subject/visual_${ver}/temp_$filename/05Yc_SR.alignment

ref=$masterpath/refTE/sequence/AluYk13.fa
env TMPDIR=$datapath/$subject exonerate --model affine:local \
   --bestn 1 --ryo "INFO: %qi %qal %pi %tS %ti %qab %qae %tab %tae %tas\n" \
   $srtemp \
   $ref \
   > $datapath/$subject/visual_${ver}/temp_$filename/05Yk_SR.alignment
fi

if [ -s "$srtemp" ] && [ $TEclass == 'SVA' ]
then
ref=$masterpath/refTE/sequence/SVA_E.fa
env TMPDIR=$datapath/$subject exonerate --model affine:local \
   --bestn 1 --ryo "INFO: %qi %qal %pi %tS %ti %qab %qae %tab %tae %tas\n" \
   $srtemp \
   $ref \
   > $datapath/$subject/visual_${ver}/temp_$filename/05SE_SR.alignment

ref=$masterpath/refTE/sequence/SVA_F.fa
env TMPDIR=$datapath/$subject exonerate --model affine:local \
   --bestn 1 --ryo "INFO: %qi %qal %pi %tS %ti %qab %qae %tab %tae %tas\n" \
   $srtemp \
   $ref \
   > $datapath/$subject/visual_${ver}/temp_$filename/05SF_SR.alignment
fi

### step4: plotting the supporting reads ### 
imgfile=$subject\_retro$ver\_strand$strand\_$TEfamily\_$chr\_$cord1\_$label\.png
matfile=$datapath/$subject/retro_v$ver\_$strand/$subject.pe.$TEclass.matrix
if [ $hg == 'hg38']
then
mapfile=$masterpath/RetroNet/hg38_100bp.bedGraph
elif [ $hg == 'b37' ]
then
mapfile=$masterpath/RetroNet/b37_100bp.bedGraph
fi

if [ $TEclass == 'LINE' ]
then
    ./LINEvis.pl $masterpath $imgfile $chr $TEfamily $datapath/$subject/visual_${ver}/${TEclass} $hg $strand $matfile $cord1 $datapath/$subject/visual_${ver}/temp_$filename $mapfile
elif [ $TEclass == 'ALU' ]
then
    ./ALUvis.pl $masterpath $imgfile $chr $TEfamily $datapath/$subject/visual_${ver}/${TEclass} $hg $strand $matfile $cord1 $datapath/$subject/visual_${ver}/temp_$filename $mapfile
elif [ $TEclass == 'SVA' ]
then
    ./SVAvis.pl $masterpath $imgfile $chr $TEfamily $datapath/$subject/visual_${ver}/${TEclass} $hg $strand $matfile $cord1 $datapath/$subject/visual_${ver}/temp_$filename $mapfile
fi
