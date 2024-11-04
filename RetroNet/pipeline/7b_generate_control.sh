#!/bin/bash

usage="$(basename "$0") [-h] [-o j m v l] -- Pairing the supporting reads for the same MEI

where:
    -h  show this help text
    -o  output folder path
    -j  subject ID
    -m  masterpath
    -v  version control for RetroSom (default 3)
    -l  number of the libraries"

while getopts ":ho:j:m:v:c:l:" opt; do
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

sub=$subject\_Combined

if [ "$readgroup" != 1 ]
then
    sub2=$subject\-
    startlib=1
    for (( c=startlib; c<=$readgroup; c++ ))
       do       
         grep "L1HS" $outpath/$sub2$c/retro_v$ver/$sub2$c.0.discover | awk '{if ($8 >= 90) print}' | awk '{if ($11 <= 6000 || $12 <=6000) print $1"\t"$2"\t"$3"\t"$5}' | sort -u >> /$outpath/$sub/retro_v$ver\_0/LINE/$sub.LINE.nofilter.0.bed
         grep "L1HS" $outpath/$sub2$c/retro_v$ver/$sub2$c.1.discover | awk '{if ($8 >= 90) print}' | awk '{if ($11 <= 6000 || $12 <=6000) print $1"\t"$2"\t"$3"\t"$5}' | sort -u >> /$outpath/$sub/retro_v$ver\_1/LINE/$sub.LINE.nofilter.1.bed
         grep "Alu" $outpath/$sub2$c/retro_v$ver/$sub2$c.0.discover | awk '{if ($8 >= 90) print $1"\t"$2"\t"$3"\t"$5}' | sort -u >> /$outpath/$sub/retro_v$ver\_0/ALU/$sub.ALU.nofilter.0.bed
         grep "Alu" $outpath/$sub2$c/retro_v$ver/$sub2$c.1.discover | awk '{if ($8 >= 90) print $1"\t"$2"\t"$3"\t"$5}' | sort -u >> /$outpath/$sub/retro_v$ver\_1/ALU/$sub.ALU.nofilter.1.bed
         grep "SVA" $outpath/$sub2$c/retro_v$ver/$sub2$c.0.discover | awk '{if ($8 >= 90) print $1"\t"$2"\t"$3"\t"$5}' | sort -u >> /$outpath/$sub/retro_v$ver\_0/SVA/$sub.SVA.nofilter.0.bed
         grep "SVA" $outpath/$sub2$c/retro_v$ver/$sub2$c.1.discover | awk '{if ($8 >= 90) print $1"\t"$2"\t"$3"\t"$5}' | sort -u >> /$outpath/$sub/retro_v$ver\_1/SVA/$sub.SVA.nofilter.1.bed
       done
    awk '{if (length($1) < 6) print}' /$outpath/$sub/retro_v$ver\_0/LINE/$sub.LINE.nofilter.0.bed | sort -k1,1V -k2,2n -u -o /$outpath/$sub/retro_v$ver\_0/LINE/$sub.LINE.nofilter.temp.bed
    awk '{if (length($1) < 6) print}' /$outpath/$sub/retro_v$ver\_1/LINE/$sub.LINE.nofilter.1.bed | sort -k1,1V -k2,2n -u -o /$outpath/$sub/retro_v$ver\_1/LINE/$sub.LINE.nofilter.temp.bed
    awk '{if (length($1) < 6) print}' /$outpath/$sub/retro_v$ver\_0/ALU/$sub.ALU.nofilter.0.bed | sort -k1,1V -k2,2n -u -o /$outpath/$sub/retro_v$ver\_0/ALU/$sub.ALU.nofilter.temp.bed
    awk '{if (length($1) < 6) print}' /$outpath/$sub/retro_v$ver\_1/ALU/$sub.ALU.nofilter.1.bed | sort -k1,1V -k2,2n -u -o /$outpath/$sub/retro_v$ver\_1/ALU/$sub.ALU.nofilter.temp.bed
    awk '{if (length($1) < 6) print}' /$outpath/$sub/retro_v$ver\_0/SVA/$sub.SVA.nofilter.0.bed | sort -k1,1V -k2,2n -u -o /$outpath/$sub/retro_v$ver\_0/SVA/$sub.SVA.nofilter.temp.bed
    awk '{if (length($1) < 6) print}' /$outpath/$sub/retro_v$ver\_1/SVA/$sub.SVA.nofilter.1.bed | sort -k1,1V -k2,2n -u -o /$outpath/$sub/retro_v$ver\_1/SVA/$sub.SVA.nofilter.temp.bed
    mv $outpath/$sub/retro_v$ver\_0/LINE/$sub.LINE.SR.PE.calls $outpath/$sub/retro_v$ver\_0/LINE/$sub.LINE.SR.PE.old.calls
    awk '{if ($5 > 1) print}' $outpath/$sub/retro_v$ver\_0/LINE/$sub.LINE.SR.PE.old.calls > $outpath/$sub/retro_v$ver\_0/LINE/$sub.LINE.SR.PE.calls
    bedtools merge -i /$outpath/$sub/retro_v$ver\_0/LINE/$sub.LINE.nofilter.temp.bed -c 1 -o count | awk '{if ($4 >1) print}' >> $outpath/$sub/retro_v$ver\_0/LINE/$sub.LINE.SR.PE.calls
    mv $outpath/$sub/retro_v$ver\_1/LINE/$sub.LINE.SR.PE.calls $outpath/$sub/retro_v$ver\_1/LINE/$sub.LINE.SR.PE.old.calls
    awk '{if ($5 > 1) print}' $outpath/$sub/retro_v$ver\_1/LINE/$sub.LINE.SR.PE.old.calls > $outpath/$sub/retro_v$ver\_1/LINE/$sub.LINE.SR.PE.calls
    bedtools merge -i /$outpath/$sub/retro_v$ver\_1/LINE/$sub.LINE.nofilter.temp.bed -c 1 -o count | awk '{if ($4 >1) print}' >> $outpath/$sub/retro_v$ver\_1/LINE/$sub.LINE.SR.PE.calls
    mv $outpath/$sub/retro_v$ver\_0/ALU/$sub.ALU.SR.PE.calls $outpath/$sub/retro_v$ver\_0/ALU/$sub.ALU.SR.PE.old.calls
    awk '{if ($5 > 1) print}' $outpath/$sub/retro_v$ver\_0/ALU/$sub.ALU.SR.PE.old.calls > $outpath/$sub/retro_v$ver\_0/ALU/$sub.ALU.SR.PE.calls
    bedtools merge -i /$outpath/$sub/retro_v$ver\_0/ALU/$sub.ALU.nofilter.temp.bed -c 1 -o count | awk '{if ($4 >1) print}' >> $outpath/$sub/retro_v$ver\_0/ALU/$sub.ALU.SR.PE.calls
    mv $outpath/$sub/retro_v$ver\_1/ALU/$sub.ALU.SR.PE.calls $outpath/$sub/retro_v$ver\_1/ALU/$sub.ALU.SR.PE.old.calls
    awk '{if ($5 > 1) print}' $outpath/$sub/retro_v$ver\_1/ALU/$sub.ALU.SR.PE.old.calls > $outpath/$sub/retro_v$ver\_1/ALU/$sub.ALU.SR.PE.calls
    bedtools merge -i /$outpath/$sub/retro_v$ver\_1/ALU/$sub.ALU.nofilter.temp.bed -c 1 -o count | awk '{if ($4 >1) print}' >> $outpath/$sub/retro_v$ver\_1/ALU/$sub.ALU.SR.PE.calls
    mv $outpath/$sub/retro_v$ver\_0/SVA/$sub.SVA.SR.PE.calls $outpath/$sub/retro_v$ver\_0/SVA/$sub.SVA.SR.PE.old.calls
    awk '{if ($5 > 1) print}' $outpath/$sub/retro_v$ver\_0/SVA/$sub.SVA.SR.PE.old.calls > $outpath/$sub/retro_v$ver\_0/SVA/$sub.SVA.SR.PE.calls
    bedtools merge -i /$outpath/$sub/retro_v$ver\_0/SVA/$sub.SVA.nofilter.temp.bed -c 1 -o count | awk '{if ($4 >1) print}' >> $outpath/$sub/retro_v$ver\_0/SVA/$sub.SVA.SR.PE.calls
    mv $outpath/$sub/retro_v$ver\_1/SVA/$sub.SVA.SR.PE.calls $outpath/$sub/retro_v$ver\_1/SVA/$sub.SVA.SR.PE.old.calls
    awk '{if ($5 > 1) print}' $outpath/$sub/retro_v$ver\_1/SVA/$sub.SVA.SR.PE.old.calls > $outpath/$sub/retro_v$ver\_1/SVA/$sub.SVA.SR.PE.calls
    bedtools merge -i /$outpath/$sub/retro_v$ver\_1/SVA/$sub.SVA.nofilter.temp.bed -c 1 -o count | awk '{if ($4 >1) print}' >> $outpath/$sub/retro_v$ver\_1/SVA/$sub.SVA.SR.PE.calls
    sort -k1,1V -k2,2n -u -o $outpath/$sub/retro_v$ver\_0/LINE/$sub.LINE.SR.PE.calls $outpath/$sub/retro_v$ver\_0/LINE/$sub.LINE.SR.PE.calls
    sort -k1,1V -k2,2n -u -o $outpath/$sub/retro_v$ver\_1/LINE/$sub.LINE.SR.PE.calls $outpath/$sub/retro_v$ver\_1/LINE/$sub.LINE.SR.PE.calls
    sort -k1,1V -k2,2n -u -o $outpath/$sub/retro_v$ver\_0/ALU/$sub.ALU.SR.PE.calls $outpath/$sub/retro_v$ver\_0/ALU/$sub.ALU.SR.PE.calls
    sort -k1,1V -k2,2n -u -o $outpath/$sub/retro_v$ver\_1/ALU/$sub.ALU.SR.PE.calls $outpath/$sub/retro_v$ver\_1/ALU/$sub.ALU.SR.PE.calls
    sort -k1,1V -k2,2n -u -o $outpath/$sub/retro_v$ver\_0/SVA/$sub.SVA.SR.PE.calls $outpath/$sub/retro_v$ver\_0/SVA/$sub.SVA.SR.PE.calls
    sort -k1,1V -k2,2n -u -o $outpath/$sub/retro_v$ver\_1/SVA/$sub.SVA.SR.PE.calls $outpath/$sub/retro_v$ver\_1/SVA/$sub.SVA.SR.PE.calls
else
    grep "L1HS" $outpath/$subject/retro_v$ver/$subject.0.discover | awk '{if ($8 >= 90) print}' | awk '{if ($11 <= 6000 || $12 <=6000) print $1"\t"$2"\t"$3"\t"$5}' | awk '{if (length($1) < 6) print}' | sort -k1,1V -k2,2n -u -o /$outpath/$subject/retro_v$ver\_0/LINE/$subject.LINE.nofilter.0.bed
    grep "L1HS" $outpath/$subject/retro_v$ver/$subject.1.discover | awk '{if ($8 >= 90) print}' | awk '{if ($11 <= 6000 || $12 <=6000) print $1"\t"$2"\t"$3"\t"$5}' | awk '{if (length($1) < 6) print}' | sort -k1,1V -k2,2n -u -o /$outpath/$subject/retro_v$ver\_1/LINE/$subject.LINE.nofilter.1.bed
    grep "Alu" $outpath/$subject/retro_v$ver/$subject.0.discover | awk '{if ($8 >= 90) print $1"\t"$2"\t"$3"\t"$5}' | awk '{if (length($1) < 6) print}' | sort -k1,1V -k2,2n -u -o /$outpath/$subject/retro_v$ver\_0/ALU/$subject.ALU.nofilter.0.bed
    grep "Alu" $outpath/$subject/retro_v$ver/$subject.1.discover | awk '{if ($8 >= 90) print $1"\t"$2"\t"$3"\t"$5}' | awk '{if (length($1) < 6) print}' | sort -k1,1V -k2,2n -u -o /$outpath/$subject/retro_v$ver\_1/ALU/$subject.ALU.nofilter.1.bed
    grep "SVA" $outpath/$subject/retro_v$ver/$subject.0.discover | awk '{if ($8 >= 90) print $1"\t"$2"\t"$3"\t"$5}' | awk '{if (length($1) < 6) print}' | sort -k1,1V -k2,2n -u -o /$outpath/$subject/retro_v$ver\_0/SVA/$subject.SVA.nofilter.0.bed
    grep "SVA" $outpath/$subject/retro_v$ver/$subject.1.discover | awk '{if ($8 >= 90) print $1"\t"$2"\t"$3"\t"$5}' | awk '{if (length($1) < 6) print}' | sort -k1,1V -k2,2n -u -o /$outpath/$subject/retro_v$ver\_1/SVA/$subject.SVA.nofilter.1.bed
    mv $outpath/$subject/retro_v$ver\_0/LINE/$subject.LINE.SR.PE.calls $outpath/$subject/retro_v$ver\_0/LINE/$subject.LINE.SR.PE.old.calls
    awk '{if ($5 > 1) print}' $outpath/$subject/retro_v$ver\_0/LINE/$subject.LINE.SR.PE.old.calls > $outpath/$subject/retro_v$ver\_0/LINE/$subject.LINE.SR.PE.calls
    bedtools merge -i /$outpath/$subject/retro_v$ver\_0/LINE/$subject.LINE.nofilter.0.bed -c 1 -o count | awk '{if ($4 >1) print}' >> $outpath/$subject/retro_v$ver\_0/LINE/$subject.LINE.SR.PE.calls
    mv $outpath/$subject/retro_v$ver\_1/LINE/$subject.LINE.SR.PE.calls $outpath/$subject/retro_v$ver\_1/LINE/$subject.LINE.SR.PE.old.calls
    awk '{if ($5 > 1) print}' $outpath/$subject/retro_v$ver\_1/LINE/$subject.LINE.SR.PE.old.calls > $outpath/$subject/retro_v$ver\_1/LINE/$subject.LINE.SR.PE.calls
    bedtools merge -i /$outpath/$subject/retro_v$ver\_1/LINE/$subject.LINE.nofilter.1.bed -c 1 -o count | awk '{if ($4 >1) print}' >> $outpath/$subject/retro_v$ver\_1/LINE/$subject.LINE.SR.PE.calls
    mv $outpath/$subject/retro_v$ver\_0/ALU/$subject.ALU.SR.PE.calls $outpath/$subject/retro_v$ver\_0/ALU/$subject.ALU.SR.PE.old.calls
    awk '{if ($5 > 1) print}' $outpath/$subject/retro_v$ver\_0/ALU/$subject.ALU.SR.PE.old.calls > $outpath/$subject/retro_v$ver\_0/ALU/$subject.ALU.SR.PE.calls
    bedtools merge -i /$outpath/$subject/retro_v$ver\_0/ALU/$subject.ALU.nofilter.0.bed -c 1 -o count | awk '{if ($4 >1) print}' >> $outpath/$subject/retro_v$ver\_0/ALU/$subject.ALU.SR.PE.calls
    mv $outpath/$subject/retro_v$ver\_1/ALU/$subject.ALU.SR.PE.calls $outpath/$subject/retro_v$ver\_1/ALU/$subject.ALU.SR.PE.old.calls
    awk '{if ($5 > 1) print}' $outpath/$subject/retro_v$ver\_1/ALU/$subject.ALU.SR.PE.old.calls > $outpath/$subject/retro_v$ver\_1/ALU/$subject.ALU.SR.PE.calls
    bedtools merge -i /$outpath/$subject/retro_v$ver\_1/ALU/$subject.ALU.nofilter.1.bed -c 1 -o count | awk '{if ($4 >1) print}' >> $outpath/$subject/retro_v$ver\_1/ALU/$subject.ALU.SR.PE.calls
    mv $outpath/$subject/retro_v$ver\_0/SVA/$subject.SVA.SR.PE.calls $outpath/$subject/retro_v$ver\_0/SVA/$subject.SVA.SR.PE.old.calls
    awk '{if ($5 > 1) print}' $outpath/$subject/retro_v$ver\_0/SVA/$subject.SVA.SR.PE.old.calls > $outpath/$subject/retro_v$ver\_0/SVA/$subject.SVA.SR.PE.calls
    bedtools merge -i /$outpath/$subject/retro_v$ver\_0/SVA/$subject.SVA.nofilter.0.bed -c 1 -o count | awk '{if ($4 >1) print}' >> $outpath/$subject/retro_v$ver\_0/SVA/$subject.SVA.SR.PE.calls
    mv $outpath/$subject/retro_v$ver\_1/SVA/$subject.SVA.SR.PE.calls $outpath/$subject/retro_v$ver\_1/SVA/$subject.SVA.SR.PE.old.calls
    awk '{if ($5 > 1) print}' $outpath/$subject/retro_v$ver\_1/SVA/$subject.SVA.SR.PE.old.calls > $outpath/$subject/retro_v$ver\_1/SVA/$subject.SVA.SR.PE.calls
    bedtools merge -i /$outpath/$subject/retro_v$ver\_1/SVA/$subject.SVA.nofilter.1.bed -c 1 -o count | awk '{if ($4 >1) print}' >> $outpath/$subject/retro_v$ver\_1/SVA/$subject.SVA.SR.PE.calls
    sort -k1,1V -k2,2n -u -o $outpath/$subject/retro_v$ver\_0/LINE/$subject.LINE.SR.PE.calls $outpath/$subject/retro_v$ver\_0/LINE/$subject.LINE.SR.PE.calls
    sort -k1,1V -k2,2n -u -o $outpath/$subject/retro_v$ver\_1/LINE/$subject.LINE.SR.PE.calls $outpath/$subject/retro_v$ver\_1/LINE/$subject.LINE.SR.PE.calls
    sort -k1,1V -k2,2n -u -o $outpath/$subject/retro_v$ver\_0/ALU/$subject.ALU.SR.PE.calls $outpath/$subject/retro_v$ver\_0/ALU/$subject.ALU.SR.PE.calls
    sort -k1,1V -k2,2n -u -o $outpath/$subject/retro_v$ver\_1/ALU/$subject.ALU.SR.PE.calls $outpath/$subject/retro_v$ver\_1/ALU/$subject.ALU.SR.PE.calls
    sort -k1,1V -k2,2n -u -o $outpath/$subject/retro_v$ver\_0/SVA/$subject.SVA.SR.PE.calls $outpath/$subject/retro_v$ver\_0/SVA/$subject.SVA.SR.PE.calls
    sort -k1,1V -k2,2n -u -o $outpath/$subject/retro_v$ver\_1/SVA/$subject.SVA.SR.PE.calls $outpath/$subject/retro_v$ver\_1/SVA/$subject.SVA.SR.PE.calls
fi


