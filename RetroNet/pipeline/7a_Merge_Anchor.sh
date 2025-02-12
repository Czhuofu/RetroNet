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

# Activate RetroSom virtual env #
source /opt/miniconda3/bin/activate RetroSom

sub=$subject\_Combined
sub2=$subject\-
startlib=1

### combine the anchor ###
for (( strand=0; strand<=1; strand++ ))
     do
	    echo -e "read\tdirect\tinsert\tseg\talign" > $outpath/$sub/retro_v$ver\_$strand/LINE/split-reads-anchor.txt
	    echo -e "read\tdirect\tinsert\tseg\talign" > $outpath/$sub/retro_v$ver\_$strand/ALU/split-reads-anchor.txt
            echo -e "read\tdirect\tinsert\tseg\talign" > $outpath/$sub/retro_v$ver\_$strand/SVA/split-reads-anchor.txt
            head -n 1 $outpath/$sub2\1/retro_v$ver\_$strand/$sub2\1.pe.LINE.matrix > $outpath/$sub/retro_v$ver\_$strand/$sub.pe.$TE.matrix
	    head -n 1 $outpath/$sub2\1/retro_v$ver\_$strand/$sub2\1.sr.LINE.matrix > $outpath/$sub/retro_v$ver\_$strand/$sub.sr.$TE.matrix
	    head -n 1 $outpath/$sub2\1/retro_v$ver\_$strand/$sub2\1.pe.ALU.matrix > $outpath/$sub/retro_v$ver\_$strand/$sub.pe.$TE.matrix
	    head -n 1 $outpath/$sub2\1/retro_v$ver\_$strand/$sub2\1.sr.ALU.matrix > $outpath/$sub/retro_v$ver\_$strand/$sub.sr.$TE.matrix
     done




for (( c=startlib; c<=$readgroup; c++ ))
do
  for type in LINE ALU SVA
  do
    for strand in 0 1
    do
      awk '{if (!/^read/) print $0}' $outpath/$sub2$c/retro_v${ver}_$strand/$type/split-reads-anchor.txt | \
      awk -v lib=$c '{
        split($1, names, ":");
        output = names[1]":lib"lib"_";
        for (i = 2; i <= length(names); i++) {
          output = output names[i];
          if (i < length(names)) {
            output = output ":";
          }
        }
        print output "\t" $2 "\t" $3 "\t" $4 "\t" $5
      }' >> $outpath/$sub/retro_v${ver}_$strand/$type/split-reads-anchor.txt
    done
  done
done


for (( c=startlib; c<=$readgroup; c++ ))
     do
       for  (( strand=0; strand<=1; strand++ ))
            do
       		awk -v lib=$c '$4="lib"lib"_"$4' OFS='\t' $outpath/$sub2$c/retro_v$ver\_$strand/$sub2$c.pe.LINE.matrix | awk 'FNR > 1' >> $outpath/$sub/retro_v$ver\_$strand/$sub.pe.LINE.matrix
       		awk -v lib=$c '$4="lib"lib"_"$4' OFS='\t' $outpath/$sub2$c/retro_v$ver\_$strand/$sub2$c.sr.LINE.matrix | awk 'FNR > 1' >> $outpath/$sub/retro_v$ver\_$strand/$sub.sr.LINE.matrix
       		awk -v lib=$c '$4="lib"lib"_"$4' OFS='\t' $outpath/$sub2$c/retro_v$ver\_$strand/$sub2$c.pe.ALU.matrix | awk 'FNR > 1' >> $outpath/$sub/retro_v$ver\_$strand/$sub.pe.ALU.matrix
      		awk -v lib=$c '$4="lib"lib"_"$4' OFS='\t' $outpath/$sub2$c/retro_v$ver\_$strand/$sub2$c.sr.ALU.matrix | awk 'FNR > 1' >> $outpath/$sub/retro_v$ver\_$strand/$sub.sr.ALU.matrix
            done     
     done
