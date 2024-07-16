#!/bin/sh

module load samtools
module load exonerate
module load bcftools
module load bedtools

workpath="/home/xwzhu/levinson/retro_camal2/common"

if [ "$4" == "ALU" ]; then
	matchTE="Alu" 
elif [ "$4" == "LINE" ]; then
	matchTE="L1"
fi

sed 's/\;/\|/g' $2/$3/$4/$1.$4.SR.calls > $2/$3/$4/$1.$4.merge50.bed

awk '{printf $1"\t"$4"\t"$5"\tpe,"$9","$6","$7","; for (i=10; i<=NF; ++i) {printf "%s,", $i; if ($i == $NF) printf "\n"}}' $2/$3/$4/$1.$4.PE.refine.calls >> $2/$3/$4/$1.$4.merge50.bed   

sort -k1,1 -k2,3n  $2/$3/$4/$1.$4.merge50.bed | \
 mergeBed -d 50 -c 4 -o distinct -delim ";" \
 -i stdin \
 > $2/$3/$4/$1.$4.SR.PE.calls.test50

perl $workpath/20_parse.SR.PE.pl $1 $2 $3 $4

#rm $2/$3/$4/$1.$4.merge50.bed
#rm $2/$3/$4/$1.$4.SR.PE.calls.test
