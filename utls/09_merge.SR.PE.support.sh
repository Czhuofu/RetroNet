#!/bin/sh

workpath=$5

sed 's/\;/\|/g' $2/$3/$4/$1.$4.SR.calls > $2/$3/$4/$1.$4.merge.bed

awk '{printf $1"\t"$4"\t"$5"\tpe,"$9","$6","$7","; for (i=10; i<=NF; ++i) {printf "%s,", $i; if ($i == $NF) printf "\n"}}' $2/$3/$4/$1.$4.PE.refine.calls >> $2/$3/$4/$1.$4.merge.bed   

sort -k1,1 -k2,3n  $2/$3/$4/$1.$4.merge.bed | \
 mergeBed -d 50 -c 4 -o distinct -delim ";" \
 -i stdin \
 > $2/$3/$4/$1.$4.SR.PE.calls.test

perl $workpath/utls/10_parse.SR.PE.pl $1 $2 $3 $4

rm $2/$3/$4/$1.$4.merge.bed
rm $2/$3/$4/$1.$4.SR.PE.calls.test
