masterpath=/home/xwzhu/levinson/RetroSom
$masterpath/exonerate/exonerate \
   --bestn 1 --model affine:local -s 10 \
   --ryo "INFO: %qi %qal %pi %tS %ti %qab %qae %tab %tae %tas\n" \
   $1 $2 > $3

