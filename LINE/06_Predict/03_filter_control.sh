# module load bedtools
workfolder=$1
ver=$2
sub=$3
cont=$4
cutoff=$5
strand=$6
masterpath=$7
hg=$8
Retro=retro_v$ver\_$strand
TE=LINE

windowBed -v \
   -w 100 \
   -a $workfolder/$sub/$Retro/$TE/RS$cutoff\.$strand\strand.bed \
   -b $workfolder/$cont/$Retro/$TE/$cont\.$TE.SR.PE.calls \
   | windowBed -v \
       -w 100 \
       -a stdin \
       -b $masterpath/refTE/1kgenome/$TE.1kgenome.$hg.$strand.bed \
       | sort -T $workfolder/$sub/$Retro -k1,1 -k2,3n \
           > $workfolder/$sub/$Retro/$TE/noCONT.cutoff$cutoff.calls
