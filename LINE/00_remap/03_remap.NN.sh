### do not allow gaps in ALU alignments ###
### allow gaps in L1 alignments ###
sub=$1
datapath=$2
Retro=$3
hgX=$4
sr_length=50
geno=$5
cord1=$6
cord2=$7
masterpath=$8
echo "$sub $datapath $Retro $hgX"

#echo before comment
#: <<'END'
date '+%m/%d/%y %H:%M:%S'
echo "Filtering support reads"

mv $datapath/$Retro/$sub.discover \
   $datapath/$Retro/$sub.alleles.discover

mv $datapath/$Retro/$sub.sr.discover \
   $datapath/$Retro/$sub.sr.alleles.discover

### L1HS ###
awk -v c1=$cord1 -v c2=$cord2 '{if ($8>85 && $11<c1 && $12>=c2 && ($4 == "L1HS")) {split($13, str, ""); print $0"\t+\t"str[c1-$11]str[c2-$11];}}' \
   $datapath/$Retro/$sub.alleles.discover \
   > $datapath/$Retro/$sub.L1HS.discover

awk -v c1=$cord1 -v c2=$cord2 '{if ($8>85 && $11>=c2 && $12<c1 && ($4 == "L1HS")) {split($13, str, ""); print $0"\t-\t"str[$11-c2-1]str[$11-c1-1];}}' \
   $datapath/$Retro/$sub.alleles.discover \
   >> $datapath/$Retro/$sub.L1HS.discover

awk -v c1=$cord1 -v c2=$cord2 '{if ($9>85 && $14<c1 && $15>=c2 && ($4 == "L1HS")) {split($17, str, ""); print $0"\t+\t"str[c1-$14]str[c2-$14];}}' \
   $datapath/$Retro/$sub.sr.alleles.discover \
   > $datapath/$Retro/$sub.sr.L1HS.discover

awk -v c1=$cord1 -v c2=$cord2 '{if ($9>85 && $14>=c2 && $15<c1 && ($4 == "L1HS")) {split($17, str, ""); print $0"\t-\t"str[$14-c2-1]str[$14-c1-1];}}' \
   $datapath/$Retro/$sub.sr.alleles.discover \
   >> $datapath/$Retro/$sub.sr.L1HS.discover

### print the fasta file ###
awk '{print ">pe_"$4"_"$5"_"$8"\n"$13}' $datapath/$Retro/$sub.L1HS.discover \
   > $datapath/$Retro/$sub.L1HS.fa

awk '{print ">sr_"$4"_"$5"_"$9"\n"$17}' $datapath/$Retro/$sub.sr.L1HS.discover \
   >> $datapath/$Retro/$sub.L1HS.fa

env TMPDIR=$datapath/$Retro exonerate \
    --bestn 1 --model affine:local --showalignment 1 \
    --ryo "INFO: %qi %qal %pi %tS %ti %qab %qae %tab %tae %qs\n" \
    $datapath/$Retro/$sub.L1HS.fa \
    $masterpath/refTE/sequence/L1HS.fa \
    > $datapath/$Retro/$sub.L1HS.alignout

$masterpath/LINE/00_remap/02_realign.L1HS.pl $sub $Retro $datapath $geno $cord1 $cord2

