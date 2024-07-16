### do not allow gaps in ALU alignments ###
### allow gaps in L1 alignments ###
sub=$1
datapath=$2
Retro=$3
echo "$sub $datapath $Retro"

#echo before comment
#: <<'END'
date '+%m/%d/%y %H:%M:%S'
echo "Filtering support reads"
mv $datapath/$Retro/$sub.discover \
   $datapath/$Retro/$sub.alleles.discover
mv $datapath/$Retro/$sub.sr.discover \
   $datapath/$Retro/$sub.sr.alleles.discover

### PE reads, no AluYa5a2 ####
awk '{if ($11<253 && $12>258 && ($4 ~ /AluYb/)) {split($13, str, ""); print $0"\t+\t"str[253-$11+$9]str[254-$11+$9]str[255-$11+$9]str[256-$11+$9]str[257-$11+$9]str[258-$11+$9];}}' \
   $datapath/$Retro/$sub.alleles.discover \
   > $datapath/$Retro/$sub.Alu.discover
awk '{if ($11>258 && $12<253 && ($4 ~ /AluYb/)) {split($13, str, ""); print $0"\t-\t"str[$11+$9-257]str[$11+$9-256]str[$11+$9-255]str[$11+$9-254]str[$11+$9-253]str[$11+$9-252];}}' \
   $datapath/$Retro/$sub.alleles.discover \
   >> $datapath/$Retro/$sub.Alu.discover

awk '{if (($27 == "+") && ($28 == "CAGTCC")) print}' $datapath/$Retro/$sub.Alu.discover \
   > $datapath/$Retro/$sub.discover
awk '{if (($27 == "-") && ($28 == "GGACTG")) print}' $datapath/$Retro/$sub.Alu.discover \
   >> $datapath/$Retro/$sub.discover

### split reads, no AluYa5a2 ###
awk '{if ($14<253 && $15>258 && ($4 ~ /AluYb/)) {split($17, str, ""); print $0"\t+\t"str[253-$14-$12]str[254-$14-$12]str[255-$14-$12]str[256-$14-$12]str[257-$14-$12]str[258-$14-$12];}}' \
   $datapath/$Retro/$sub.sr.alleles.discover \
   > $datapath/$Retro/$sub.sr.Alu.discover
awk '{if ($14>258 && $15<253 && ($4 ~ /AluYb/)) {split($17, str, ""); print $0"\t-\t"str[$14+$12-257]str[$14+$12-256]str[$14+$12-255]str[$14+$12-254]str[$14+$12-253]str[$14+$12-252];}}' \
   $datapath/$Retro/$sub.sr.alleles.discover \
   >> $datapath/$Retro/$sub.sr.Alu.discover

awk '{if (($27=="+") && ($28 == "CAGTCC")) print}' $datapath/$Retro/$sub.sr.Alu.discover \
   > $datapath/$Retro/$sub.sr.discover
awk '{if (($27=="-") && ($28 == "GGACTG")) print}' $datapath/$Retro/$sub.sr.Alu.discover \
   >> $datapath/$Retro/$sub.sr.discover

