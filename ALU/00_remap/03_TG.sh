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
awk '{if ($11<236 && $12>239 && ($4 ~ /Alu/) && (!($4 == "AluYa5a2"))) {split($13, str, ""); print $0"\t+\t"str[236-$11+$9]str[237-$11+$9]str[238-$11+$9]str[239-$11+$9]}}' \
   $datapath/$Retro/$sub.alleles.discover \
   > $datapath/$Retro/$sub.Alu.discover
awk '{if ($11>239 && $12<236 && ($4 ~ /Alu/) && (!($4 == "AluYa5a2"))) {split($13, str, ""); print $0"\t-\t"str[$11+$9-238]str[$11+$9-237]str[$11+$9-236]str[$11+$9-235];}}' \
   $datapath/$Retro/$sub.alleles.discover \
   >> $datapath/$Retro/$sub.Alu.discover

### PE reads, AluYa5a2 ###
awk '{if ($11<237 && $12>240 && ($4 == "AluYa5a2")) {split($13, str, ""); print $0"\t+\t"str[237-$11+$9]str[238-$11+$9]str[239-$11+$9]str[240-$11+$9];}}' \
   $datapath/$Retro/$sub.alleles.discover \
   >>  $datapath/$Retro/$sub.Alu.discover
awk '{if ($11>240 && $12<237 && ($4 == "AluYa5a2")) {split($13, str, ""); print $0"\t-\t"str[$11+$9-239]str[$11+$9-238]str[$11+$9-237]str[$11+$9-236];}}' \
   $datapath/$Retro/$sub.alleles.discover \
   >>  $datapath/$Retro/$sub.Alu.discover

awk '{if (($27 == "+") && ($28 == "TGCG")) print}' $datapath/$Retro/$sub.Alu.discover \
   > $datapath/$Retro/$sub.discover
awk '{if (($27 == "-") && ($28 == "CGCA")) print}' $datapath/$Retro/$sub.Alu.discover \
   >> $datapath/$Retro/$sub.discover

### split reads, no AluYa5a2 ###
awk '{if ($14<236 && $15>239 && ($4 ~ /Alu/) && (!($4 == "AluYa5a2"))) {split($17, str, ""); print $0"\t+\t"str[236-$14-$12]str[237-$14-$12]str[238-$14-$12]str[239-$14-$12]}}' \
   $datapath/$Retro/$sub.sr.alleles.discover \
   > $datapath/$Retro/$sub.sr.Alu.discover
awk '{if ($14>239 && $15<236 && ($4 ~ /Alu/) && (!($4 == "AluYa5a2"))) {split($17, str, ""); print $0"\t-\t"str[$14+$12-238]str[$14+$12-237]str[$14+$12-236]str[$14+$12-235];}}' \
   $datapath/$Retro/$sub.sr.alleles.discover \
   >> $datapath/$Retro/$sub.sr.Alu.discover

### SR reads, AluYa5a2 ###
awk '{if ($14<237 && $15>240 && ($4 == "AluYa5a2")) {split($17, str, ""); print $0"\t+\t"str[237-$14-$12]str[238-$14-$12]str[239-$14-$12]str[240-$14-$12];}}' \
   $datapath/$Retro/$sub.sr.alleles.discover \
   >> $datapath/$Retro/$sub.sr.Alu.discover
awk '{if ($14>240 && $15<237 && ($4 == "AluYa5a2")) {split($17, str, ""); print $0"\t-\t"str[$14+$12-239]str[$14+$12-238]str[$14+$12-237]str[$14+$12-236];}}' \
   $datapath/$Retro/$sub.sr.alleles.discover \
   >> $datapath/$Retro/$sub.sr.Alu.discover

awk '{if (($27=="+") && ($28 == "TGCG")) print}' $datapath/$Retro/$sub.sr.Alu.discover \
   > $datapath/$Retro/$sub.sr.discover
awk '{if (($27=="-") && ($28 == "CGCA")) print}' $datapath/$Retro/$sub.sr.Alu.discover \
   >> $datapath/$Retro/$sub.sr.discover
