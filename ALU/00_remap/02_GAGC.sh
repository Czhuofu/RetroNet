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
awk '{if ($11<215 && $12>218 && ($4 ~ /Alu/) && (!($4 == "AluYa5a2"))) {split($13, str, ""); print $0"\t+\t"str[215-$11+$9]str[216-$11+$9]str[217-$11+$9]str[218-$11+$9];}}' \
   $datapath/$Retro/$sub.alleles.discover \
   > $datapath/$Retro/$sub.Alu.discover
awk '{if ($11>218 && $12<215 && ($4 ~ /Alu/) && (!($4 == "AluYa5a2"))) {split($13, str, ""); print $0"\t-\t"str[$11+$9-217]str[$11+$9-216]str[$11+$9-215]str[$11+$9-214];}}' \
   $datapath/$Retro/$sub.alleles.discover \
   >> $datapath/$Retro/$sub.Alu.discover

### PE reads, AluYa5a2 ###
awk '{if ($11<216 && $12>219 && ($4 == "AluYa5a2")) {split($13, str, ""); print $0"\t+\t"str[216-$11+$9]str[217-$11+$9]str[218-$11+$9]str[219-$11+$9];}}' \
   $datapath/$Retro/$sub.alleles.discover \
   >>  $datapath/$Retro/$sub.Alu.discover
awk '{if ($11>219 && $12<216 && ($4 == "AluYa5a2")) {split($13, str, ""); print $0"\t-\t"str[$11+$9-218]str[$11+$9-217]str[$11+$9-216]str[$11+$9-215];}}' \
   $datapath/$Retro/$sub.alleles.discover \
   >>  $datapath/$Retro/$sub.Alu.discover

awk '{if (($27 == "+") && ($28 == "GAGC")) print}' $datapath/$Retro/$sub.Alu.discover \
   > $datapath/$Retro/$sub.discover
awk '{if (($27 == "-") && ($28 == "GCTC")) print}' $datapath/$Retro/$sub.Alu.discover \
   >> $datapath/$Retro/$sub.discover

### split reads, no AluYa5a2 ###
awk '{if ($14<215 && $15>218 && ($4 ~ /Alu/) && (!($4 == "AluYa5a2"))) {split($17, str, ""); print $0"\t+\t"str[215-$14-$12]str[216-$14-$12]str[217-$14-$12]str[218-$14-$12];}}' \
   $datapath/$Retro/$sub.sr.alleles.discover \
   > $datapath/$Retro/$sub.sr.Alu.discover

awk '{if ($14>218 && $15<215 && ($4 ~ /Alu/) && (!($4 == "AluYa5a2"))) {split($17, str, ""); print $0"\t-\t"str[$14+$12-217]str[$14+$12-216]str[$14+$12-215]str[$14+$12-214];}}' \
   $datapath/$Retro/$sub.sr.alleles.discover \
   >> $datapath/$Retro/$sub.sr.Alu.discover

### SR reads, AluYa5a2 ###
awk '{if ($14<216 && $15>219 && ($4 == "AluYa5a2")) {split($17, str, ""); print $0"\t+\t"str[216-$14-$12]str[217-$14-$12]str[218-$14-$12]str[219-$14-$12];}}' \
   $datapath/$Retro/$sub.sr.alleles.discover \
   >> $datapath/$Retro/$sub.sr.Alu.discover

awk '{if ($14>219 && $15<216 && ($4 == "AluYa5a2")) {split($17, str, ""); print $0"\t-\t"str[$14+$12-218]str[$14+$12-217]str[$14+$12-216]str[$14+$12-215];}}' \
   $datapath/$Retro/$sub.sr.alleles.discover \
   >> $datapath/$Retro/$sub.sr.Alu.discover

awk '{if (($27=="+") && ($28 == "GAGC")) print}' $datapath/$Retro/$sub.sr.Alu.discover \
   > $datapath/$Retro/$sub.sr.discover
awk '{if (($27=="-") && ($28 == "GCTC")) print}' $datapath/$Retro/$sub.sr.Alu.discover \
   >> $datapath/$Retro/$sub.sr.discover

