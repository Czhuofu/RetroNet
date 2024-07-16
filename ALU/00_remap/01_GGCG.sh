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
awk '{if ($11<196 && $12>199 && ($4 ~ /Alu/) && (!($4 == "AluYa5a2"))) {split($13, str, ""); print $0"\t+\t"str[196-$11+$9]str[197-$11+$9]str[198-$11+$9]str[199-$11+$9];}}' \
   $datapath/$Retro/$sub.alleles.discover \
   > $datapath/$Retro/$sub.Alu.discover
awk '{if ($11>199 && $12<196 && ($4 ~ /Alu/) && (!($4 == "AluYa5a2"))) {split($13, str, ""); print $0"\t-\t"str[$11+$9-198]str[$11+$9-197]str[$11+$9-196]str[$11+$9-195];}}' \
   $datapath/$Retro/$sub.alleles.discover \
   >> $datapath/$Retro/$sub.Alu.discover

### PE reads, AluYa5a2 ###
awk '{if ($11<197 && $12>200 && ($4 == "AluYa5a2")) {split($13, str, ""); print $0"\t+\t"str[197-$11+$9]str[198-$11+$9]str[199-$11+$9]str[200-$11+$9];}}' \
   $datapath/$Retro/$sub.alleles.discover \
   >>  $datapath/$Retro/$sub.Alu.discover
awk '{if ($11>200 && $12<197 && ($4 == "AluYa5a2")) {split($13, str, ""); print $0"\t-\t"str[$11+$9-199]str[$11+$9-198]str[$11+$9-197]str[$11+$9-196];}}' \
   $datapath/$Retro/$sub.alleles.discover \
   >>  $datapath/$Retro/$sub.Alu.discover

awk '{if (($27 == "+") && ($28 == "GGCG")) print}' $datapath/$Retro/$sub.Alu.discover \
   > $datapath/$Retro/$sub.discover
awk '{if (($27 == "-") && ($28 == "CGCC")) print}' $datapath/$Retro/$sub.Alu.discover \
   >> $datapath/$Retro/$sub.discover

### split reads, no AluYa5a2 ###
awk '{if ($14<196 && $15>199 && ($4 ~ /Alu/) && (!($4 == "AluYa5a2"))) {split($17, str, ""); print $0"\t+\t"str[196-$14-$12]str[197-$14-$12]str[198-$14-$12]str[199-$14-$12];}}' \
   $datapath/$Retro/$sub.sr.alleles.discover \
   > $datapath/$Retro/$sub.sr.Alu.discover

awk '{if ($14>199 && $15<196 && ($4 ~ /Alu/) && (!($4 == "AluYa5a2"))) {split($17, str, ""); print $0"\t-\t"str[$14+$12-198]str[$14+$12-197]str[$14+$12-196]str[$14+$12-195];}}' \
   $datapath/$Retro/$sub.sr.alleles.discover \
   >> $datapath/$Retro/$sub.sr.Alu.discover

### SR reads, AluYa5a2 ###
awk '{if ($14<197 && $15>200 && ($4 == "AluYa5a2")) {split($17, str, ""); print $0"\t+\t"str[197-$14-$12]str[198-$14-$12]str[199-$14-$12]str[200-$14-$12];}}' \
   $datapath/$Retro/$sub.sr.alleles.discover \
   >> $datapath/$Retro/$sub.sr.Alu.discover

awk '{if ($14>200 && $15<197 && ($4 == "AluYa5a2")) {split($17, str, ""); print $0"\t-\t"str[$14+$12-199]str[$14+$12-198]str[$14+$12-197]str[$14+$12-196];}}' \
   $datapath/$Retro/$sub.sr.alleles.discover \
   >> $datapath/$Retro/$sub.sr.Alu.discover

awk '{if (($27=="+") && ($28 == "GGCG")) print}' $datapath/$Retro/$sub.sr.Alu.discover \
   > $datapath/$Retro/$sub.sr.discover
awk '{if (($27=="-") && ($28 == "CGCC")) print}' $datapath/$Retro/$sub.sr.Alu.discover \
   >> $datapath/$Retro/$sub.sr.discover
