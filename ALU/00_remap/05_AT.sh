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
awk '{if ($11<76 && $12>86 && ($4 ~ /Alu/) ) {split($13, str, ""); print $0"\t+\t"str[76-$11+$9]str[86-$11+$9];}}' \
   $datapath/$Retro/$sub.alleles.discover \
   > $datapath/$Retro/$sub.Alu.discover
awk '{if ($11>86 && $12<76 && ($4 ~ /Alu/) ) {split($13, str, ""); print $0"\t-\t"str[$11+$9-85]str[$11+$9-75];}}' \
   $datapath/$Retro/$sub.alleles.discover \
   >> $datapath/$Retro/$sub.Alu.discover

awk '{if (($27 == "+") && ($28 == "AT")) print}' $datapath/$Retro/$sub.Alu.discover \
   > $datapath/$Retro/$sub.discover
awk '{if (($27 == "-") && ($28 == "AT")) print}' $datapath/$Retro/$sub.Alu.discover \
   >> $datapath/$Retro/$sub.discover

### split reads, no AluYa5a2 ###
awk '{if ($14<76 && $15>86 && ($4 ~ /Alu/)) {split($17, str, ""); print $0"\t+\t"str[76-$14+$12]str[86-$14+$12];}}' \
   $datapath/$Retro/$sub.sr.alleles.discover \
   > $datapath/$Retro/$sub.sr.Alu.discover

awk '{if ($14>86 && $15<76 && ($4 ~ /Alu/) ) {split($17, str, ""); print $0"\t-\t"str[$14+$12-85]str[$14+$12-75];}}' \
   $datapath/$Retro/$sub.sr.alleles.discover \
   >> $datapath/$Retro/$sub.sr.Alu.discover

awk '{if (($27=="+") && ($28 == "AT")) print}' $datapath/$Retro/$sub.sr.Alu.discover \
   > $datapath/$Retro/$sub.sr.discover
awk '{if (($27=="-") && ($28 == "AT")) print}' $datapath/$Retro/$sub.sr.Alu.discover \
   >> $datapath/$Retro/$sub.sr.discover

