sub=$1
workfolder=$2
Retro=retro_$3
chr1=`grep -P "^chr1\t" $workfolder/$sub/$Retro/$sub.discover | wc -l`
chr11=`grep -P "^chr11\t" $workfolder/$sub/$Retro/$sub.discover | wc -l`
chr22=`grep -P "^chr22\t" $workfolder/$sub/$Retro/$sub.discover | wc -l`
chr7=`grep -P "^chr7\t" $workfolder/$sub/$Retro/$sub.discover | wc -l`
chr8=`grep -P "^chr8\t" $workfolder/$sub/$Retro/$sub.discover | wc -l`
chr9=`grep -P "^chr9\t" $workfolder/$sub/$Retro/$sub.discover | wc -l`
total=`wc -l $workfolder/$sub/$Retro/$sub.discover | awk '{print $1}'`
rat1=$(echo "scale=4; $chr1/$total" | bc)
rat11=$(echo "scale=4; $chr11/$total" | bc)
rat22=$(echo "scale=4; $chr22/$total" | bc)
rat7=$(echo "scale=4; $chr7/$total" | bc)
rat8=$(echo "scale=4; $chr8/$total" | bc)
rat9=$(echo "scale=4; $chr9/$total" | bc)
echo "chr1=$rat1"
echo "chr7=$rat7"
echo "chr8=$rat8"
echo "chr9=$rat9"
echo "chr11=$rat11"
echo "chr22=$rat22"
