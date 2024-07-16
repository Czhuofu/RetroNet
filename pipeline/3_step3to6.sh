#!/bin/bash

usage="$(basename "$0") [-h] [-o j m v g] -- Putative MEIs & Remapping L1HS or AluY specific Alleles

where:
    -h  show this help text
    -o  output folder path
    -j  subject ID
    -m  masterpath
    -v  version control for RetroSom (default 1)
    -g  reference genome (default hg38, supporting hg38, hg19 and b37)"

while getopts ":ho:j:m:v:g:" opt; do
  case $opt in
    h) echo "$usage"
       exit
       ;;
    o) outpath=$(readlink -f $OPTARG)
       ;;
    j) sub="$OPTARG"
       ;;
    m) masterpath=$(readlink -f $OPTARG)
       ;;
    v) ver="$OPTARG"
       ;;
    g) hg="$OPTARG"
       ;;
    :) printf "missing argument for -%s\n" "$OPTARG" >&2
       echo "$usage" >&2
       exit 1
       ;;
    \?) printf "illegal option: -%s\n" "$OPTARG" >&2
       echo "$usage" >&2
       exit 1
       ;;
  esac
done

##load env
source /opt/miniconda3/bin/activate RetroSom
export TMPDIR=$outpath/$sub/temp

# Preprocessing #
bash 3a_Preprocess.sh -o $outpath -j $sub -m $masterpath -v $ver -g $hg
# May not need this wait command #
wait

### Parallel running 3b, 4a, 4b ###
# insertions in -strand #
bash 3b_CallPutative_MEI.sh -o $outpath -j $sub -m $masterpath -v $ver -g $hg -s 0 -f 1 &
# insertions in +strand #
bash 3b_CallPutative_MEI.sh -o $outpath -j $sub -m $masterpath -v $ver -g $hg -s 1 -f 1 &
# Realign L1 supporting reads #
bash 4a_remapL1.sh -o $outpath -j $sub -m $masterpath -v $ver -g $hg &
# Realign Alu supporting reads #
bash 4b_remapAlu.sh -o $outpath -j $sub -m $masterpath -v $ver -g $hg &

wait
##run sva combine
bash 3c_runSVA.sh -o $outpath -i $sub -m $masterpath -r $ver -g $hg -s 0 -f 1 
bash 3c_runSVA.sh -o $outpath -i $sub -m $masterpath -r $ver -g $hg -s 1 -f 1

wait

### Parallel running 5a0, 5b0, 5c0, 5d0 ###
bash 5a_matrix_gen_L1PE.sh -o $outpath -j $sub -m $masterpath -v $ver -g $hg -s 0 &
bash 5b_matrix_gen_L1SR.sh -o $outpath -j $sub -m $masterpath -v $ver -g $hg -s 0 &
bash 5c_matrix_gen_AluPE.sh -o $outpath -j $sub -m $masterpath -v $ver -g $hg -s 0 &
bash 5d_matrix_gen_AluSR.sh -o $outpath -j $sub -m $masterpath -v $ver -g $hg -s 0 &

wait

### Parallel running 5a1, 5b1, 5c1, 5d1 ###
bash 5a_matrix_gen_L1PE.sh -o $outpath -j $sub -m $masterpath -v $ver -g $hg -s 1 &
bash 5b_matrix_gen_L1SR.sh -o $outpath -j $sub -m $masterpath -v $ver -g $hg -s 1 &
bash 5c_matrix_gen_AluPE.sh -o $outpath -j $sub -m $masterpath -v $ver -g $hg -s 1 &
bash 5d_matrix_gen_AluSR.sh -o $outpath -j $sub -m $masterpath -v $ver -g $hg -s 1 &

wait
TE_all=("SVA")
#TE_all=("LINE" "ALU" "SVA")
for ((i=0; i<"${#TE_all[@]}"; i++)); do
   for ((k=0; k<2; k++)); do
       $masterpath/pipeline/3d_merge_anchor.pl $outpath $sub $ver $k $hg $masterpath ${TE_all[i]}
   done
done

### Parallel running 6a, 6b, 6c, 6d ###
#srun --ntasks=1 bash 6a_Level1_L1PE.sh -o $outpath -j $sub -m $masterpath -v $ver &
#srun --ntasks=1 bash 6b_Level1_L1SR.sh -o $outpath -j $sub -m $masterpath -v $ver &
#srun --ntasks=1 bash 6c_Level1_AluPE.sh -o $outpath -j $sub -m $masterpath -v $ver &
#srun --ntasks=1 bash 6d_Level1_AluSR.sh -o $outpath -j $sub -m $masterpath -v $ver &

#wait
