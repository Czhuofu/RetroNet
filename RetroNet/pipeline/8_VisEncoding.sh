#!/bin/bash

usage="$(basename "$0") [-h] [-o j m v p g c z x t f] -- Support Reads PNG Visualizing

where:
    -h  show this help text
    -o  output folder path
    -j  subject ID
    -m  masterpath
    -v  version control for RetroSom (default 3)
    -p  slurm partition name
    -g  reference genome (default hg38, supporting hg38, hg19 and b37)
    -c  control or reference group
    -z  gpu partition name (default is N)
    -x  cutoff or threshold (default 0.5)
    -t  TE class: LINE/ALU
    -f  TE family: L1HS/AluYa5"

while getopts ":ho:j:m:v:p:g:c:z:x:t:f:" opt; do
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
    p) partition_name="$OPTARG"
       ;;
    g) hg="$OPTARG"
       ;;
    c) cont="$OPTARG"
       ;;
    z) gpu_partition="$OPTARG"
       ;;
    x) cutoff="$OPTARG"
       ;;
    t) TEclass="$OPTARG"
       ;;
    f) TEfamily="$OPTARG"
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

# Activate RetroSom virtual env #
source /opt/miniconda3/bin/activate RetroSom

mkdir $outpath/$sub/visual_$ver
mkdir $outpath/$sub/calls_$ver
mkdir $outpath/$sub/visual_$ver/$TEclass

if [ "$hg" == hg38 ] && { [ "$TEclass" == LINE ] || [ "$TEclass" == SVA ]; }
then
awk '{if ($5 > 1) print}' $outpath/$sub/retro_v${ver}_1/$TEclass/$sub.$TEclass.SR.PE.calls \
  | grep -v 'chrUn' | grep -v 'random' | grep -v 'alt' \
  | windowBed -v -w 200 \
  -a stdin \
  -b $outpath/$cont/retro_v${ver}_1/$TEclass/$cont.$TEclass.SR.PE.calls \
  | windowBed -v -w 50 \
  -b $masterpath/RetroNet/hg38.low10div.TE.bed \
  -a stdin > $outpath/$sub/calls_$ver/$sub.$TEclass.1.calls

awk '{if ($5 > 1) print}' $outpath/$sub/retro_v${ver}_0/$TEclass/$sub.$TEclass.SR.PE.calls \
  | grep -v 'chrUn' | grep -v 'random' | grep -v 'alt' \
  | windowBed -v -w 200 \
  -a stdin \
  -b $outpath/$cont/retro_v${ver}_0/$TEclass/$cont.$TEclass.SR.PE.calls \
  | windowBed -v -w 50 \
  -b $masterpath/RetroNet/hg38.low10div.TE.bed \
  -a stdin > $outpath/$sub/calls_$ver/$sub.$TEclass.0.calls
fi

if [ "$hg" == hg38 ] && [ "$TEclass" == ALU ]
then
awk '{if ($5 > 1) print}' $outpath/$sub/retro_v${ver}_1/$TEclass/$sub.$TEclass.SR.PE.calls \
  | grep -v 'chrUn' | grep -v 'random' | grep -v 'alt' \
  | windowBed -v -w 200 \
  -a stdin \
  -b $outpath/$cont/retro_v${ver}_1/$TEclass/$cont.$TEclass.SR.PE.calls \
  | windowBed -v -w 50 \
  -b $masterpath/RetroNet/hg38.low10div.TE.alu20.bed \
  -a stdin > $outpath/$sub/calls_$ver/$sub.$TEclass.1.calls

awk '{if ($5 > 1) print}' $outpath/$sub/retro_v${ver}_0/$TEclass/$sub.$TEclass.SR.PE.calls \
  | grep -v 'chrUn' | grep -v 'random' | grep -v 'alt' \
  | windowBed -v -w 200 \
  -a stdin \
  -b $outpath/$cont/retro_v${ver}_0/$TEclass/$cont.$TEclass.SR.PE.calls \
  | windowBed -v -w 50 \
  -b $masterpath/RetroNet/hg38.low10div.TE.alu20.bed \
  -a stdin > $outpath/$sub/calls_$ver/$sub.$TEclass.0.calls
fi

if [ "$hg" == b37 ] && { [ "$TEclass" == LINE ] || [ "$TEclass" == SVA ]; }
then
awk '{if ($5 > 1) print}' $outpath/$sub/retro_v${ver}_1/$TEclass/$sub.$TEclass.SR.PE.calls \
  | grep -v 'chrUn' | grep -v 'random' | grep -v 'alt' \
  | windowBed -v -w 200 \
  -a stdin \
  -b $outpath/$cont/retro_v${ver}_1/$TEclass/$cont.$TEclass.SR.PE.calls \
  | windowBed -v -w 50 \
  -b $masterpath/RetroNet/b37.low10div.TE.bed \
  -a stdin > $outpath/$sub/calls_$ver/$sub.$TEclass.1.calls

awk '{if ($5 > 1) print}' $outpath/$sub/retro_v${ver}_0/$TEclass/$sub.$TEclass.SR.PE.calls \
  | grep -v 'chrUn' | grep -v 'random' | grep -v 'alt' \
  | windowBed -v -w 200 \
  -a stdin \
  -b $outpath/$cont/retro_v${ver}_0/$TEclass/$cont.$TEclass.SR.PE.calls \
  | windowBed -v -w 50 \
  -b $masterpath/RetroNet/b37.low10div.TE.bed \
  -a stdin > $outpath/$sub/calls_$ver/$sub.$TEclass.0.calls
fi

if [ "$hg" == b37 ] && [ "$TEclass" == ALU ]
then
awk '{if ($5 > 1) print}' $outpath/$sub/retro_v${ver}_1/$TEclass/$sub.$TEclass.SR.PE.calls \
  | grep -v 'chrUn' | grep -v 'random' | grep -v 'alt' \
  | windowBed -v -w 200 \
  -a stdin \
  -b $outpath/$cont/retro_v${ver}_1/$TEclass/$cont.$TEclass.SR.PE.calls \
  | windowBed -v -w 50 \
  -b $masterpath/RetroNet/b37.low10div.TE.alu20.bed \
  -a stdin > $outpath/$sub/calls_$ver/$sub.$TEclass.1.calls

awk '{if ($5 > 1) print}' $outpath/$sub/retro_v${ver}_0/$TEclass/$sub.$TEclass.SR.PE.calls \
  | grep -v 'chrUn' | grep -v 'random' | grep -v 'alt' \
  | windowBed -v -w 200 \
  -a stdin \
  -b $outpath/$cont/retro_v${ver}_0/$TEclass/$cont.$TEclass.SR.PE.calls \
  | windowBed -v -w 50 \
  -b $masterpath/RetroNet/b37.low10div.TE.alu20.bed \
  -a stdin > $outpath/$sub/calls_$ver/$sub.$TEclass.0.calls
fi

# --------------------------------------------------------------------------- #
#slurm_sc="-o %x.%A.output -e %x.%A.output -p $partition_name --mem=16gb --time=100:00:00"
#slurm_cpu="-o %x.%A.output -e %x.%A.output -p $partition_name --nodes=1 --ntasks-per-node=1 --cpus-per-task=8 --mem=64gb --time=100:00:00"
#slurm_gpu="-o %x.%A.output -e %x.%A.output -p $gpu_partition --gres=gpu:1 --mem=32gb --time=100:00:00"

###########################################
### Step9: Draw PNG image of Supporting ###
###########################################
## Create PNG plot Job ##
#split_row=1000
# Strand 1 #
#jid9=
#split -l $split_row $outpath/$sub/calls_$ver/$sub.$TEclass.1.calls $sub.$TEclass.calls_$ver.1
#for file in $sub.$TEclass.calls_$ver.1*; do
#   jid_tmp=$(echo -e '#!/bin/sh\n'"singularity exec -B $outpath,$tmppath,$masterpath $masterpath/pipeline/RetroNet.sif ./9_DrawPNG.sh -o $outpath -j $sub -m $masterpath -g $hg -i $file -c $conda_path -t $TEclass -f $TEfamily -v $ver -s 1 -l candidate" | sbatch -J PNG_${TEclass} $slurm_sc | awk '{print $4}')
#   jid9=${jid9}${jid_tmp}
#done
# Strand 0 #
#split -l $split_row $outpath/$sub/calls_$ver/$sub.$TEclass.0.calls $sub.$TEclass.calls_$ver.0
#for file in $sub.$TEclass.calls_$ver.0*; do
#   jid_tmp=$(echo -e '#!/bin/sh\n'"singularity exec -B $outpath,$tmppath,$masterpath $masterpath/pipeline/RetroNet.sif ./9_DrawPNG.sh -o $outpath -j $sub -m $masterpath -g $hg -i $file -c $conda_path -t $TEclass -f $TEfamily -v $ver -s 0 -l candidate" | sbatch -J PNG_${TEclass} $slurm_sc | awk '{print $4}')
#   jid9=${jid9}${jid_tmp}
#done

###############################
### Step10: RetroNet Predict ##
###############################
#if [ "$gpu_partition" != "N" ]
#then
#    jid10=$(echo -e '#!/bin/sh\n' "singularity exec -B $outpath,$tmppath,$masterpath $masterpath/pipeline/RetroNet.sif ./10_RunRetroNet.sh -o $outpath -j $sub -t $TEclass -v $ver -m $masterpath" | sbatch -J ${TEclass}_RetroNet $slurm_gpu --dependency=afterok:$jid9 | awk '{print $4}')
#else
#    jid10=$(echo -e '#!/bin/sh\n' "singularity exec -B $outpath,$tmppath,$masterpath $masterpath/pipeline/RetroNet.sif ./10_RunRetroNet.sh -o $outpath -j $sub -t $TEclass -v $ver -m $masterpath" | sbatch -J ${TEclass}_RetroNet $slurm_cpu --dependency=afterok:$jid9 | awk '{print $4}')
#fi
