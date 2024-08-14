#!/bin/bash

usage="$(basename "$0") [-h] [-o j m v p g s l c z x] -- Combine Library and Deep-Learning Model Prediction

where:
    -h  show this help text
    -o  output folder path
    -j  subject ID
    -m  masterpath
    -v  version control for RetroSom (default 3)
    -p  slurm partition name
    -g  reference genome (default hg38, supporting hg38, hg19 and b37)
    -l  number of the libraries (default is 1)
    -c  control or reference group
    -z  gpu partition name (default is N)
    -x  cutoff or threshold (default 0.99)"

readgroup=1
ver=3
hg=hg38
gpu_partition="N"
cutoff=0.99
while getopts ":ho:j:m:v:p:g:l:c:z:x:" opt; do
  case $opt in
    h) echo "$usage"
       exit
       ;;
    o) outpath=$(readlink -f $OPTARG)
       ;;
    j) subject="$OPTARG"
       ;;
    m) masterpath=$(readlink -f $OPTARG)
       ;;
    v) ver="$OPTARG"
       ;;
    p) partition_name="$OPTARG"
       ;;
    g) hg="$OPTARG"
       ;;
    l) readgroup="$OPTARG"
       ;;
    c) cont="$OPTARG"
       ;;
    z) gpu_partition="$OPTARG"
       ;;
    x) cutoff="$OPTARG"
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

### updating the location of the reference sequences ###
echo -e "Alu\t$masterpath/refTE/sequence/ALU.fa" > $masterpath/refTE/TE_ALHS.bed
echo -e "L1\t$masterpath/refTE/sequence/L1HS.fa" >> $masterpath/refTE/TE_ALHS.bed
echo -e "HERVK\t$masterpath/refTE/sequence/HERVK.fa" >> $masterpath/refTE/TE_ALHS.bed
echo -e "HERVH\t$masterpath/refTE/sequence/HERVH.fa" >> $masterpath/refTE/TE_ALHS.bed
echo -e "SVA\t$masterpath/refTE/sequence/SVA.fa" >> $masterpath/refTE/TE_ALHS.bed

slurm_sc="-o %x.%A.output -e %x.%A.output -p $partition_name --mem=16gb --time=100:00:00"
slurm_sc_parallel="-o %x.%A.output -e %x.%A.output -p $partition_name --ntasks=4 --cpus-per-task=1 --mem-per-cpu=16gb --time=48:00:00"
slurm_mc="-o %x.%A.output -e %x.%A.output -p $partition_name --time=96:00:00 --nodes=1 --ntasks-per-node=1 --cpus-per-task=10 --mem=150GB"
export TMPDIR=$outpath/$sub/script

###########################################
### Step7: Pairing the supporting reads ###
###########################################
if [ "$readgroup" != 1 ]
then
    # Creating folders #
    sub=${subject}_Combined
    mkdir $outpath/$sub
    mkdir $outpath/$sub/reads
    mkdir $outpath/$sub/align
    mkdir $outpath/$sub/QC
    mkdir $outpath/$sub/temp
    mkdir $outpath/$sub/script
    mkdir $outpath/$sub/retro_v$ver
    cd $outpath/$sub/script
    cp $masterpath/RetroNet/pipeline/* $outpath/$sub/script
    # submit jobs #
    jid7a=$(echo -e '#!/bin/sh\n'"singularity exec -B $outpath,$tmppath,$masterpath $masterpath/pipeline/RetroNet.sif ./7_CombineLib.sh -o $outpath -j $subject -m $masterpath -v $ver -s 0 -l $readgroup -g $hg" | sbatch -J CombineStrand0 $slurm_sc | awk '{print $4}')
    jid7b=$(echo -e '#!/bin/sh\n'"singularity exec -B $outpath,$tmppath,$masterpath $masterpath/pipeline/RetroNet.sif ./7_CombineLib.sh -o $outpath -j $subject -m $masterpath -v $ver -s 1 -l $readgroup -g $hg" | sbatch -J CombineStrand1 $slurm_sc | awk '{print $4}')
    jid7c=$(echo -e '#!/bin/sh\n'"singularity exec -B $outpath,$tmppath,$masterpath $masterpath/pipeline/RetroNet.sif ./7a_Merge_Anchor.sh -o $outpath -j $subject -m $masterpath -v $ver -l $readgroup" | sbatch -J CombineAnchor $slurm_sc --dependency=afterok:$jid7a:$jid7b | awk '{print $4}')
fi

#############################################################
### Step8: Draw PNG image of Supporting Reads and Predict ###
#############################################################
if [ "$readgroup" != 1 ]
then
    jid8a=$(echo -e '#!/bin/sh\n'"singularity exec -B $outpath,$tmppath,$masterpath $masterpath/pipeline/RetroNet.sif ./8_VisEncoding.sh -o $outpath -j $sub -m $masterpath -v $ver -p $partition_name -g $hg -z $gpu_partition -c $cont -x $cutoff -t LINE -f L1HS" | sbatch -J Vis_L1 $slurm_sc --dependency=afterok:$jid7c | awk '{print $4}')
    jid8b=$(echo -e '#!/bin/sh\n'"singularity exec -B $outpath,$tmppath,$masterpath $masterpath/pipeline/RetroNet.sif ./8_VisEncoding.sh -o $outpath -j $sub -m $masterpath -v $ver -p $partition_name -g $hg -z $gpu_partition -c $cont -x $cutoff -t ALU -f AluY" | sbatch -J Vis_Alu $slurm_sc --dependency=afterok:$jid7c| awk '{print $4}')
    jid8c=$(echo -e '#!/bin/sh\n'"singularity exec -B $outpath,$tmppath,$masterpath $masterpath/pipeline/RetroNet.sif ./8_VisEncoding.sh -o $outpath -j $sub -m $masterpath -v $ver -p $partition_name -g $hg -z $gpu_partition -c $cont -x $cutoff -t SVA -f SVA" | sbatch -J Vis_SVA $slurm_sc --dependency=afterok:$jid7c| awk '{print $4}')
    jid8cLINE=$(echo -e '#!/bin/sh\n' "bash ./8c_submit910.sh -o $outpath -j $sub -m $masterpath -v $ver -p $partition_name -g $hg -z $gpu_partition -c $cont -x $cutoff -t LINE -f L1HS" | sbatch -J L19_10 $slurm_sc --dependency=afterok:$jid8a | awk '{print $4}')
    jid8cAlu=$(echo -e '#!/bin/sh\n'"bash ./8c_submit910.sh -o $outpath -j $sub -m $masterpath -v $ver -p $partition_name -g $hg -z $gpu_partition -c $cont -x $cutoff -t ALU -f AluY" | sbatch -J Alu9_10 $slurm_sc --dependency=afterok:$jid8b | awk '{print $4}')
    jid8cSVA=$(echo -e '#!/bin/sh\n'"bash ./8c_submit910.sh -o $outpath -j $sub -m $masterpath -v $ver -p $partition_name -g $hg -z $gpu_partition -c $cont -x $cutoff -t SVA -f SVA" | sbatch -J SVA9_10 $slurm_sc --dependency=afterok:$jid8c| awk '{print $4}')
else
    # Copy script #
    cd $outpath/$subject/script
    cp $masterpath/RetroNet/pipeline/* $outpath/$subject/script
    #
    jid8a=$(echo -e '#!/bin/sh\n'"singularity exec -B $outpath,$tmppath,$masterpath $masterpath/pipeline/RetroNet.sif ./8_VisEncoding.sh -o $outpath -j $subject -m $masterpath -v $ver -p $partition_name -g $hg -z $gpu_partition -c $cont -x $cutoff -t LINE -f L1HS" | sbatch -J Vis_L1 $slurm_sc | awk '{print $4}')
    jid8b=$(echo -e '#!/bin/sh\n'"singularity exec -B $outpath,$tmppath,$masterpath $masterpath/pipeline/RetroNet.sif ./8_VisEncoding.sh -o $outpath -j $subject -m $masterpath -v $ver -p $partition_name -g $hg -z $gpu_partition -c $cont -x $cutoff -t ALU -f AluY" | sbatch -J Vis_Alu $slurm_sc | awk '{print $4}')    
    jid8c=$(echo -e '#!/bin/sh\n'"singularity exec -B $outpath,$tmppath,$masterpath $masterpath/pipeline/RetroNet.sif ./8_VisEncoding.sh -o $outpath -j $subject -m $masterpath -v $ver -p $partition_name -g $hg -z $gpu_partition -c $cont -x $cutoff -t SVA -f SVA" | sbatch -J Vis_SVA $slurm_sc | awk '{print $4}')
    jid8cLINE=$(echo -e '#!/bin/sh\n' "bash ./8c_submit910.sh -o $outpath -j $subject -m $masterpath -v $ver -p $partition_name -g $hg -z $gpu_partition -c $cont -x $cutoff -t LINE -f L1HS" | sbatch -J L19_10 $slurm_sc --dependency=afterok:$jid8a | awk '{print $4}')
    jid8cAlu=$(echo -e '#!/bin/sh\n'"bash ./8c_submit910.sh -o $outpath -j $subject -m $masterpath -v $ver -p $partition_name -g $hg -z $gpu_partition -c $cont -x $cutoff -t ALU -f AluY" | sbatch -J Alu9_10 $slurm_sc --dependency=afterok:$jid8b | awk '{print $4}')
    jid8cSVA=$(echo -e '#!/bin/sh\n'"bash ./8c_submit910.sh -o $outpath -j $subject -m $masterpath -v $ver -p $partition_name -g $hg -z $gpu_partition -c $cont -x $cutoff -t SVA -f SVA" | sbatch -J SVA9_10 $slurm_sc --dependency=afterok:$jid8c| awk '{print $4}')
fi

