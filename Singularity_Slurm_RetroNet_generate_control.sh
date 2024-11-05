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
    -l  number of the libraries (default is 1)"

readgroup=1
ver=3
hg=hg38

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

###############################################################
### Step7: Pairing the supporting reads and getting control ###
###############################################################
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
    mkdir $outpath/$sub/retro_v$ver\_0
    mkdir $outpath/$sub/retro_v$ver\_1
    mkdir $outpath/$sub/retro_v$ver\_0/LINE
    mkdir $outpath/$sub/retro_v$ver\_1/LINE
    mkdir $outpath/$sub/retro_v$ver\_0/ALU
    mkdir $outpath/$sub/retro_v$ver\_1/ALU
    mkdir $outpath/$sub/retro_v$ver\_0/SVA
    mkdir $outpath/$sub/retro_v$ver\_1/SVA
    cd $outpath/$sub/script
    cp $masterpath/RetroNet/pipeline/* $outpath/$sub/script
    # submit jobs #
    jid7a=$(echo -e '#!/bin/sh\n'"singularity exec -B $outpath,$tmppath,$masterpath $masterpath/pipeline/RetroNet.sif ./7_CombineLib.sh -o $outpath -j $subject -m $masterpath -v $ver -s 0 -l $readgroup -g $hg" | sbatch -J CombineStrand0 $slurm_sc | awk '{print $4}')
    jid7b=$(echo -e '#!/bin/sh\n'"singularity exec -B $outpath,$tmppath,$masterpath $masterpath/pipeline/RetroNet.sif ./7_CombineLib.sh -o $outpath -j $subject -m $masterpath -v $ver -s 1 -l $readgroup -g $hg" | sbatch -J CombineStrand1 $slurm_sc | awk '{print $4}')
    jid7c=$(echo -e '#!/bin/sh\n'"singularity exec -B $outpath,$tmppath,$masterpath $masterpath/pipeline/RetroNet.sif ./7a_Merge_Anchor.sh -o $outpath -j $subject -m $masterpath -v $ver -l $readgroup" | sbatch -J CombineAnchor $slurm_sc --dependency=afterok:$jid7a:$jid7b | awk '{print $4}')
    jid7d=$(echo -e '#!/bin/sh\n'"singularity exec -B $outpath,$tmppath,$masterpath $masterpath/pipeline/RetroNet.sif ./7b_generate_control.sh -o $outpath -j $subject -m $masterpath -v $ver -l $readgroup" | sbatch -J Generatecontrol $slurm_sc --dependency=afterok:$jid7c | awk '{print $4}')
else
    cd $outpath/$sub/script
    cp $masterpath/RetroNet/pipeline/* $outpath/$sub/script
    jid7d=$(echo -e '#!/bin/sh\n'"singularity exec -B $outpath,$tmppath,$masterpath $masterpath/pipeline/RetroNet.sif ./7b_generate_control.sh -o $outpath -j $subject -m $masterpath -v $ver -l $readgroup" | sbatch -J Generatecontrol $slurm_sc | awk '{print $4}')
fi

