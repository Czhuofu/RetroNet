#! /bin/bash
usage="$(basename "$0") [-h] [-o j m v p i b c g] -- Discovering somatic MEI insertions supported with >=2 reads.

where:
    -h  show this help text
    -o  output folder path
    -j  subject ID
    -m  masterpath
    -v  version control for RetroSom (default 3)
    -p  slurm partition name
    -i  input type (1=sort.bam; 2=CleanBAM_and_ready_to_RetroDiscover)
    -b  BAM file (Input the path of sort.bam)
    -c  conda path (e.g. ~/miniconda3 or ~/Anaconda3)
    -g  reference genome (default hg38, supporting hg38, hg19 and b37)"

ver=3
hg=hg38
maxread=100
while getopts ":ho:j:m:v:p:i:b:c:g:n:" opt; do
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
    i) datatype="$OPTARG"
       ;;
    b) bam=$(readlink -f $OPTARG)
       ;;
    c) conda_path=$(readlink -f $OPTARG)
       ;;
    g) hg="$OPTARG"
       ;;
    n) maxreads="$OPTARG"
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

######################################
### Step0: Creating output folders ###
######################################
mkdir $outpath/$sub
mkdir $outpath/$sub/reads
mkdir $outpath/$sub/align
mkdir $outpath/$sub/QC
mkdir $outpath/$sub/temp
mkdir $outpath/$sub/script
mkdir $outpath/$sub/retro_v$ver
cd $outpath/$sub/script
cp $masterpath/pipeline/*sh $outpath/$sub/script

### updating the location of the reference sequences ###
echo -e "Alu\t$masterpath/refTE/sequence/ALU.fa" > $masterpath/refTE/TE_ALHS.bed
echo -e "L1\t$masterpath/refTE/sequence/L1HS.fa" >> $masterpath/refTE/TE_ALHS.bed
echo -e "HERVK\t$masterpath/refTE/sequence/HERVK.fa" >> $masterpath/refTE/TE_ALHS.bed
echo -e "HERVH\t$masterpath/refTE/sequence/HERVH.fa" >> $masterpath/refTE/TE_ALHS.bed
echo -e "SVA\t$masterpath/refTE/sequence/SVA.fa" >> $masterpath/refTE/TE_ALHS.bed

slurm_sc="-o %x.%A.output -e %x.%A.output -p $partition_name --mem=48gb --time=100:00:00"
slurm_sc_parallel="-o %x.%A.output -e %x.%A.output -p $partition_name --ntasks=1 --cpus-per-task=4 --mem=64gb --time=100:00:00"
slurm_mc="-o %x.%A.output -e %x.%A.output -p $partition_name --time=120:00:00 --ntasks=1 --cpus-per-task=10 --mem-per-cpu=16gb"
tmppath=$outpath/$sub/temp

##################################
### Step1: CleanBAM and QC_BAM ###
##################################
if [ "$datatype" == 1 ]
then
   jid1=$(echo -e '#!/bin/sh\n'"singularity exec -B $outpath,$tmppath,$masterpath,$bam $masterpath/pipeline/RetroNet.sif ./1_CleanBAM_QCBAM.sh -o $outpath -j $sub -b $bam" | sbatch -J CleanBAM_QCBAM $slurm_sc | awk '{print $4}')
fi

##################################################
### Step2: Discover candidate supporting reads ###
##################################################
if [ "$datatype" == 2 ]
then
   jid2=$(echo -e '#!/bin/sh\n'"singularity exec -B $outpath,$tmppath,$masterpath,$bam $masterpath/pipeline/RetroNet.sif ./2_RetroDiscover.sh -o $outpath -j $sub -m $masterpath -v $ver" | sbatch -J RetroDiscover $slurm_mc | awk '{print $4}')
else
   jid2=$(echo -e '#!/bin/sh\n'"singularity exec -B $outpath,$tmppath,$masterpath,$bam $masterpath/pipeline/RetroNet.sif ./2_RetroDiscover.sh -o $outpath -j $sub -m $masterpath -v $ver" | sbatch -J RetroDiscover $slurm_mc --dependency=afterok:$jid1 | awk '{print $4}')
fi

########################################################
### Step3: > Putative MEIs                           ###
###        > Remapping L1HS or AluY specific Alleles ###
###        > Matricies for L1/Alu supporting reads   ###
###        > Level1 prediction with RF, NB and LR    ###
########################################################
jid3=$(echo -e '#!/bin/sh\n'"singularity exec -B $outpath,$tmppath,$masterpath $masterpath/pipeline/RetroNet.sif ./3_step3to6.sh -o $outpath -j $sub -m $masterpath -v $ver -g $hg" | sbatch -J Step3to6 $slurm_sc_parallel --dependency=afterok:$jid2 | awk '{print $4}')
