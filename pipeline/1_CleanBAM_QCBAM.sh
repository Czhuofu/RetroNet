#!/bin/bash
usage="$(basename "$0") [-h] [-o j b] -- Remove PCR duplicates, secondary alignment and supplementary alignment. Check the alignment statistics, and average insert size.

where:
    -h  show this help text
    -o  output folder path
    -j  subject ID
    -b  BAM file (Input the path of sort.bam)"

while getopts ":ho:j:b:" opt; do
  case $opt in
    h) echo "$usage"
       exit
       ;;
    o) outpath=$(readlink -f $OPTARG)
       ;;
    j) sub="$OPTARG"
       ;;
    b) bam=$(readlink -f $OPTARG)
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

echo $sub
# Activate RetroNet virtual env #
source /opt/miniconda3/bin/activate RetroSom

#### CleanBAM ####
### remove PCR duplicates, secondary alignment and supplementary alignment ###
java -Xmx10G -jar /opt/miniconda3/envs/RetroSom/share/picard-2.13-1/picard.jar MarkDuplicates \
   I=$bam \
   O=/dev/stdout \
   ASSUME_SORTED=true \
   TMP_DIR=$outpath/$sub/align/ \
   MAX_RECORDS_IN_RAM=2000000 \
   REMOVE_DUPLICATES=true \
   VALIDATION_STRINGENCY=SILENT \
   METRICS_FILE=$outpath/$sub/QC/$sub.dedupp.metrics | \
   samtools view -bh -F 0x0100 -F 0x400 -F 0x800 - \
   > $outpath/$sub/align/$sub.final.bam
### index the BAM file ###
samtools index \
   $outpath/$sub/align/$sub.final.bam

#### QC_BAM ####
### check the alignment statistics ###
samtools flagstat \
    $outpath/$sub/align/$sub.final.bam \
    > $outpath/$sub/QC/$sub.flagstat
### check the average insert size ###
java -Xmx10G -jar /opt/miniconda3/envs/RetroSom/share/picard-2.13-1/picard.jar CollectInsertSizeMetrics \
   I=$outpath/$sub/align/$sub.final.bam \
   O=$outpath/$sub/QC/$sub.insert.txt \
   H=$outpath/$sub/QC/$sub.histo.pdf \
   M=0.001
sed -n 8,8p $outpath/$sub/QC/$sub.insert.txt | awk '{print $1"\t"$2}' > $outpath/$sub/insert.$sub.txt
