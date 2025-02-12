# Overview
- [Overview](#overview)
- [ILMN-PCR-free-2 (I2) DNA WGS Data - Sequenced by BCM](#ilmn-pcr-free-2-i2-dna-wgs-data---sequenced-by-bcm)
  - [Download BAM Files](#download-bam-files)
  - [Analyze I2 Control Sample - HG008-N-D](#analyze-i2-control-sample---hg008-n-d)
  - [Analyze I2 Tumor Sample - HG008-T](#analyze-i2-tumor-sample---hg008-t)
- [ILMN-PCR-free-3 (I3) DNA WGS Data - Sequenced by NYGC](#ilmn-pcr-free-3-i3-dna-wgs-data---sequenced-by-nygc)
  - [Download BAM Files](#download-bam-files-1)
  - [Analyze I3 Control Sample - HG008-N-D](#analyze-i3-control-sample---hg008-n-d)
  - [Analyze I3 Tumor Sample - HG008-T](#analyze-i3-tumor-sample---hg008-t)
  - [(Split BAM) Analyze I3 Tumor Sample - HG008-T](#split-bam-analyze-i3-tumor-sample---hg008-t)
<p id="1"></p>

# ILMN-PCR-free-2 (I2) DNA WGS Data - Sequenced by BCM
<p id="2"></p>

## Download BAM Files
````bash
# Download HG008-N-D (Duodenal Solid Tissue Normal)
wget https://42basepairs.com/download/web/giab/data_somatic/HG008/Liss_lab/BCM_Illumina-WGS_20240313/HG008-N-D_Illumina_169x_GRCh38-GIABv3.bam
wget https://42basepairs.com/download/web/giab/data_somatic/HG008/Liss_lab/BCM_Illumina-WGS_20240313/HG008-N-D_Illumina_169x_GRCh38-GIABv3.bam.bai
# Download HG008-T (Tumor cells)
wget https://42basepairs.com/download/web/giab/data_somatic/HG008/Liss_lab/BCM_Illumina-WGS_20240313/HG008-T_Illumina_195x_GRCh38-GIABv3.bam
wget https://42basepairs.com/download/web/giab/data_somatic/HG008/Liss_lab/BCM_Illumina-WGS_20240313/HG008-T_Illumina_195x_GRCh38-GIABv3.bam.bai
````
<p id="3"></p>

## Analyze I2 Control Sample - HG008-N-D
````bash
# Initialization
outpath=/home/zhinanlin/scratch/HG008_BCM
subject_control=HG008-N-D
yourdownloadpath=/home/zhinanlin/scratch/RetroNet
slurm_cpu_partition=batch
input_type=1
BAM_file=/home/zhinanlin/scratch/HG008_BCM/HG008-N-D_Illumina_169x_GRCh38-GIABv3.bam
num_control_BAM=1
# Run step1 for control
$yourdownloadpath/Singularity_Slurm_RetroNet_step1.sh \
-o $outpath \
-j $subject_control \
-m $yourdownloadpath \
-g hg38 \
-p $slurm_cpu_partition \
-i $input_type \
-b $BAM_file \
-n 100
# Wait for step1 finish...
$yourdownloadpath/Singularity_Slurm_RetroNet_generate_control.sh \
-o $outpath \
-j $subject_control \
-m $yourdownloadpath \
-g hg38 \
-p $slurm_cpu_partition \
-l $num_control_BAM
````
<p id="4"></p>

## Analyze I2 Tumor Sample - HG008-T
**Noted:** Please ensure you have finished all of `Singularity_Slurm_RetroNet_step1.sh` and `Singularity_Slurm_RetroNet_generate_control.sh`
````bash
# Initialization
outpath=/home/zhinanlin/scratch/HG008_BCM
subject_case=HG008-T
yourdownloadpath=/home/zhinanlin/scratch/RetroNet
slurm_cpu_partition=batch
slurm_gpu_partition=N
input_type=1
BAM_file=/home/zhinanlin/scratch/HG008_BCM/HG008-T_Illumina_195x_GRCh38-GIABv3.bam
num_case_BAM=1
# Run step1 for tumor
$yourdownloadpath/Singularity_Slurm_RetroNet_step1.sh \
-o $outpath \
-j $subject_case \
-m $yourdownloadpath \
-g hg38 \
-p $slurm_cpu_partition \
-i $input_type \
-b $BAM_file \
-n 100
# Wait for step1 finish...
$yourdownloadpath/Singularity_Slurm_RetroNet_step2.sh \
-o $outpath \
-j $subject_case \
-m $yourdownloadpath \
-g hg38 \
-p $slurm_cpu_partition \
-c $subject_control \
-z $slurm_gpu_partition \
-x 0.95 \
-l $num_case_BAM
````
<p id="5"></p>

# ILMN-PCR-free-3 (I3) DNA WGS Data - Sequenced by NYGC
<p id="6"></p>

## Download BAM Files
````bash
# Download HG008-N-D (Duodenal Solid Tissue Normal)
wget https://42basepairs.com/download/web/giab/data_somatic/HG008/Liss_lab/NYGC_Illumina-WGS_20231023/HG008-N-D_Illumina_118x_GRCh38-GIABv3.bam
wget https://42basepairs.com/download/web/giab/data_somatic/HG008/Liss_lab/NYGC_Illumina-WGS_20231023/HG008-N-D_Illumina_118x_GRCh38-GIABv3.bam.bai
# Download HG008-T (Tumor cells)
wget https://42basepairs.com/download/web/giab/data_somatic/HG008/Liss_lab/NYGC_Illumina-WGS_20231023/HG008-T_Illumina_161x_GRCh38-GIABv3.bam
wget https://42basepairs.com/download/web/giab/data_somatic/HG008/Liss_lab/NYGC_Illumina-WGS_20231023/HG008-T_Illumina_161x_GRCh38-GIABv3.bam.bai
````
<p id="7"></p>

## Analyze I3 Control Sample - HG008-N-D
````bash
# Initialization
outpath=/home/zhinanlin/scratch/HG008_NYGC
subject_control=HG008-N-D
yourdownloadpath=/home/zhinanlin/scratch/RetroNet
slurm_cpu_partition=batch
input_type=1
BAM_file=/home/zhinanlin/scratch/HG008_NYGC/HG008-N-D_Illumina_118x_GRCh38-GIABv3.bam
num_control_BAM=1
# Run step1 for control
$yourdownloadpath/Singularity_Slurm_RetroNet_step1.sh \
-o $outpath \
-j $subject_control \
-m $yourdownloadpath \
-g hg38 \
-p $slurm_cpu_partition \
-i $input_type \
-b $BAM_file \
-n 100
# Wait for step1 finish...
$yourdownloadpath/Singularity_Slurm_RetroNet_generate_control.sh \
-o $outpath \
-j $subject_control \
-m $yourdownloadpath \
-g hg38 \
-p $slurm_cpu_partition \
-l $num_control_BAM
````
<p id="8"></p>

## Analyze I3 Tumor Sample - HG008-T
````bash
# Initialization
outpath=/home/zhinanlin/scratch/HG008_NYGC
subject_case=HG008-T
yourdownloadpath=/home/zhinanlin/scratch/RetroNet
slurm_cpu_partition=batch
slurm_gpu_partition=N
input_type=1
BAM_file=/home/zhinanlin/scratch/HG008_NYGC/HG008-T_Illumina_161x_GRCh38-GIABv3.bam
num_case_BAM=1
# Run step1 for tumor
$yourdownloadpath/Singularity_Slurm_RetroNet_step1.sh \
-o $outpath \
-j $subject_case \
-m $yourdownloadpath \
-g hg38 \
-p $slurm_cpu_partition \
-i $input_type \
-b $BAM_file \
-n 100
# Wait for tumor step1 and control analysis finish...
$yourdownloadpath/Singularity_Slurm_RetroNet_step2.sh \
-o $outpath \
-j $subject_case \
-m $yourdownloadpath \
-g hg38 \
-p $slurm_cpu_partition \
-c $subject_control \
-z $slurm_gpu_partition \
-x 0.95 \
-l $num_case_BAM
````
<p id="9"></p>

## (Split BAM) Analyze I3 Tumor Sample - HG008-T
**Split HG008-T BAM by DNA libraries**
````bash
# If your BAM file consists multiple DNA libraries
subjectID=HG008-T
BAM_file=/home/zhinanlin/scratch/HG008_NYGC/HG008-T_Illumina_161x_GRCh38-GIABv3.bam
samtools split -f ${subjectID}-%#.bam $BAM_file
# Rename the splited BAM -> HG008-T-1.bam, HG008-T-2.bam, HG008-T-3.bam...
for file in ${subjectID}-*.bam; do
  num=$(echo "$file" | grep -oP '(?<=-)\d+(?=\.bam)')
  new_num=$((num + 1))
  new_name=${subjectID}-${new_num}.bam
  mv $file $new_name
done
````
**Analyze**
````bash
# Initialization
outpath=/home/zhinanlin/scratch/HG008_NYGC
yourdownloadpath=/home/zhinanlin/scratch/RetroNet
slurm_cpu_partition=batch
slurm_gpu_partition=N
input_type=1
# Get number of splited BAM and Run step1 for each BAM
subject_case=HG008-T
num_case_BAM=$(ls ${subject_case}-*.bam | wc -l)
count=1
for bam_split in ${subject_case}-*.bam; do
subject_case_split=${subject_case}-${count}
# Run step1 for each splited tumor BAM
$yourdownloadpath/Singularity_Slurm_RetroNet_step1.sh \
-o $outpath \
-j $subject_case_split \
-m $yourdownloadpath \
-g hg38 \
-p $slurm_cpu_partition \
-i $input_type \
-b $bam_split \
-n 100
count=$((count + 1))
done
# Wait for all tumor step1 and control analysis finish...
subject_control=${subject_control}_Combined
$yourdownloadpath/Singularity_Slurm_RetroNet_step2.sh \
-o $outpath \
-j $subject_case \
-m $yourdownloadpath \
-g hg38 \
-p $slurm_cpu_partition \
-c $subject_control \
-z $slurm_gpu_partition \
-x 0.95 \
-l $num_case_BAM
````

