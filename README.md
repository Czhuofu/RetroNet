
# Overview
- [RetroNet](#1)
- [System Requirements](#2)
- [Installation Guide](#3)
- [Usage](#4)
- [Demo](#5)
- [Analyze HG008 Cancer Cell](#6)
- [Extended Data and Files](#7)
<p id="1"></p>

# RetroNet
`RetroNet` is a computational tool that detects somatic mobile element insertions (MEIs) in human genomes using deep learning. 
By encoding sequencing reads into images, it identifies `L1`, `Alu`, and `SVA` insertions with high precision, even at low frequencies.
### Workflow
![Retro3 2](https://github.com/user-attachments/assets/10da5f78-0252-4d01-8145-b3f5c36cfca3)
<p id="2"></p>

# System Requirements
**Cluster Management**:
- `RetroNet` requires a cluster environment managed by `SLURM`.
- Memory Requirements:
   - Single-core tasks: Require 48 GB.
   - Multi-core tasks: Require 16 GB per core.

**Containerization**:
- `Singularity` is required to execute the containerized environment.
   - `RetroNet.sif` was created using `Singularity v3.7.3`, 
   and has been tested on `Singularity v3.8.6`, `v1.2.2-1.el8`, `v3.8.7-1.el7`
<br/>
<br/>
<p id="3"></p>

# Installation Guide
**Install from Github**
````bash
git clone https://github.com/Czhuofu/RetroNet.git
````
Due to GitHub's file size limitations, some large files and folders are hosted on Google Drive. 
Please download them using the link provided below and place them in the appropriate locations. 
For instance, after clone the `RetroNet` files, set the directory path as `$yourdownloadpath`

- [RetroNet.sif](https://drive.google.com/file/d/1OUg-L2sQ7ucaNsTXolIus0uCT5VG29CC/view?usp=drive_link) : `mv RetroNet.sif $yourdownloadpath/pipeline/`
- [hg38_100bp.bedGraph](https://drive.google.com/file/d/1IhiktWmqZSTcrPg2p9OIb5vtcQ5GLLjh/view?usp=sharing) : `mv hg38_100bp.bedGraph $yourdownloadpath/RetroNet/`
- [position](https://drive.google.com/drive/folders/1L-XxCCGRMnNShd7ysbeM2kFIxkQANI9D?usp=sharing) : This is a folder `mv -r ./position $yourdownloadpath/refTE/`
- [b37_100bp.bedGraph](https://drive.google.com/file/d/14eOmzhz0pMYpfuLU5spLgwuZJv8_n75R/view?usp=drive_link) **(Optional)** : If you need to analyse b37 bam files `mv b37_100bp.bedGraph $yourdownloadpath/RetroNet/`
<br/>
<br/>
<p id="4"></p>

# Usage
**This pipeline consists of two main steps:**
1. Analyze control sample
   - Run [`Singularity_Slurm_RetroNet_step1.sh`](#1-analyze-control-sample) screens candidate MEI sites and extracts supporting reads from the input `BAM` file.
   - Run [`Singularity_Slurm_RetroNet_generate_control.sh`](#1-analyze-control-sample) construct germline MEIs in control sample.
2.  Analyze case sample 
    - Run [`Singularity_Slurm_RetroNet_step1.sh`](#2-analyze-case-sample) screens candidate MEI sites and extracts supporting reads from the input `BAM` file.
    - Run [`Singularity_Slurm_RetroNet_step2.sh`](#2-analyze-case-sample) excludes germline MEIs by comparing with control samples. 
It also encodes the supporting reads into `PNG` images for deep learning model predictions and human observation.

**Note 1:**

If the input control or tissue `BAM` files exceed 100X sequence depth (approximately larger than 200 GB), please split them into smaller parts and append a suffix in the format `-number` to the filenames.
- Example: `bigbam.bam` → `bigbam-1.bam`, `bigbam-2.bam`, `bigbam-3.bam`, ... (starting from `-1`).
````bash
# If your BAM file consists multiple DNA readgroups, you can follow these to split your bam files
subjectID=bigbam
samtools split -f "%*-%#.bam" ${subjectID}.bam
# Rename the splited BAM -> bigbam-1.bam, bigbam-2.bam, bigbam-3.bam...
for file in ${subjectID}-*.bam; do
  num=$(echo "$file" | grep -oP '(?<=-)\d+(?=\.bam)')
  new_num=$((num + 1))
  new_name=${subjectID}-${new_num}.bam
  mv $file $new_name
done
````
**Note 2:**

Currently, the pipeline supports both hg38 and b37 genome builds.
<br/>
## Detail code for running pipeline
### 1. Analyze Control Sample
````bash
### if you have multiple bam files, please run step 1 for each split bam
$yourdownloadpath/Singularity_Slurm_RetroNet_step1.sh \
   -o $outpath \              # Output path for your analysis
   -j $subject_control \      # Subject ID of contol sample
   -m $yourdownloadpath \     # The path where you place RetroNet
   -g hg38 \                  # Default hg38, supporting hg38 and b37
   -p $slurm_cpu_partition \  # SLURM cpu partition
   -i $input_type \           # 1=sort.bam; 2=Cleaned BAM (excludes duplicates, supplementary alignment, secondary alignment)
   -b $BAM_file \             # BAM file of the subject
   -n 100                     # Maximum number of supporting reads analysis at a time, default 100
````
````bash
### run after finish every step 1 task
$yourdownloadpath/Singularity_Slurm_RetroNet_generate_control.sh \
   -o $outpath \              # Output path for your analysis
   -j $subject_control \      # Subject ID of contol sample
   -m $yourdownloadpath \     # The path where you place RetroNet
   -g hg38 \                  # Default hg38, supporting hg38 and b37
   -p $slurm_cpu_partition \  # SLURM cpu partition
   -l $num_control_BAM        # number of control BAM you split
````
### 2. Analyze Case Sample
**Noted**: Please ensure use same output path for your `control` and `case` sample
````bash
$yourdownloadpath/Singularity_Slurm_RetroNet_step1.sh \
   -o $outpath \              # Output path for your analysis (same as your control sample)
   -j $subject_case \         # Subject ID of case sample
   -m $yourdownloadpath \     # The path where you place RetroNet
   -g hg38 \                  # Default hg38, supporting hg38 and b37
   -p $slurm_cpu_partition \  # SLURM cpu partition
   -i input_type \            # 1=sort.bam; 2=Cleaned BAM (excludes duplicates, supplementary alignment, secondary alignment)
   -b $BAM_file \             # BAM file of the subject
   -n 100                     # Maximum number of supporting reads analysis at a time, default 100
````
**Noted**: Please ensure you have finished all of `Singularity_Slurm_RetroNet_step1.sh` and `Singularity_Slurm_RetroNet_generate_control.sh`
```bash
### If you split multiple BAM for control
### Your control will be subject_control=${subject_control}_Combined
$yourdownloadpath/Singularity_Slurm_RetroNet_step2.sh \
   -o $outpath \              # Output path for your analysis (same as your control sample)
   -j $subject_case \         # Subject ID of case sample
   -m $yourdownloadpath \     # The path where you place RetroNet
   -g hg38 \                  # Default hg38, supporting hg38 and b37
   -p $slurm_cpu_partition \  # SLURM cpu partition
   -c $subject_control \      # If you have multiple BAM file of control, please use
   -z $slurm_gpu_partition \  # (Optional parameter) SLURM gpu partition if you have, default N
   -x 0.95 \                  # Probability cutoff default 0.95  
   -l $num_case_BAM           # number of case BAM you split
```
## Expected Output 

- If you do no split `BAM` file of case sample:
   - Detailed somatic MEIs infomation
      - `$outpath/$subject_case/retro_v3/${subject_case}.LINE.bed`
      - `$outpath/$subject_case/retro_v3/${subject_case}.ALU.bed`
      - `$outpath/$subject_case/retro_v3/${subject_case}.SVA.bed`
   - Prediction of `PNG` images pass probability `cutoff`
      - `$outpath/$subject_case/RetroNet/LINE_Inspected_${subject_case}_cut0.95.txt`
      - `$outpath/$subject_case/RetroNet/ALU_Inspected_${subject_case}_cut0.95.txt`
      - `$outpath/$subject_case/RetroNet/SVA_Inspected_${subject_case}_cut0.95.txt`
   - Plot `svg` images of somatic MEIs using `RetroVis`
      - `$outpath/$subject_case/visual/`: somatic MEIs of `L1`, `Alu`, and `SVA`
      
- If you split your case `BAM` file, the output folder will be `subject_case=${subject_case}_Combined`
<p id="5"></p>

# Demo
This demo has been tested under the `Burgundy HPC` in City University of Hong Kong.
Please download the demo output, and input `BAM` files use this link: [Demo.tar.gz](https://drive.google.com/file/d/1oErysrV78xyxqi_0H9PtMRRjVWGLXVkN/view?usp=sharing)
- Input `BAM` files: `Demo/normal.bam` and `Demo/tumour.bam`
- Output folder: `Demo/normal` and `Demo/tumour`
### 1. Analyze the normal tissue
````bash
# Initialize the input of normal tissue
outpath=/gpfs1/scratch/zhinanlin2/Demo_GitHub_RetroNet/Demo
subject_control=normal
yourdownloadpath=/gpfs1/scratch/zhinanlin2/Demo_GitHub_RetroNet/RetroNet
slurm_cpu_partition=tiny
input_type=1
BAM_file=/gpfs1/scratch/zhinanlin2/Demo_GitHub_RetroNet/Demo/normal.bam
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
**Consumed time for normal tissue:**
- `Singularity_Slurm_RetroNet_step1.sh` 6m12s
- `Singularity_Slurm_RetroNet_generate_control.sh` 4s
### 2. Analyze the tumour cell line
````bash
# Initialize the input of tumor sample
outpath=/gpfs1/scratch/zhinanlin2/Demo_GitHub_RetroNet/Demo
subject_case=tumour
yourdownloadpath=/gpfs1/scratch/zhinanlin2/Demo_GitHub_RetroNet/RetroNet
slurm_cpu_partition=tiny
slurm_gpu_partition=N
input_type=1
BAM_file=/gpfs1/scratch/zhinanlin2/Demo_GitHub_RetroNet/Demo/tumour.bam
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
**Consumed time for tumour cell:**
- `Singularity_Slurm_RetroNet_step1.sh` 6m28s
- `Singularity_Slurm_RetroNet_step2.sh` 2h47m30s
<p id="6"></p>

# Example: Analysis of HG008 Cancer Cell Line
Please see the [HG008_analyze_log.md](https://github.com/Czhuofu/RetroNet/blob/master/HG008_analyze_log.md)
<p id="7"></p>

# Extended Data and Supplementary Materials

This part provides extended data files related to our manuscript, include manual validations, visualizations, and raw calls of other tools (PALMER and xTea long). All files are publicly accessible and are intended to enhance transparency and reproducibility.

### Manual check of Germline MEI discordant with PacBio xTea long

We manually validated germline MEIs from 11 Illumina DNA WGS that were not recovered by **Pangenome xTea long** calls from PacBio HiFi data:

- 📊 **Validation of labels:**
<div align="center">
<img src="https://raw.githubusercontent.com/RLin1998/user_file/2d6a1bbd39739ec25da8d9393547262019433d6a/summary_MEI_short_and_Pacbio.png" width="400">
</div>

- **Extended File 1**: Manual check of LINE-1 germline insertions not identified by PacBio xTea_long [👉 View File](https://drive.google.com/file/d/1NPqNDS65UXYq3qETa2U8fkxAVQgdz3Ct/view?usp=sharing)

- **Extended File 2**: Manual check of *Alu* germline insertions not identified by PacBio xTea_long [👉 View File](https://drive.google.com/file/d/1OTmkHFIbf8eypNMO7uLULvT6kcUfZdP1/view?usp=sharing)

- **Extended File 3**: Manual check of SVA germline insertions not identified by PacBio xTea_long or the polymorphic SVA dataset [👉 View File](https://drive.google.com/file/d/1pzh1s4_Ee_ck-NgVgp5gJD-tiBVSwVR8/view?usp=sharing)

---
### Manual check for benchmarking in HG008-T and Patient DTB-205 cfDNA

- **Extended File 4**: Detailed visualization of **tumor somatic L1 insertions** in **HG008-T** [👉 View File](https://drive.google.com/file/d/1-CZgrNEfO88l8kr6kF74KE0ot83ei-XQ/view?usp=sharing)

- **Extended File 5**: Summary of **tumor somatic L1 insertions** in HG008‑T and their detection status by **PALMER** or **xTea_long** [👉 View File](https://drive.google.com/file/d/1sxCzCEjxUk2F8b4kjsOSSGoGuTdDb9V2/view?usp=drive_link)

- **Extended File 6**: **Tumor somatic *Alu* insertions** in HG008‑T detected by PALMER or xTea_long [👉 View File](https://drive.google.com/file/d/1GAEFjkO0ayvMEMAh0TN6hHauUmTxkuwu/view?usp=sharing)

- **Extended File 7**: **Tumor somatic SVA insertions** in HG008‑T detected by PALMER or xTea_long [👉 View File](https://drive.google.com/file/d/1dX8Boc_6eiWR6AxSuw2aXKYEZxpz4J5q/view?usp=sharing)

- **Extended File 8**: Detailed visualization of patient DTB-205 **cfDNA somatic L1 insertions** identified by RetroNet with low mosaicism in the tumor [👉 View File](https://drive.google.com/file/d/1YJbWhk1fXZkVHtlTEFzjQxUymS4p9q4b/view?usp=sharing)

---

### 📂 Raw MEI calls from PALMER and xTea:

- **LINE-1, *Alu*, and SVA insertions** in **HG008-T** and **HG008**, detected by **xTea_long** and **PALMER** [👉 View Folder](https://drive.google.com/drive/folders/1LJL7DCN2-GzOGHCTqxASpAxEyBnStO61?usp=sharing)

# Licence
RetroNet is licensed under the MIT License, approved by the Open Source Initiative (OSI). This license allows users to use, copy, modify, and distribute the software, provided proper attribution is given to the authors. The full license text is available in the [LICENSE](https://github.com/Czhuofu/RetroNet/blob/master/LICENSE).
