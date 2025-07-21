# Germline MEI Detection with RetroNet (Demo: HG00096)
This repository provides a step-by-step demonstration of using RetroNet, originally developed for somatic mobile element insertion (MEI) detection, to identify germline MEIs from Illumina short-read whole-genome sequencing (WGS) data.

Here, we showcase how to apply RetroNet on a representative sample, HG00096, from the HGSVC (Human Genome Structural Variation Consortium) dataset. Although RetroNet was primarily designed for somatic MEI detection, it can be extended to call non-reference germline MEIs by leveraging an unrelated high-coverage Illumina sample (e.g., NA12878) as a control to suppress common background signals.

## System Requirements
**Cluster Management**:
- `RetroNet` requires a cluster environment managed by `SLURM`.
- Memory Requirements:
   - Single-core tasks: Require 48 GB.
   - Multi-core tasks: Require 16 GB per core.

**Containerization**:
- `Singularity` is required to execute the containerized environment.
   - `RetroNet.sif` was created using `Singularity v3.7.3`, 
   and has been tested on `Singularity v3.8.6`, `v1.2.2-1.el8`, `v3.8.7-1.el7`

## Installation Guide
**Install from Github**
````bash
git clone https://github.com/Czhuofu/RetroNet.git
````
Due to GitHub's file size limitations, some large files and folders are hosted on Google Drive. 
Please download them using the link provided below and place them in the appropriate locations. 
For instance, after clone the `RetroNet` files, set the directory path as `$yourdownloadpath`

- [RetroNet.sif](https://drive.google.com/file/d/1OUg-L2sQ7ucaNsTXolIus0uCT5VG29CC/view?usp=drive_link) : `mv RetroNet.sif $yourdownloadpath/RetroNet/pipeline/`
- [hg38_100bp.bedGraph](https://drive.google.com/file/d/1IhiktWmqZSTcrPg2p9OIb5vtcQ5GLLjh/view?usp=sharing) : `mv hg38_100bp.bedGraph $yourdownloadpath/RetroNet/`
- [position](https://drive.google.com/drive/folders/1L-XxCCGRMnNShd7ysbeM2kFIxkQANI9D?usp=sharing) : This is a folder `mv ./position $yourdownloadpath/RetroNet/refTE/`
- [b37_100bp.bedGraph](https://drive.google.com/file/d/14eOmzhz0pMYpfuLU5spLgwuZJv8_n75R/view?usp=drive_link) **(Optional)** : If you need to analyse b37 bam files `mv b37_100bp.bedGraph $yourdownloadpath/RetroNet/`

```bash
# Download NA12878 filter as control
wget https://github.com/Czhuofu/RetroNet/raw/Germline-MEI-detection-with-RetroNet/NA12878_100X_hg38_filter2.tar.gz
tar -zxvf NA12878_100X_hg38_filter2.tar.gz
mv NA12878_100X_hg38_filter2 $yourdownloadpath

# Replace a modified script
wget https://github.com/Czhuofu/RetroNet/raw/Germline-MEI-detection-with-RetroNet/generate_bed.py
mv generate_bed.py $yourdownloadpath/RetroNet/pipeline/
```

## Usage (HG00096 Demo)
### Download HG00096 from HGSVC
`0_Download_and_Preprocess_HG00096.sh`
```bash
#!/bin/bash

# 0) Path to Workfolder (Replace with your downlaod path)
Workfolder=$yourdownloadpath

# 1) Downlad HG00096 WGS from HGSVC
cd ${Workfolder}
wget ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR324/ERR3240114/HG00096.final.cram
# md5sum d3354f61a055adfcfc988470bc507b2d

# 2) Convert CRAM to BAM file
mkdir HG00096
samtools view -bh -F 0x0100 -F 0x400 -F 0x800 \
    -T ./GRCh38DH/GRCh38_full_analysis_set_plus_decoy_hla.fa \
    HG00096.final.cram \
    -o ${Workfolder}/HG00096/HG00096.clean.bam
samtools index ${Workfolder}/HG00096/HG00096.clean.bam
```

### Run RetroNet step 1
`1_Run_RetroNet_step1_for_HG00096.sh`
```bash
#!/bin/bash

# 0) Path to Workfolder
Workfolder=$yourdownloadpath

# 1) Run Step 1 for HG00096
sub=HG00096
bash ${Workfolder}/RetroNet/Singularity_Slurm_RetroNet_step1.sh \
   -o ${Workfolder} \
   -j ${sub} \
   -m ${Workfolder}/RetroNet \
   -g hg38 \
   -p tiny \
   -i 1 \
   -b ${Workfolder}/HG00096/HG00096.clean.bam \
   -n 100
```

### Run RetroNet step 2
`2_Run_RetroNet_step2_by_using_NA12878_as_control.sh`
```bash
#!/bin/bash

# 0) Path to Workfolder
Workfolder=$yourdownloadpath

# 1) Run Step 2 for HG00096 by using NA12878 as control
sub=HG00096
bash ${Workfolder}/RetroNet/Singularity_Slurm_RetroNet_step2.sh \
   -o ${Workfolder} \
   -j ${sub} \
   -m ${Workfolder}/RetroNet \
   -g hg38 \
   -p tiny \
   -c ${Workfolder}/RetroNet/NA12878_100X_hg38_filter2 \
   -x 0.95 \
   -l 1
```

## Expected Output 
- Detailed somatic MEIs infomation `sub=HG00096`
    - `$yourdownloadpath/$sub/retro_v3/$sub.LINE.bed`
    - `$yourdownloadpath/$sub/retro_v3/$sub.ALU.bed`
    - `$yourdownloadpath/$sub/retro_v3/$sub.SVA.bed`