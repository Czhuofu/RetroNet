# A brief introduction for RetroSomV3
## Installation instructions
Basically, all the files were storage in the GITHUB.
Due to GitHub's file size restrictions, some large files/folder have been placed in Google Drive. You should download them from the link below and place them in the corresponding location.
For example, set the directory you download the RetroSomV2 as $yourdownloadpath

1. [RetroNet.sif](https://drive.google.com/file/d/1OwZvYDdvbTUrQXMidWQkYp9vbx8YjFrQ/view?usp=sharing) : Move to $yourdownloadpath/pipeline/
2. [hg38_100bp.bedGraph](https://drive.google.com/file/d/1IhiktWmqZSTcrPg2p9OIb5vtcQ5GLLjh/view?usp=sharing) : Move to $yourdownloadpath/RetroVis/
3. [position](https://drive.google.com/drive/folders/1L-XxCCGRMnNShd7ysbeM2kFIxkQANI9D?usp=sharing) : This is a folder, please download and move it to $yourdownloadpath/refTE/
4. [02_PE_level1](https://drive.google.com/drive/folders/197ogIPePEDBNah-Ff1IjNSKq7F3SMGzr?usp=sharing) : This is a folder, please download and move it to $yourdownloadpath/LINE/
5. [02_PE_level1](https://drive.google.com/drive/folders/18kA4IrlP7OKStuReX4koZ8dwx4sbjqzS?usp=sharing) : This is a folder, please download and move it to $yourdownloadpath/ALU/

## How to use
### Workflow


### Usage
This pipeline is divided into two parts, Singularity_Slurm_RetroSomV2.6.sh, including steps 1-6 in the workflow above; and Singularity_Slurm-RetroSom.step2.sh, which includes Step7-8 in the workflow above.
The input file for this pipeline is a Bam file. If the size of the Bam file you input is greater than 200G, please split it and name it with a suffix of - and a number.

For example: bigbam.bam &rarr; bigbam-1.bam bigbam-2.bam bigbam-3.bam ... (start from 1)
Recommend using samtools split for split. If you only have one Bam file, please also add the -1 suffix to the bam file.

### Analyze control first 

```
$yourdownloadpath/Singularity_Slurm_RetroSomV2.6.sh \

   -o /directory_path_for_output \

   -i ControlID \

   -m $yourdownloadpath \

   -r 1 \

   -g b37 (default hg38, supporting hg38, hg19 and b37)\

   -t 3 \

   -c /ControlID.bam \

   -n 150 (maximum number of supporting reads to be considered as a putative soamtic insertion) \

   -p 0.1 p_value cutoff (default p<0.1) \

   -e ControlID
```

**Then run the step2 to merge the result of Control:**

```
$yourdownloadpath/Singularity_Slurm_RetroSom.step2.sh \
  
   -o /directory_path_for_output \
  
   -i ControlID \
  
   -m $yourdownloadpath \
  
   -r 1 \
  
   -g b37 (default hg38, supporting hg38, hg19 and b37) \
  
   -t 0 \
  
   -c /ControlID.bam \
  
   -n 150 (maximum number of supporting reads to be considered as a putative soamtic insertion) \
  
   -p 0.1 p_value cutoff (default p<0.1) \
   
   -e ControlID \
  
   -l 1 (number of bam you split)
```

### Analyze case next

**Please ensure that the case and control use the same output path.**

```
$yourdownloadpath/Singularity_Slurm_RetroSomV2.6.sh \

   -o /directory_path_for_output \

   -i CaseID \

   -m $yourdownloadpath \

   -r 1 \

   -g b37 (default hg38, supporting hg38, hg19 and b37)\

   -t 3 \

   -c /CaseID.bam \

   -n 150 (maximum number of supporting reads to be considered as a putative soamtic insertion) \

   -p 0.1 p_value cutoff (default p<0.1) \

   -e ControlID_NoModel
```

**Then run the step2 to merge the result of Case:**

```
$yourdownloadpath/Singularity_Slurm_RetroSom.step2.sh \
  
   -o /directory_path_for_output \
  
   -i CaseID \
  
   -m $yourdownloadpath \
  
   -r 1 \
  
   -g b37 (default hg38, supporting hg38, hg19 and b37) \
  
   -t 0 \
  
   -c /CaseID.bam \
  
   -n 150 (maximum number of supporting reads to be considered as a putative soamtic insertion) \
  
   -p 0.1 p_value cutoff (default p<0.1) \
   
   -e ControlID_NoModel \
  
   -l 1 (number of bam you split)
```

A folder called CaseID_NoModel will be generated, result svg will be in CaseID_NoModel/visual/
