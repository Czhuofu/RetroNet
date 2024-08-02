# A brief introduction for RetroSomV3
## Installation instructions
Basically, all the files were storage in the GITHUB.
Due to GitHub's file size restrictions, some large files/folder have been placed in Google Drive. You should download them from the link below and place them in the corresponding location.
For example, set the directory you download the RetroSomV3 as $yourdownloadpath.

1. [RetroNet.sif](https://drive.google.com/file/d/1OwZvYDdvbTUrQXMidWQkYp9vbx8YjFrQ/view?usp=sharing) : Move to $yourdownloadpath/pipeline/
2. [hg38_100bp.bedGraph](https://drive.google.com/file/d/1IhiktWmqZSTcrPg2p9OIb5vtcQ5GLLjh/view?usp=sharing) : Move to $yourdownloadpath/RetroNet/
3. [position](https://drive.google.com/drive/folders/1L-XxCCGRMnNShd7ysbeM2kFIxkQANI9D?usp=sharing) : This is a folder, please download and move it to $yourdownloadpath/refTE/
4. [02_PE_level1](https://drive.google.com/drive/folders/197ogIPePEDBNah-Ff1IjNSKq7F3SMGzr?usp=sharing) : This is a folder, please download and move it to $yourdownloadpath/LINE/
5. [02_PE_level1](https://drive.google.com/drive/folders/18kA4IrlP7OKStuReX4koZ8dwx4sbjqzS?usp=sharing) : This is a folder, please download and move it to $yourdownloadpath/ALU/
6. [b37_100bp.bedGraph](https://drive.google.com/file/d/14eOmzhz0pMYpfuLU5spLgwuZJv8_n75R/view?usp=drive_link) : If you need to analyse b37 bam files, please download and move to $yourdownloadpath/RetroNet/

## How to use
### Workflow

![Retro3 2](https://github.com/user-attachments/assets/10da5f78-0252-4d01-8145-b3f5c36cfca3)

### Usage
This pipeline is divided into two steps, you can run Singularity_Slurm_RetroNet_step1.sh and Singularity_Slurm_RetroNet_step2.sh to finish the analysis process.
The input file for this pipeline is a Bam file. If the size of the Bam file you input is greater than 200G, please split it and name it with a suffix of - and a number.

For example: bigbam.bam &rarr; bigbam-1.bam bigbam-2.bam bigbam-3.bam ... (start from 1)
Recommend using samtools split for split. If you only have one Bam file, please also add the -1 suffix to the bam file.
Now, the pipeline can support hg38 and b37.

### Analyze control first 

```
$yourdownloadpath/Singularity_Slurm_RetroNet_step1.sh \

   -o /directory_path_for_output \

   -j ControlID \

   -m $yourdownloadpath \

   -v 3 (version control for RetroSom, default 3) \

   -g hg38 (default hg38, supporting hg38, hg19 and b37)\

   -p your_slurm_partition_name \

   -i input_type (1=sort.bam; 2=CleanBAM_and_ready_to_RetroDiscover)

   -b /ControlID.bam \

   -n 100 (maximum number of supporting reads analysis at a time, default 100) 

```

**Then run the step2 to merge the result of Control:**

```
$yourdownloadpath/Singularity_Slurm_RetroNet_step2.sh \
  
   -o /directory_path_for_output \
  
   -j ControlID \
  
   -m $yourdownloadpath \
  
   -v 3 (version control for RetroSom, default 3) \
  
   -g hg38 (default hg38, supporting hg38, hg19 and b37) \
  
   -p your_slurm_partition_name \
  
   -c ControlID \
  
   -z "N" (gpu_partition name, default "N") \
  
   -x 0.99 ( probability cutoff default 0.99) \
  
   -l 1 (number of bam you split)
```

### Analyze case next

**Please ensure that the case and control use the same output path.**

```
$yourdownloadpath/Singularity_Slurm_RetroNet_step1.sh \

   -o /directory_path_for_output \

   -j CaseID \

   -m $yourdownloadpath \

   -v 3 (version control for RetroSom, default 3) \

   -g hg38 (default hg38, supporting hg38, hg19 and b37)\

   -p your_slurm_partition_name \

   -i input_type (1=sort.bam; 2=CleanBAM_and_ready_to_RetroDiscover)

   -b /CaseID.bam \

   -n 100 (maximum number of supporting reads analysis at a time, default 100) 

```

**Then run the step2 to merge the result of Case:**

```
$yourdownloadpath/Singularity_Slurm_RetroNet_step2.sh \
  
   -o /directory_path_for_output \
  
   -j CaseID \
  
   -m $yourdownloadpath \
  
   -v 3 (version control for RetroSom, default 3) \
  
   -g hg38 (default hg38, supporting hg38, hg19 and b37) \
  
   -p your_slurm_partition_name \
  
   -c ControlID_Combined \
  
   -z "N" (gpu_partition name, default "N") \
  
   -x 0.99 ( probability cutoff default 0.99) \
  
   -l 1 (number of bam you split)
```

A folder called CaseID_NoModel will be generated, there will be three kinds of output, the probability for each MEI will be in CaseID_Combined/RetroNet; the bed files for detailed MEI information will be in CaseID_Combined/retro_v3; the result svg will be in CaseID_Combined/visual/
