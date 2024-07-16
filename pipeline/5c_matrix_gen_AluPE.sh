#!/bin/bash
usage="$(basename "$0") [-h] [-o j m v g s f] -- Create feature matrix for the Alu PE supporting reads. Requires one core, and 8gb memory.

where:
    -h  show this help text
    -o  output folder path
    -j  subject ID
    -m  masterpath
    -v  version control for RetroSom (default 3)
    -g  reference genome (default hg38, supporting hg38, hg19 and b37)
    -s  strand (0, -strand; 1, +strand)"

while getopts ":ho:j:m:v:g:s:" opt; do
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
    g) hg="$OPTARG"
       ;;
    s) strand="$OPTARG"
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
#

echo "$sub ran on $(hostname -s)"

$masterpath/ALU/01_PE_matrix/01_AluPE_matrix.pl $outpath $sub $ver $strand $hg $masterpath
