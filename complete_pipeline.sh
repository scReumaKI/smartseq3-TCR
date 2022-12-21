#!/bin/bash
# ============================================================================ #
# complete_pipeline.sh                                                         #
# Author: Juan Sebastian Diaz Boada                                            #
# Creation Date: 21/12/2022                                                    #
# ============================================================================ #
# Define the help function
function help {
  # Print the usage message
  echo "Usage: $0 [PLATE_NAME][NODES]"
  echo "Runs the Smart-seq3 TCR extraction pipeline on a the sequencing data of a plate."
  echo "A wrapper of Pysam, samtools, TrimGalore and TraCeR"
  echo "This script assumes the existence of a directory structure as in the repo"
  echo "https://github.com/scReumaKI/smartseq3-TCR"
  echo "and the existence of the singularity container images (.sif)."
  echo ""
  # Print a description of the script's parameters
  echo "Parameters:"
  echo "  PLATE_NAME     Name of the plate as prefix of sequencing files and folder names."
  echo "  NODES          Number of nodes to use in parallel computations."
  echo ""
  # Print the list of options
  echo "Options:"
  echo "  -h, --help        display this help and exit"
  # Exit with a success status code
  exit 0
}

# Parse the options and arguments
if [ "$1" = "-h" ] || [ "$1" = "--help" ]; then
  help
fi

PLATE_NAME=$1
NODES=$2
# 00. SPLIT BAM files
echo "================================================================================="
./env/figlet.sif "00 BAM file splitting with pysam"
echo "================================================================================="
# Split Aligned reads file
if [ ! -d data/01_SS3_splitted_bams/${PLATE_NAME}/Aligned/ ];then
  mkdir data/01_SS3_splitted_bams/${PLATE_NAME}/Aligned/
  echo "Created folder for Aligned reads"
echo "./env/00_split_bam_SS3.sif data/00_SS3_raw_data/${PLATE_NAME}/${PLATE_NAME}.filtered.tagged.Aligned.out.bam \ 
data/00_SS3_raw_data/${PLATE_NAME}/${PLATE_NAME}.barcodes.csv \
data/01_SS3_splitted_bams/${PLATE_NAME}/Aligned/"
# Split unmapped reads file
if [ ! -d data/01_SS3_splitted_bams/${PLATE_NAME}/unmapped/ ];then
  mkdir data/01_SS3_splitted_bams/${PLATE_NAME}/unmapped/
  echo "Created folder for unmapped reads"
echo "./env/00_split_bam_SS3.sif data/00_SS3_raw_data/${PLATE_NAME}/${PLATE_NAME}.filtered.tagged.unmapped.bam \
data/00_SS3_raw_data/${PLATE_NAME}/${PLATE_NAME}.barcodes.csv \
data/01_SS3_splitted_bams/${PLATE_NAME}/unmapped/"
# 01. Translate and merge
echo "================================================================================="
./env/figlet.sif "01 Translate to fastq and merge with samtools"
echo "================================================================================="
echo "./env/01_merge_fastq.sif data/01_SS3_splitted_bams/${PLATE_NAME} \
data/02_SS3_merged_fastq/${PLATE_NAME}/ $NODES"
# 02. Trim adapters
echo "================================================================================="
./env/figlet.sif "02 Trim adapters with TrimGalore!"
echo "================================================================================="
echo "./env/02_trim_adapters.sif data/02_SS3_merged_fastq/${PLATE_NAME}/ \
data/03_SS3_trimmed_fastq/${PLATE_NAME}/ $NODES"
# 03. TCR assemble
echo "================================================================================="
./env/figlet.sif "03 Assemble TCR with TraCeR"
echo "================================================================================="
echo "./env/03_assemble_TCR.sif data/02_SS3_merged_fastq/${PLATE_NAME}/ \
data/03_SS3_trimmed_fastq/${PLATE_NAME}/ $NODES 'AB'"
echo "./env/03_assemble_TCR.sif data/02_SS3_merged_fastq/${PLATE_NAME}/ \
data/03_SS3_trimmed_fastq/${PLATE_NAME}/ $NODES 'AB'"
