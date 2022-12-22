#!/bin/bash
# ============================================================================ #
# bam2fastq.sh                                                                 #
# Author: Juan Sebastian Diaz Boada                                            #
# Creation Date: 8/12/2022                                                     #
# ============================================================================ #

# Define the help function
function help {
  # Print the usage message
  echo "Usage: $0 [INPUT_DIR][OUTPUT_DIR][NODES]"
  echo "Converts all of the .bam files in INPUT_DIR into fastq files for the
        Aligned and unmapped directories. Then concatenates the 2 fastqfiles of the
        same cell."
  echo ""
  # Print a description of the script's parameters
  echo "Parameters:"
  echo "  INPUT_DIR     Relative path to the directory containing the single-cell .bam files."
  echo "  OUTPUT_DIR    Relative path to the directory where the fastq files will be saved."
  echo "  NODES         Number of nodes for samtools."
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

INPUT_DIR=$1
OUTPUT_DIR=$2
NODES=$3

TEMP_DIR=${OUTPUT_DIR}__temporary_fastq__
mkdir $TEMP_DIR
mkdir $TEMP_DIR/__Aligned__
mkdir $TEMP_DIR/__unmapped__
echo "Created temporary folders"
for FILE in ${INPUT_DIR}Aligned/*.bam;do
  NAME=$(basename $FILE .bam)
  # Create fastq files for the aligned bams
  echo "Creating aligned fastq files for cell "$NAME
  samtools sort -n $FILE | samtools fastq --threads $NODES \
  -1 $TEMP_DIR/__Aligned__/${NAME}_R1.fastq.gz \
  -2 $TEMP_DIR/__Aligned__/${NAME}_R2.fastq.gz \
  -0 /dev/null -s /dev/null
  # Create fastq files for the unmapped files
  echo "Creating unmapped fastq files for cell "$NAME
  samtools sort -n ${INPUT_DIR}unmapped/${NAME}.bam | samtools fastq --threads $NODES \
  -1 $TEMP_DIR/__unmapped__/${NAME}_R1.fastq.gz \
  -2 $TEMP_DIR/__unmapped__/${NAME}_R2.fastq.gz \
  -0 /dev/null -s /dev/null
  # Concatenate Aligned and unmapped fastq files
  mkdir ${OUTPUT_DIR}${NAME}
  cat $TEMP_DIR/__Aligned__/${NAME}_R1.fastq.gz > ${OUTPUT_DIR}${NAME}/${NAME}_R1.fastq.gz
  cat $TEMP_DIR/__Aligned__/${NAME}_R2.fastq.gz > ${OUTPUT_DIR}${NAME}/${NAME}_R2.fastq.gz
  cat $TEMP_DIR/__unmapped__/${NAME}_R1.fastq.gz >> ${OUTPUT_DIR}${NAME}/${NAME}_R1.fastq.gz
  cat $TEMP_DIR/__unmapped__/${NAME}_R2.fastq.gz >> ${OUTPUT_DIR}${NAME}/${NAME}_R2.fastq.gz
  echo "Concatenated Aligned and unmapped for cell "$NAME
  echo "======================================================================"
  rm $TEMP_DIR/__Aligned__/${NAME}_R1.fastq.gz
  rm $TEMP_DIR/__Aligned__/${NAME}_R2.fastq.gz
  rm $TEMP_DIR/__unmapped__/${NAME}_R1.fastq.gz
  rm $TEMP_DIR/__unmapped__/${NAME}_R2.fastq.gz
done
rm -rf $TEMP_DIR
echo "Deleted temporary folders"
