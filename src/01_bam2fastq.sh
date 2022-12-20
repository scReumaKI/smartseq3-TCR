#!/bin/bash
# ============================================================================ #
# bam2fastq_concat.sh                                                          #


# Author: Juan Sebastian Diaz Boada                                            #
# Creation Date: 8/12/2022                                                     #
# ============================================================================ #
INPUT_DIR=$1
OUTPUT_DIR=$2
NODES=$3

mkdir __temporary_fastq__
mkdir __temporary_fastq__/__Aligned__
mkdir __temporary_fastq__/__unmapped__
echo "Created temporary folders"
for FILE in ${INPUT_DIR}Aligned/*.bam
do
  NAME=$(basename $FILE .bam)
  # Create fastq files for the aligned bams
  echo "Creating aligned fastq files for cell "$NAME
  samtools sort -n $FILE | samtools fastq --threads $NODES \
  -1 __temporary_fastq__/__Aligned__/${NAME}_R1.fastq.gz \
  -2 __temporary_fastq__/__Aligned__/${NAME}_R2.fastq.gz \
  -0 /dev/null -s /dev/null
  # Create fastq files for the unmapped files
  echo "Creating unmapped fastq files for cell "$NAME
  samtools sort -n ${INPUT_DIR}unmapped/${NAME}.bam | samtools fastq --threads $NODES \
  -1 __temporary_fastq__/__unmapped__/${NAME}_R1.fastq.gz \
  -2 __temporary_fastq__/__unmapped__/${NAME}_R2.fastq.gz \
  -0 /dev/null -s /dev/null

  cat __temporary_fastq__/__Aligned__/${NAME}_R1.fastq.gz > ${OUTPUT_DIR}${NAME}_R1.fastq.gz
  cat __temporary_fastq__/__Aligned__/${NAME}_R2.fastq.gz > ${OUTPUT_DIR}${NAME}_R2.fastq.gz
  cat __temporary_fastq__/__unmapped__/${NAME}_R1.fastq.gz >> ${OUTPUT_DIR}${NAME}_R1.fastq.gz
  cat __temporary_fastq__/__unmapped__/${NAME}_R2.fastq.gz >> ${OUTPUT_DIR}${NAME}_R2.fastq.gz
  echo "Concatenated Aligned and unmapped for cell "$NAME
  echo "======================================================================"
  rm __temporary_fastq__/__Aligned__/${NAME}_R1.fastq.gz
  rm __temporary_fastq__/__Aligned__/${NAME}_R2.fastq.gz
  rm __temporary_fastq__/__unmapped__/${NAME}_R1.fastq.gz
  rm __temporary_fastq__/__unmapped__/${NAME}_R2.fastq.gz
done
rm -rf __temporary_fastq__
echo "Deleted temporary folders"
