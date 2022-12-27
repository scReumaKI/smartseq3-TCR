#!/bin/bash
# ============================================================================ #
# 04_assemble_trimmed_cells.sh                                                    #
# Author: Juan Sebastian Diaz Boada                                            #
# Creation Date: 13/12/2022                                                    #
# ============================================================================ #
# Define the help function
function help {
  # Print the usage message
  echo "Usage: $0 [INPUT_DIR][OUTPUT_DIR][NODES][LOCI]"
  echo "Runs TraCeR over all trimmed fastq files."
  echo ""
  # Print a description of the script's parameters
  echo "Parameters:"
  echo "  INPUT_DIR     Relative path to the directory containing the trimmed fastq files."
  echo "  OUTPUT_DIR    Relative path to the directory where the TCR files will be saved."
  echo "  NODES         Number of nodes for TraCeR."
  echo "  LOCI          'AB' for assembling alpha-beta chains and 'GD' for assembling gamma-delta chains."
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

if [ $4 == "AB" ]; then
  LOCI="A B"
elif [ $4 == "GD" ]; then
  LOCI="G D"
else
  echo "Invalid loci. It has to be either AB or GD"
  exit 1
fi

for DIR in ${INPUT_DIR}*;do
  CELL=$(basename $DIR)
  R1=${DIR}/${CELL}_R1_val_1.fq.gz
  R2=${DIR}/${CELL}_R2_val_2.fq.gz
  echo "================================================================="
  echo "Running TraCeR for cell $CELL"
  echo "================================================================="
  tracer assemble --loci $LOCI -p $NODES -s Hsap $R1 $R2 $CELL $OUTPUT_DIR
done
