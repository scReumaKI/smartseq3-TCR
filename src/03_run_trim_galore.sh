#!/bin/bash
# ============================================================================ #
# 03_run_trim_galore.sh                                                        #
# Author: Juan Sebastian Diaz Boada                                            #
# Creation Date: 13/12/2022                                                    #
# ============================================================================ #
# Define the help function
function help {
  # Print the usage message
  echo "Usage: $0 [INPUT_DIR][OUTPUT_DIR][NODES]"
  echo "Runs TrimGalore over all cell fastq files."
  echo ""
  # Print a description of the script's parameters
  echo "Parameters:"
  echo "  INPUT_DIR     Relative path to the directory containing the single-cell folders with fastq files."
  echo "  OUTPUT_DIR    Relative path to the directory where the trimmed fastq files will be saved."
  echo "  NODES         Number of nodes for TrimGalore."
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
# ---------------------------------------------------------------------------- #
INPUT_DIR=$1
OUTPUT_DIR=$2
NODES=$3

for DIR in ${INPUT_DIR}/*;do
  CELL=$(basename $DIR)
  echo "======================================================================="
  # Create cell folder if not existing
  if [ ! -d $OUTPUT_DIR$CELL ];then
    mkdir $OUTPUT_DIR$CELL
    echo "Created folder for cell $CELL"
  else
    echo "Writting on existing directory $OUTPUT_DIR$CELL"
  fi
  trim_galore --paired ${INPUT_DIR}/${CELL}/${CELL}_R1.fastq.gz ${INPUT_DIR}/${CELL}/${CELL}_R2.fastq.gz -o $OUTPUT_DIR$CELL --cores $NODES
  echo "Trimmed cell $CELL successfully"
done
