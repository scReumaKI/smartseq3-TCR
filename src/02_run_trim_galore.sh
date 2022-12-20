#!/bin/bash
# ============================================================================ #
# run_trim_galore.sh                                                           #


# Author: Juan Sebastian Diaz Boada                                            #
# Creation Date: 13/12/2022                                                    #
# ============================================================================ #
INPUT_DIR=$1
OUTPUT_DIR=$2
NODES=$3

for FILE in ${INPUT_DIR}*.fastq.gz
do
  TEMP=$(basename $FILE)
  CELL="${TEMP%_*}"
  echo "======================================================================="
  # Create cell folder if not existing
  if [ ! -d $OUTPUT_DIR$CELL ];then
    mkdir $OUTPUT_DIR$CELL
    echo "Created folder for cell $CELL"
  else
    echo "Writting on existing directory $OUTPUT_DIR$CELL"
  fi
  ./env/trimgalore6.7.sif trim_galore $FILE -o $OUTPUT_DIR$CELL --cores $NODES
  echo "Trimmed cell $CELL successfully"
done
#trim_galore --cores $NODES -o /ou6put/folder/
