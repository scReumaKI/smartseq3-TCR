#!/bin/bash
# ============================================================================ #
# assemble_trimmed_cells.sh                                                    #
# Author: Juan Sebastian Diaz Boada                                            #
# Creation Date: 13/12/2022                                                    #
# ============================================================================ #
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
  R1=${DIR}/${CELL}_R1_trimmed.fq.gz
  R2=${DIR}/${CELL}_R2_trimmed.fq.gz

  # if [ ! -d $OUTPUT_DIR$CELL ];then
  #   mkdir $OUTPUT_DIR$CELL
  #   echo "Created folder for cell $CELL"
  # else
  #   echo "Writting on existing directory $OUTPUT_DIR$CELL"
  # fi
  tracer assemble --loci $LOCI -p $NODES -s Hsap $R1 $R2 $CELL $OUTPUT_DIR
done
