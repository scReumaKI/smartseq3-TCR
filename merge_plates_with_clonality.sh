#!/bin/bash
# ============================================================================ #
# merge_plates_with_clonality.sh                                               #
# Author: Juan Sebastian Diaz Boada                                            #
# Creation Date: 29/12/2022                                                    #
# ============================================================================ #
# Define the help function
function help {
  # Print the usage message
  echo "Usage: $0"
  echo "Converts all of the .bam files in INPUT_DIR into fastq files for the
        Aligned and unmapped directories. Then concatenates the 2 fastqfiles of the
        same cell."
  echo ""
  # Print a description of the script's parameters
  echo "Parameters:"
  echo "None"
  echo ""
  # Print the list of options
  echo "Options:"
  echo "  --input_dir       Relative path to the directory containing the Tracer output for each plate. Defaults to 'data/05_SS3_collected_TCRs'."
  echo "  --out_file        Relative path to the output dataset file. Defaults to 'results/TCR_clonality.tsv'."
  echo "  -h, --help        display this help and exit"
  # Exit with a success status code
  exit 0
}
# Set default values for the input folder and output file
input_dir="data/05_SS3_collected_TCRs"
out_file="results/TCR_clonality.tsv"
# Process the command-line arguments
while [ $# -gt 0 ]; do
  case "$1" in
    --input_dir)
      input_dir=$2
      shift 2
      ;;
    --out_file)
      out_file=$2
      shift 2
      ;;
    --help)
      # Print usage message and exit
      help
      ;;
    *)
      echo "Invalid argument: $1" >&2
      exit 1
      ;;
  esac
done
# ---------------------------------------------------------------------------- #
singularity exec env/01_pysam_SS3.sif ./src/06_clonality_analysis.py $input_dir $out_file
