#!/bin/bash
# ============================================================================ #
# complete_pipeline.sh                                                         #
# Author: Juan Sebastian Diaz Boada                                            #
# Creation Date: 21/12/2022                                                    #
# ============================================================================ #
# Handling incomplete parameters
if [ $# -lt 1 ]; then
  echo "Error: No plate name was given." >&2
  exit 1
fi
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
  echo "  NODES          Number of nodes to use in parallel computations. Defaults to 10."
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
shift
NODES=${1:-10}
echo "Nodes = $NODES"
# Handling missing data
if [ ! -d data/00_SS3_raw_data/${PLATE_NAME}/ ];then
  echo "There is no raw data for plate ${PLATE_NAME} in data/00_SS3_raw_data/" >&2
  exit 1
fi
# -------------------------------------------------------------------------------------
# 00. Container building
if [ ! -f env/figlet.sif ];then
  singularity build --fakeroot env/figlet.sif env/figlet.def
  echo "Built env/figlet.sif singularity container."
else
  echo "Container env/figlet.sif already exists. Using existing image..."
fi
echo "================================================================================="
./env/figlet.sif "0. Singularity"
echo "================================================================================="
if [ ! -f env/01_pysam_SS3.sif ];then
  singularity build --fakeroot env/01_pysam_SS3.sif env/01_pysam_SS3.def
  echo "Built env/01_pysam_SS3.sif singularity container."
else
  echo "Container env/01_pysam_SS3.sif already exists. Using existing image..."
fi
if [ ! -f env/02_samtools_SS3.sif ];then
  singularity build --fakeroot env/02_samtools_SS3.sif env/02_samtools_SS3.def
  echo "Built env/02_samtools_SS3.sif singularity container."
else
  echo "Container env/02_samtools_SS3.sif already exists. Using existing image..."
fi
if [ ! -f env/03_trimgalore_SS3.sif ];then
  singularity build --fakeroot env/03_trimgalore_SS3.sif env/03_trimgalore_SS3.def
  echo "Built env/03_trimgalore_SS3.sif singularity container."
else
  echo "Container env/03_trimgalore_SS3.sif already exists. Using existing image..."
fi
if [ ! -f env/04_tracer_SS3.sif ];then
  singularity build --fakeroot env/04_tracer_SS3.sif env/04_tracer_SS3.def
  echo "Built env/04_tracer_SS3.sif singularity container."
else
  echo "Container env/04_tracer_SS3.sif already exists. Using existing image..."
fi
# ------------------------------------------------------------------------------------
# 01. SPLIT BAM files
echo "================================================================================="
./env/figlet.sif "1. Pysam"
echo "================================================================================="
# Split Aligned reads file
if [ ! -d data/01_SS3_splitted_bams/${PLATE_NAME}/Aligned/ ];then
  mkdir -p data/01_SS3_splitted_bams/${PLATE_NAME}/Aligned/
  echo "Created folder for Aligned reads"
fi
./env/01_pysam_SS3.sif \
data/00_SS3_raw_data/${PLATE_NAME}/${PLATE_NAME}.filtered.tagged.Aligned.out.bam \
data/00_SS3_raw_data/${PLATE_NAME}/${PLATE_NAME}.barcodes.csv \
data/01_SS3_splitted_bams/${PLATE_NAME}/Aligned/ \
--condition_tag_col Barcode --condition_name_col Name --bam_tag_flag BC
# Split unmapped reads file
if [ ! -d data/01_SS3_splitted_bams/${PLATE_NAME}/unmapped/ ];then
  mkdir -p data/01_SS3_splitted_bams/${PLATE_NAME}/unmapped/
  echo "Created folder for unmapped reads"
fi
./env/01_pysam_SS3.sif \
data/00_SS3_raw_data/${PLATE_NAME}/${PLATE_NAME}.filtered.tagged.unmapped.bam \
data/00_SS3_raw_data/${PLATE_NAME}/${PLATE_NAME}.barcodes.csv \
data/01_SS3_splitted_bams/${PLATE_NAME}/unmapped/ \
--condition_tag_col Barcode --condition_name_col Name --bam_tag_flag BC
# ------------------------------------------------------------------------------------
# 01. Translate and merge
echo "================================================================================="
./env/figlet.sif "2. Samtools"
echo "================================================================================="
if [ ! -d data/02_SS3_merged_fastq/${PLATE_NAME}/ ];then
  mkdir -p data/02_SS3_merged_fastq/${PLATE_NAME}/
fi
./env/02_samtools_SS3.sif data/01_SS3_splitted_bams/${PLATE_NAME}/ \
data/02_SS3_merged_fastq/${PLATE_NAME}/ $NODES
# -------------------------------------------------------------------------------------
# 02. Trim adapters
echo "================================================================================="
./env/figlet.sif "3. TrimGalore!"
echo "================================================================================="
if [ ! -d data/03_SS3_trimmed_fastq/${PLATE_NAME}/ ];then
  mkdir -p data/03_SS3_trimmed_fastq/${PLATE_NAME}/
fi
./env/03_trimgalore_SS3.sif data/02_SS3_merged_fastq/${PLATE_NAME}/ \
data/03_SS3_trimmed_fastq/${PLATE_NAME}/ 8
# -------------------------------------------------------------------------------------
# 03. TCR assemble
echo "================================================================================="
./env/figlet.sif "4. TraCeR"
echo "================================================================================="
# Assemble alpha-beta
if [ ! -d data/04_SS3_Tracer_assembled_cells/${PLATE_NAME}/AB/ ];then
  mkdir -p data/04_SS3_Tracer_assembled_cells/${PLATE_NAME}/AB/
  echo "Created folder for AB TCRs"
fi
./env/04_tracer_SS3.sif data/03_SS3_trimmed_fastq/${PLATE_NAME}/ \
data/04_SS3_Tracer_assembled_cells/${PLATE_NAME}/AB $NODES 'AB'
# Assemble gamma-delta
if [ ! -d data/04_SS3_Tracer_assembled_cells/${PLATE_NAME}/GD/ ];then
  mkdir -p data/04_SS3_Tracer_assembled_cells/${PLATE_NAME}/GD/
  echo "Created folder for GD TCRs"
fi
./env/04_tracer_SS3.sif data/03_SS3_trimmed_fastq/${PLATE_NAME}/ \
data/04_SS3_Tracer_assembled_cells/${PLATE_NAME}/GD $NODES 'GD'
# -------------------------------------------------------------------------------------
# 04. TCR collection
echo "================================================================================="
./env/figlet.sif "5. TCR collection"
echo "================================================================================="
if [ ! -d data/05_SS3_collected_TCRs/${PLATE_NAME}/ ];then
  mkdir -p data/05_SS3_collected_TCRs/${PLATE_NAME}/
fi
singularity exec env/01_pysam_SS3.sif ./src/04_collect_assemble.py \
data/04_SS3_Tracer_assembled_cells/${PLATE_NAME}/ \
data/05_SS3_collected_TCRs/${PLATE_NAME}/${PLATE_NAME}.tsv
