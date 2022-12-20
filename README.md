# TCR extraction from Smart-seq3 sequencing data
Pipeline for extracting TCR using [TraCeR](https://github.com/Teichlab/tracer) from Smart-seq3 sequencing data processed by [zUMIs](https://github.com/sdparekh/zUMIs).
### Requisites
A HPC cluster running a version of Linux with [Singularity](https://sylabs.io/singularity/) installed.

### Expected directory structure (for 3 plates)
```bash
├── data
│   ├── 00_SS3_raw_data
│   │   ├── Plate_1
│   │   │   ├── Plate_1.filtered.tagged.Aligned.out.bam
│   │   │   ├── Plate_1.filtered.tagged.unmapped.bam
│   │   │   └── Plate_1_barcodes.csv
│   │   ├── Plate_2
│   │   │   ├── Plate_2.filtered.tagged.Aligned.out.bam
│   │   │   ├── Plate_2.filtered.tagged.unmapped.bam
│   │   │   └── Plate_2_barcodes.csv
│   │   └── Plate_3
│   │   │   ├── Plate_3.filtered.tagged.Aligned.out.bam
│   │   │   ├── Plate_3.filtered.tagged.unmapped.bam
│   │   │   └── Plate_3_barcodes.csv
│   ├── 01_SS3_splitted_bams
│   │   ├── Plate_1
│   │   │   ├── Aligned
│   │   │   └── unmapped
│   │   ├── Plate_2
│   │   │   ├── Aligned
│   │   │   └── unmapped
│   │   └── Plate_3
│   │       ├── Aligned
│   │       └── unmapped
│   ├── 02_SS3_merged_fastq
│   │   ├── Plate_1
│   │   ├── Plate_2
│   │   └── Plate_3
│   └── 03_SS3_trimmed_fastq
│   │   ├── Plate_1
│   │   ├── Plate_2
│   │   └── Plate_3
│   └── 04_SS3_Tracer_assembled_cells
│       ├── Plate_1
│       ├── Plate_2
│       └── Plate_3
├── env
│   ├── 00_split_bam_SS3.def
│   ├── 01_merge_fastq.def
│   ├── 02_trim_adapters.def
│   └── 03_assemble_TCR.def
├── README.md
├── results
└── src
    ├── 00_split_bam_by_tag_and_condition_file.py
    ├── 01_bam2fastq.sh
    ├── 02_run_trim_galore.sh
    └── test_individual_bam.ipynb
```

# Smartseq3
## 0. Build singularity containers
After cloning this repository, build the singularity images in the [`env`](env/) folder using the singularity definition files:
```bash
singularity build --fakeroot env/00_split_bam_SS3.sif env/00_split_bam_SS3.def
singularity build --fakeroot env/01_merge_fastq.sif env/01_merge_fastq.def
singularity build --fakeroot env/02_trim_adapters.sif env/02_trim_adapters.def
singularity build --fakeroot env/03_assemble_TCR.sif env/03_assemble_TCR.def
```
Then move your smartseq3 data to the folder [`data/00_SS3_raw_data/`](data/00_SS3_raw_data/) and to a sub-folder corresponding to a plate.
## 1. Split the bam file
### Context
The output files from the sequence facility are two big `.bam` file per plate:
```bash
<plate>.filtered.tagged.Aligned.out.bam
<plate>.filtered.tagged.unmapped.bam
```
which are located in `data/00_SS3_raw_data/<plate>/`
This step extracts one `.bam` file per cell from each big `.bam` given the barcode and name of each cell.
#### How to run
```
./env/00_split_bam_SS3.sif path/to/bam/ path/to/barcode/dataframe.csv path/to/output --condition_tag_col barcode_column --condition_name_col name_column --bam_tag_flag BC
```
The container calls the python script `/src/00_split_bam_by_tag_and_condition_file.py` internally.

For detailed help, type `./env/00_split_bam_SS3.sif --help` or `singularity run-help env/00_split_bam_SS3.sif`

Example:
```bash
> pwd
/srv/shared/Common_New/Transcriptomics/TCR-myositis/smartseq3
> ./env/00_split_bam_SS3.sif data/00_SS3_raw_data/SS3_21_231 SS3_21_231.filtered.tagged.Aligned.out.bam data/00_SS3_raw_data/SS3_21_231/P231_barcodes.csv data/01_SS3_splitted_bams/Aligned/SS3_21_231/ Barcode Name BC
```
#### Considerations
+ Execution time on Reuma: ~1 hour for a plate of 384 cells.
+ The previous procedure has to be done for the `.Aligned.out.bam` file and for the `.unmapped.bam` file.
+ Multiple instances of this code can be run in parallel on the same machine without interfering with each other

### 2. Convert to fastq and concatenate:
#### Context
The individual `.bam` files from the *Aligned* and *unmapped* files have to be converted to `fastq.gz` and concatenated.
#### How to run
This section uses a [`samtools`](https://github.com/samtools/samtools) singularity container (`samtools_1.16.sif`) executing the bash script `bam2fastq.sh`
```bash
singularity exec container.sif bash_script.sh <parameters_of_bas_script>
```
Inputs to the bash script
+ INPUT_DIR: Directory with the aligned single cell bam files
+ OUTPUT_DIR: Directory where the fastq.gz files will be written
+ NODES: Number of nodes to use in `samtools fastq`

```bash
singularity exec env/samtools_1.16.sif src/bam2fastq.sh data/01_SS3_splitted_bams/SS3_21_231/ data/02_SS3_merged_fastq/SS3_21_231/ 40
```
#### Considerations
+ Execution time on Reuma:
+ The script `bam2fastq.sh` creates and deletes temporary folders. Avoid running multiple instances (for different batches, for example) of this script simultaneously on the same machine as one job may delete files of another.

### 3. Trimming adaptors
#### Context
The untrimmed fastq files of the cells are all saved on the directory `02_SS3_merged_fastq`. This steps trim the adaptors and saves the output of each cell in a separate folder named as the cell.
#### How to run
This section runs a bash script (`.sh`) that calls a container with TrimGalore 0.6.7 and creates the output in a new folder named as the cell.
```bash
run_trim_galore.sh path/to/merged output/path n_nodes
```
Inputs to the bash script
+ INPUT_DIR: Directory with the merged fastq files
+ OUTPUT_DIR: Directory where the trimmed fastq.gz files will be written
+ NODES: Number of nodes to use in `trim_galore`
```bash
./src/run_trim_galore.sh data/02_SS3_merged_fastq/SS3_21_233/ data/03_SS3_trimmed_fastq/SS3_21_233 8
```
#### Considerations
+ Execution time in Reuma: < 2h per plate
+ PLates can be executed in parallel
+ This script will **overwrite** any previous trimming of the input fastq files cell-wise. This means that if there is an existing trimming of different cells, these files would remain untouched, it would only overwrite the files of the cells inputed to the script.
+ Apparently trim_galore does not accept more than 8 cores. If provided more, it will truncate to 8
### 4. TCR assembling
