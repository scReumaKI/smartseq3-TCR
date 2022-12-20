# TCR extraction from Smart-seq3 sequencing data
Instructions fro running smart-seq3 TCR pipeline.

### Expected directory structure
```bash
├── data
│   ├── 00_SS3_output
│   │   ├── SS3_21_231
│   │   │   ├── SS3_21_231.filtered.tagged.Aligned.out.bam
│   │   │   ├── SS3_21_231.filtered.tagged.unmapped.bam
│   │   ├── SS3_21_233
│   │   │   ├── SS3_21_233.filtered.tagged.Aligned.out.bam
│   │   │   ├── SS3_21_233.filtered.tagged.unmapped.bam
│   │   └── SS3_21_235
│   │   │   ├── SS3_21_235.filtered.tagged.Aligned.out.bam
│   │   │   ├── SS3_21_235.filtered.tagged.unmapped.bam
│   ├── 01_SS3_splitted_bams
│   │   ├── SS3_21_231
│   │   │   ├── Aligned
│   │   │   └── unmapped
│   │   └── SS3_21_233
│   │       ├── Aligned
│   │       └── unmapped
│   ├── 02_SS3_merged_fastq
│   │   ├── SS3_21_231
│   │   └── SS3_21_233
│   └── 03_SS3_trimmed_fastq
│       ├── SS3_21_231
│       └── SS3_21_233
├── env
│   ├── 00_split_bam_SS3.def
│   ├── 00_split_bam_SS3.sif
│   ├── samtools_1.16.sif
│   └── trimgalore6.7.sif
└── src
    ├── bam2fastq.sh
    ├── run_trim_galore.sh
    └── split_bam_by_tag_and_condition_file.py
```

## Smartseq3
### 1. Split the bam file
#### Context
The output files from the sequence facility are two big `.bam` file per plate:
```bash
<plate>.filtered.tagged.Aligned.out.bam
<plate>.filtered.tagged.unmapped.bam
```
which are located in `data/00_SS3_output/<plate>/`
This step extracts one `.bam` file per cell from each big `.bam` given the barcode and name of each cell.
#### How to run
```
./container.sif path/to/bam/ BC path/to/barcode/dataframe.csv barcode_column name_column path/to/output
```
The container is `00_split_bam_SS3.sif` which calls the python script `split_bam_by_tag_and_condition_file.py` internally.

Inputs:
+ Path to the big bam file
+ String literal 'BC' (for barcode. For further information see [here](https://samtools.github.io/hts-specs/SAMv1.pdf))
+ Path to the `.csv` with the mapping between barcodes (in one columns) and cell names (in another column).
+ Name of the barcode column
+ Name of the cell name column
+ Output directory

Example:
```bash
> pwd
/srv/shared/Common_New/Transcriptomics/TCR-myositis
> ./env/00_split_bam_SS3.sif data/00_SS3_output/SS3_21_231 SS3_21_231.filtered.tagged.Aligned.out.bam BC data/00_SS3_output/SS3_21_231/P231_barcodes.csv Barcode Name data/01_SS3_splitted_bams/Aligned/SS3_21_231
```
#### Considerations
+ Execution time on Reuma: ~1 hour for a plate of 384 cells.
+ The previous procedure has to be done for the `.Aligned.out.bam` file and for the `.unmapped.bam` file.
+ Multiple instances of this code can be run in parallel on the same machine without interfering with each other

### 2. Convert to fastq and concatenate:
#### Context
The individual `.bam` files from the *Alignes* and *unmapped* files have to be converted to `fastq.gz` and concatenated.
#### How to run
This section uses a `samtools` singularity container (`samtools_1.16.sif`) executing the bash script `bam2fastq.sh`
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
