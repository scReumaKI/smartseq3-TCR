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
│   │   │   └── Plate_1.barcodes.csv
│   │   ├── Plate_2
│   │   │   ├── Plate_2.filtered.tagged.Aligned.out.bam
│   │   │   ├── Plate_2.filtered.tagged.unmapped.bam
│   │   │   └── Plate_2.barcodes.csv
│   │   └── Plate_3
│   │   │   ├── Plate_3.filtered.tagged.Aligned.out.bam
│   │   │   ├── Plate_3.filtered.tagged.unmapped.bam
│   │   │   └── Plate_3.barcodes.csv
│   ├── 01_SS3_splitted_bams
│   ├── 02_SS3_merged_fastq
│   ├── 03_SS3_trimmed_fastq
│   ├── 04_SS3_Tracer_assembled_cells
│   └── 04_SS3_collected_TCRs
├── env
│   ├── 01_pysam_SS3.def
│   ├── 02_samtools_SS3.def
│   ├── 03_trimgalore_SS3.def
│   ├── 04_tracer_SS3.def
│   └── figlet.def
├── README.md
├── results
└── src
    ├── 01_split_bam_by_tag_and_condition_file.py
    ├── 02_bam2fastq.sh
    ├── 03_run_trim_galore.sh
    ├── 04_assemble_trimmed_cells.sh
    └── 05_collect_assemble.py
```
# TL;DR
1. Place the zUMIs output in the folder [`data/00_SS3_raw_data/`](data/00_SS3_raw_data/) named as the plate (`Plate_1` for example). Make sure that the following 3 files are in the plate subfolder:
  1.1 Bam file for aligned reads, named `Plate_1.filtered.tagged.Aligned.out.bam`.
  1.2 Bam file for unmapped reads, named `Plate_1.filtered.tagged.unmapped.bam`.
  1.3 Barcodes per cell in `.csv` format named `Plate_1.barcodes.csv`.
2. Run the pipeline as follows
```bash
./complete_pipeline.sh Plate_1
```
for all the existing plates. The TCR dataset will be saved in [`results/`](results/Plate_1.tsv).

Optional: Specify the number of nodes for parallel execution after the plate name. By default it will run on 10 nodes.

3. Clonality TO BE IMPLEMENTED




# Detailed explanation: Run one module at a time.
## 0. Build singularity containers
After cloning this repository, build the singularity images in the [`env`](env/) folder. The command `singularity build` requires administration rights or being added to a *fakeroot* list. Consult your local IT for information on how to build singularity containers. Other options include building the containers remotely or on a local device and `scp` the image to your HPC system. More information on building containers can be found [here](https://hpc.nih.gov/apps/singularity.html#create). To build the containers using the `fakeroot` option from singularity definition files:
```bash
singularity build --fakeroot env/01_pysam_SS3.sif env/01_pysam_SS3.def
singularity build --fakeroot env/02_samtools_SS3.sif env/02_samtools_SS3.def
singularity build --fakeroot env/03_trimgalore_SS3.sif env/03_trimgalore_SS3.def
singularity build --fakeroot env/04_tracer_SS3.sif env/04_tracer_SS3.def
singularity build --fakeroot env/figlet.sif env/figlet.def

```
## 1. Split the bam file
### Context
The output files from the sequence facility of interest are two big `.bam` file per plate given by zUMIs:
```bash
<plate_name>.filtered.tagged.Aligned.out.bam
<plate_name>.filtered.tagged.unmapped.bam
```
and a barcode file `<plate_name>.barcodes.csv`, a tabular file with a column for barcodes and a column with the corresponding cell name.
### How to run
Move your smartseq3 data to the folder [`data/00_SS3_raw_data/`](data/00_SS3_raw_data/) and to a sub-folder corresponding to a plate.
This section uses a [container](env/01_pysam_SS3.def) with Python 3.9 and `pysam` that calls the python script [`/src/01_split_bam_by_tag_and_condition_file.py`](/src/01_split_bam_by_tag_and_condition_file.py) internally. The script extracts one `.bam` file per cell from each big `.bam` given the barcode and name of each cell.
```bash
./env/01_pysam_SS3.sif path/to/bam path/to/barcode.csv output/dir --condition_tag_col <barcode_column> --condition_name_col <cell_name_column> --bam_tag_flag BC
```
Inputs:
| Parameter | Type | Description |
| ------ | --- | --- |
| `bam_in` | string | Relative path to the directory holding the multiplexed .bam file. |
| `condition_csv` | string | Relative path to the .csv file containing the mapping between the barcodes and the name of the cells. |
| `bam_out` | string | Relative path to the folder where the output .bam files per cell are to be saved. |
| `condition_tag_col` | string (optional) | Name of the column containing the barcodes in 'condition_csv'. Defaults to 'Barcode'. |
| `condition_name_col` | string (optional) | Name of the column containing the cell name in 'condition_csv'. Defaults to 'Name'. |
| `bam_tag_flag` | string (optional) | The tag in the bam file that contains the sample barcode. Defaults to 'BC' for zUMIs output. |
| `name_part_filer` | string (optional) | Use to limit itself to samples names that contain a particular substring. Defaults to None. |
Example:
```bash
> pwd
/srv/shared/Common_New/Transcriptomics/TCR-myositis/smartseq3
> ./env/01_pysam_SS3.sif data/00_SS3_raw_data/Plate_1/ Plate_1.filtered.tagged.Aligned.out.bam data/00_SS3_raw_data/Plate_1/Plate_1.barcodes.csv data/01_SS3_splitted_bams/Aligned/Plate_1/ --condition_tag_col Barcode --condition_name_col Name --bam_tag_flag BC
```
For detailed help, type `./env/01_pysam_SS3.sif --help` or `singularity run-help env/01_pysam_SS3.sif`
### Considerations
+ Execution time on 20 nodes: ~1:10 hour for a plate of 384 cells.
+ The previous procedure has to be done for the `.Aligned.out.bam` file and for the `.unmapped.bam` file.
+ The output directory `data/01_SS3_splitted_bams/Aligned/Plate_1/` has to be created before running the script.

## 2. Convert to fastq and concatenate:
### Context
The individual `.bam` files from the *Aligned* and *unmapped* files have to be converted to `fastq.gz` and be concatenated.
### How to run
This section uses a [container](env/02_samtools_SS3.def) with [`samtools`](https://github.com/samtools/samtools) executing the bash script [`src/02_bam2fastq.sh`](src/02_bam2fastq.sh), that translates the `.bam` files to `fastq.gz` format, and the concatenates the *Aligned* and *unmapped* per cell.
```bash
./env/02_samtools_SS3.sif path/to/splitted/bams output/dir <nodes>
```
Inputs to the bash script:
| Parameter | Type | Description |
| ------ | --- | ----- |
| `INPUT_DIR` | string | Directory with the aligned single cell bam files. |
| `OUTPUT_DIR` | string | Directory where the fastq.gz files will be written. |
| `NODES` | int | Number of nodes to use in `samtools fastq`. |
Example:
```bash
./env/02_samtools_SS3.sif data/01_SS3_splitted_bams/Plate_1/ \
data/02_SS3_merged_fastq/Plate_1/ <nodes>
```
For detailed help, type `./env/02_samtools_SS3.sif --help` or `singularity run-help env/02_samtools_SS3.sif`
### Considerations
+ Execution time on 20 nodes:
+ The script creates temporary folders that are deleted if the script terminates with exit status 0.
+ The output directory `data/02_SS3_merged_fastq/Plate_1/` has to be created before running the script.

## 3. Trimming adaptors
### Context
The untrimmed fastq files of the cells are all saved on the directory [`data/02_SS3_merged_fastq`](data/02_SS3_merged_fastq). This steps trim the adapters and saves the output of each cell in a separate folder named as the cell.
### How to run
This section uses a [container](env/03_trimgalore_SS3.def) with [TrimGalore 0.6.7](https://github.com/FelixKrueger/TrimGalore) that calls the bash script [`src/03_run_trim_galore.sh`](src/03_run_trim_galore.sh), which trims the adapters of the concatenated fastq files and creates the output in a new folder named as the cell.
```bash
./env/03_trimgalore_SS3.sif path/to/fastq/dir output/dir <nodes>
```
Inputs to the bash script:
| Parameter | Type | Description |
| ------ | --- | ----- |
| `INPUT_DIR` | string | Directory with the merged fastq files. |
| `OUTPUT_DIR` | string | Directory where the trimmed fastq.gz files will be written. |
| `NODES` | int | Number of nodes to use in `trim_galore`. |
Example:
```bash
./env/03_trimgalore_SS3.sif data/02_SS3_merged_fastq/Plate_1/ \
data/03_SS3_trimmed_fastq/Plate_1/ 8
```
For detailed help, type `./env/03_trimgalore_SS3.sif --help` or `singularity run-help env/03_trimgalore_SS3.sif`
### Considerations
+ Execution time in Reuma: < 2h for a 384 cell plate
+ Apparently trim_galore does not accept more than 8 cores. If provided more, it will truncate to 8
+ The output directory `data/03_SS3_trimmed_fastq/Plate_1/` has to be created before running the script.

## 4. TCR assembling
### Context
[TraCeR](https://github.com/Teichlab/tracer) is a package that uses [Bowtie](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml), [Trinity](https://github.com/trinityrnaseq/trinityrnaseq/wiki), [IgBlast](https://www.ncbi.nlm.nih.gov/igblast/faq.html#standalone) and [Kallisto](http://pachterlab.github.io/kallisto/), among other tools, to assemble T cell receptors (TCRs) from fastq files.
### How to run
This section uses a [container](env/04_tracer_SS3.def) with [TraCeR](https://github.com/Teichlab/tracer) that calls the bash script [`src/04_assemble_trimmed_cells.sh](`.sh`), which calls `tracer assemble` iteratively over the files of multiple cells.
```bash
./env/04_tracer_SS3.sif path/to/cells output/dir <nodes> <loci>
```
Inputs to the bash script:
| Parameter | Type | Description |
| ------ | --- | ----- |
| `INPUT_DIR` | string | Relative path to the directory containing the trimmed fastq files. |
| `OUTPUT_DIR` | string | Relative path to the directory where the TCR files will be saved. |
| `NODES` | int | Number of nodes to use in `tracer assemble`. |
| `LOCI` | str | 'AB' for assembling alpha-beta chains and 'GD' for assembling gamma-delta chains." |
Example:
```bash
./env/04_tracer_SS3.sif data/03_SS3_trimmed_fastq/Plate_1/ \
data/04_SS3_Tracer_assembled_cells/Plate_1/AB 20 'AB'
```
For detailed help, type `./env/04_tracer_SS3.sif --help` or `singularity run-help env/04_tracer_SS3.sif`
### Considerations
+ Execution time with 40 nodes: ~ 5 days for a 384 cell plate
+ The output directory `data/04_SS3_Tracer_assembled_cells/Plate_1/AB` has to be created before running the script.
+ The above commands should be run for all the loci-pairs separately, namely, for alpha-beta (`AB`) and for gamma-delta (`GD`).

## 5. Clonality Dataset
### Context

### How to run

### Considerations
+ Execution time: < 1 minute for a 384 cell plate
