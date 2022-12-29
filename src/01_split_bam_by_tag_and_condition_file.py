#!/usr/bin/env python3
# coding: utf-8
# =========================================================================== #
# 01_split_bam_by_tag_and_condition_file.py                                   #
# Authors: Daniel Ramsköld                                                    #
#          Juan Sebastian Diaz Boada                                          #
# Creation Date: 15/11/22                                                     #
# =========================================================================== #
""" Splits a multiplexed .bam file into cell-wise .bam files.

    Loads a multiplexed .bam file 'bam_in' from zUMIs smart-seq3 sequencing protocol and
    separates it into .bam files per cell using the barcodes in 'condition_tag_col' and the
    cell names 'condition_name_col' from file 'condition_csv'.


    Parameters
    ----------
    bam_in : string.
        Relative path to the directory holding the multiplexed .bam file.
    condition_csv : string
        Relative path to the .csv file containing the mapping between the barcodes and the name of the cells.
    bam_out : string
        Relative path to the folder where the output .bam files per cell are to be saved.
    condition_tag_col : string (optional)
        Name of the column containing the barcodes in 'condition_csv'. Defaults to 'Barcode'.
    condition_name_col : string (optional)
        Name of the column containing the cell name in 'condition_csv'. Defaults to 'Name'.
    bam_tag_flag : string (optional)
        The tag in the bam file that contains the sample barcode. Defaults to 'BC' for zUMIs output.
    name_part_filer : string (optional)
        Use to limit itself to samples names that contain a particular substring. Defaults to None.

"""
import pandas, argparse, pysam, tqdm

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('bam_in',type=str,help="Relative path to the directory holding the multiplexed .bam file.")
    parser.add_argument('condition_csv',type=str,help="Relative path to the .csv file containing the mapping between the barcodes and the name of the cells.")
    parser.add_argument('bam_out',type=str,help="Relative path to the folder where the output .bam files per cell are to be saved.")
    parser.add_argument('--condition_tag_col',type=str,default='Barcode',help="Name of the column containing the barcodes in 'condition_csv'. Defaults to 'Barcode'.")
    parser.add_argument('--condition_name_col',type=str,default='Name',help="Name of the column containing the cell name in 'condition_csv'. Defaults to 'Name'.")
    parser.add_argument('--bam_tag_flag',type=str,default='BC',help="The tag in the bam file that contains the sample barcode. Defaults to 'BC' for zUMIs output.")
    parser.add_argument('--name_part_filter',type=str,default=None,help="Use to limit itself to samples names that contain a particular substring. Defaults to None.")
    o = parser.parse_args()

    CT = pandas.read_csv(o.condition_csv)
    if o.name_part_filter is not None:
        CT = CT.loc[CT[o.condition_name_col].str.contains(o.name_part_filter),:]
    tag_to_name = dict(zip(CT.loc[:,o.condition_tag_col], CT.loc[:,o.condition_name_col].astype(str)))

    bam_in = pysam.AlignmentFile(o.bam_in, "rb",check_sq=False)
    name_to_bam = {name:pysam.AlignmentFile(o.bam_out + name + '.bam', "wb", template=bam_in) for name in set(tag_to_name.values())}

    for aln in tqdm.tqdm(bam_in.fetch(until_eof=True)):
        barcode = aln.get_tag(o.bam_tag_flag)
        if barcode in tag_to_name.keys():
            name_to_bam[tag_to_name[barcode]].write(aln)
