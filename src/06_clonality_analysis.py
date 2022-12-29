#!/usr/bin/env python3
# coding: utf-8
# =============================================================================================
# 06_clonality_analysis.py
# Author: Juan Sebastian Diaz Boada
# juan.sebastian.diaz.boada@ki.se
# Creation Date: 28/12/2022
# =============================================================================================
""" Calculates the TCR clonality of one or more TCR datasets.

    Merges all TCR datasets found in 'in_path', finds out the different AB and GD TCR
    instances, calculates their frequency and groups them in clones, before exporting
    the data in tabular form.

    Parameters
    ----------
    in_path : string.
        Path to the nested folder with the TCR datasets. A folder per plate and
        one dataset per folder is expected.
    out_file : string.
        Path and name of the output dataframe in .csv, tsv or xlsx form.
"""
import os
import sys
import argparse
import numpy as np
import pandas as pd

module_path = os.path.abspath('bin')
if module_path not in sys.path:
    sys.path.append(module_path)

from data_functions import read_dataframe, group_with_freq
from data_functions import generate_clone_sets, concat_seqs_in_set

parser = argparse.ArgumentParser(description='in and out paths')
parser.add_argument('in_path', type=str, help='Path of the data folder.')
parser.add_argument('out_file', type=str, help='Path of the output dataset.')
args = parser.parse_args()
in_path = args.in_path
# ---------------------------------------------------------------------------- #
# 1. MERGE DATASETS
DF = pd.DataFrame()
for plate in os.listdir(in_path):
    if plate != '.gitkeep':
        for file in os.listdir(os.path.join(in_path,plate)):
            df = read_dataframe(os.path.join(in_path,plate,file))
            df.insert(0,'Plate',plate)
            DF = pd.concat([DF,df])
del(df)
# ---------------------------------------------------------------------------- #
# 2. DATA TREATMENT
# Replace Nans for zero in productive columns
cols = DF.columns[DF.columns.str.endswith('productive')|DF.columns.str.endswith('stop_codon')|DF.columns.str.endswith('in_frame')]
for i in cols:
    DF.loc[:,i] = DF.loc[:,i].fillna(0).astype(int)
# Fill missing data
loci = ['A_1','A_2','B_1','B_2','G_1','G_2','D_1','D_2']
for l in loci:
    if not np.any(DF.columns.str.contains(l)):
        DF.insert(len(DF.columns),l+'_productive',0)
        DF.insert(len(DF.columns),l+'_TPM',np.nan)
        DF.insert(len(DF.columns),l+'_stop_codon',0)
        DF.insert(len(DF.columns),l+'_in_frame',0)
        DF.insert(len(DF.columns),l+'_ID',np.nan)
        DF.insert(len(DF.columns),l+'_CDR3nt',np.nan)
        DF.insert(len(DF.columns),l+'_CDR3aa',np.nan)
        DF.insert(len(DF.columns),l+'_V',np.nan)
        DF.insert(len(DF.columns),l+'_J',np.nan)
        if l in ['B_1','B_2','D_1','D_2']:
            DF.insert(len(DF.columns),l+'_D',np.nan)
# 1.3 Reorder columns
loci = ['A_1','A_2','B_1','B_2','G_1','G_2','D_1','D_2']
new_cols = list([DF.columns[0]])
for l in loci:
    new_cols = new_cols + list(DF.columns[DF.columns.str.startswith(l)])
DF = DF[new_cols]
# ---------------------------------------------------------------------------- #
# 3. PRODUCTIVE COLUMNS
# Productive dataframe
P = DF.loc[:,DF.columns.str.endswith('productive')]
P
# Insert productive columns per loci
loci = ['A','B','G','D']
for l in loci:
    if l in ['A','G']:
        suffix = '_2_J'
    elif l in ['B','D']:
        suffix = '_2_D'
    idx = int(np.where(DF.columns==l+suffix)[0][0])
    DF.insert(idx+1,l+'_productive',P.loc[:,P.columns.str.startswith(l)].sum(axis=1))
# Insert productive columns per loci pair
idx = int(np.where(DF.columns=='B_productive')[0][0])
DF.insert(idx+1,'AB_productive',P.iloc[:,:4].sum(axis=1))
DF.insert(len(DF.columns),'GD_productive',P.iloc[:,4:8].sum(axis=1))
# ---------------------------------------------------------------------------- #
# 4. MASKING SEQUENCES BY PRODUCTIVITY
# CDR3 dataframes
CDR3nt = DF.loc[:,DF.columns.str.endswith('CDR3nt')]
CDR3aa = DF.loc[:,DF.columns.str.endswith('CDR3aa')]
# AB
AB_CDR3nt = CDR3nt.iloc[:,:4]
AB_CDR3aa = CDR3aa.iloc[:,:4]
AB_mask = P.iloc[:,:4].astype(bool)
# Nucleotide masking
AB_mask.columns = AB_CDR3nt.columns
masked_ABnt = AB_CDR3nt.mask(~AB_mask)
#  Amino acid masking
AB_mask.columns = AB_CDR3aa.columns
masked_ABaa = AB_CDR3aa.mask(~AB_mask)
# GD
GD_CDR3nt = CDR3nt.iloc[:,4:8]
GD_CDR3aa = CDR3aa.iloc[:,4:8]
GD_mask = P.iloc[:,4:8].astype(bool)
# Nucleotide masking
GD_mask.columns = GD_CDR3nt.columns
masked_GDnt = GD_CDR3nt.mask(~GD_mask)
# Amino acid masking
GD_mask.columns = GD_CDR3aa.columns
masked_GDaa = GD_CDR3aa.mask(~GD_mask)
# ---------------------------------------------------------------------------- #
# 5. CLONE DEFINITION
# ABnt
cols = ['A_1_CDR3nt','A_2_CDR3nt','B_1_CDR3nt','B_2_CDR3nt']
seq_set_ABnt = generate_clone_sets(masked_ABnt,cols)
tcr_ABnt = concat_seqs_in_set(seq_set_ABnt)
# ABaa
cols = ['A_1_CDR3aa','A_2_CDR3aa','B_1_CDR3aa','B_2_CDR3aa']
seq_set_ABaa = generate_clone_sets(masked_ABaa,cols)
tcr_ABaa = concat_seqs_in_set(seq_set_ABaa)
# GDnt
cols = ['G_1_CDR3nt','G_2_CDR3nt','D_1_CDR3nt','D_2_CDR3nt']
seq_set_GDnt = generate_clone_sets(masked_GDnt,cols)
tcr_GDnt = concat_seqs_in_set(seq_set_GDnt)
# GDaa
cols = ['G_1_CDR3aa','G_2_CDR3aa','D_1_CDR3aa','D_2_CDR3aa']
seq_set_GDaa = generate_clone_sets(masked_GDaa,cols)
tcr_GDaa = concat_seqs_in_set(seq_set_GDaa)
# ---------------------------------------------------------------------------- #
# 6. CLONE GROUPING AND FREQUENCY CALCULATION
# ABnt
DF.insert(len(DF.columns),'TCR_AB_nt',tcr_ABnt)
DF = group_with_freq(DF,'TCR_AB_nt',group_unique=False,new_name='clone_ABnt')
# ABaa
DF.insert(len(DF.columns),'TCR_AB_aa',tcr_ABaa)
DF = group_with_freq(DF,'TCR_AB_aa',group_unique=False,new_name='clone_ABaa')
# GDnt
DF.insert(len(DF.columns),'TCR_GD_nt',tcr_GDnt)
DF = group_with_freq(DF,'TCR_GD_nt',group_unique=False,new_name='clone_GDnt')
# GDaa
DF.insert(len(DF.columns),'TCR_GD_aa',tcr_GDaa)
DF = group_with_freq(DF,'TCR_GD_aa',group_unique=False,new_name='clone_GDaa')
# ---------------------------------------------------------------------------- #
# DATA EXPORTING
file_type = args.out_file.split('.')[-1]
if file_type == 'tsv':
    DF.to_csv(args.out_file,sep='\t')
elif file_type == 'xlsx':
    DF.to_excel(args.out_file)
elif file_type == 'csv':
    DF.to_csv(args.out_file,sep=',')
else:
    raise NameError("Invalid output format. Has to be either .tsv, .csv or .xlsx.")
print("Dataset exported in {}".format(args.out_file))
