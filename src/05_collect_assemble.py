#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# =============================================================================================
# 05_collect_assemble.py
# Author: Juan Sebastian Diaz Boada
# juan.sebastian.diaz.boada@ki.se
# Creation Date: 15/03/2022
# =============================================================================================
""" Iterates over the cell folders reading TraCeR assemble files to build aa dataset.

    Parameters
    ----------
    in_path : string.
        Path to the folder holding the trimmed cells. Each cell has its own folder named under
        the new convention Plate-tissue-well. For more details on the directory structure see
        README.md in this folder.
    out_file : string.
        Path and name of the output dataframe in .csv, tsv or xlsx form.

"""
import os,sys
import argparse
import pandas as pd

module_path = os.path.abspath('bin')
if module_path not in sys.path:
    sys.path.append(module_path)

from objects import Cell, Chain
from objects import AlphaChain, BetaChain, GammaChain, DeltaChain
from objects import create_cell_from_AB, append_GD_data

parser = argparse.ArgumentParser(description='in and out paths')
parser.add_argument('in_path', type=str, help='Path of the data folder.')
parser.add_argument('out_file', type=str, help='Path of the output dataset.')
args = parser.parse_args()
in_path = args.in_path
AB_path = os.path.join(in_path,os.path.normpath('AB'))
GD_path = os.path.join(in_path,os.path.normpath('GD'))
# Loop over AB
print("#######################################################################")
print("Starting alpha-beta files reading...")
print("#######################################################################")
cells = {}
cont = 1
L = len(os.listdir(AB_path))
for folder in os.listdir(AB_path):
    file = os.path.join(AB_path,folder,os.path.normpath('filtered_TCR_seqs/filtered_TCRs.txt'))
    if os.path.exists(file):
        cell = create_cell_from_AB(file)
    else:
        cell = Cell(folder)
    print("Cell {}, {}/{}".format(cell.name,cont,L))
    cells[cell.name] = cell
    cont = cont + 1
# Loop over GD
print("#######################################################################")
print("Starting gamma-delta files reading...")
print("#######################################################################")
cont = 1
L = len(os.listdir(GD_path))
for folder in os.listdir(GD_path):
    file = os.path.join(GD_path,folder,os.path.normpath('filtered_TCR_seqs/filtered_TCRs.txt'))
    if os.path.exists(file):
        append_GD_data(file,cells)
    print("Cell {}, {}/{}".format(folder,cont,L))
    cont = cont +1
# Dataframe generation
print("#######################################################################")
print("Generating dataframe...")
print("#######################################################################")
DF = pd.DataFrame(index=cells.keys())
cont = 1
L = len(cells)
for name,cell in cells.items():
    print("Cell {}, {}/{}".format(name,cont,L))
    chains = cell.A_chains + cell.B_chains + cell.G_chains + cell.D_chains
    for chain in chains:
        meta_colnames = [chain.allele + '_' + s for s in list(chain.__dict__.keys())][2:]
        values = list(chain.__dict__.values())[2:]
        DF.loc[name,meta_colnames] = values
    cont = cont +1
# Export dataset
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
