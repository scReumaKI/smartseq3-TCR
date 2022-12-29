""" Specific functions for TCR dataset treatment in pandas.

        * read_dataframe
        * group_with_freq
        * generate_clone_sets
        * group_sets
        * concat_seqs_in_set

    Authors: Juan Sebastian Diaz Boada
             juan.sebastian.diaz.boada@ki.se

    25/03/22
"""
import numpy as np
import pandas as pd
from collections import defaultdict
from itertools import product
#---------------------------------------------------------------------------------------------------#
def read_dataframe(in_file):
    """ Generalizes imports in pandas, independent of the input's type or extension.

        A wrapper function to import a generalized tabular dataset to pandas without
        the need to specify its type. It accepts `.tsv`, `.csv` and `.xlsx` formats.

        Parameters
        ----------
        in_file : string
            Path to the dataset in `.tsv`, `.csv` or `.xlsx` formats.

        Returns

        -------
        pd.DataFrame
            Imported dataframe.
    """
    file_type = in_file.split('.')[-1]
    if file_type == 'tsv':
        return pd.read_csv(in_file,sep='\t',index_col=0)
    elif file_type == 'xlsx':
        return pd.read_excel(in_file,index_col=0)
    elif file_type == 'csv':
        return pd.read_csv(in_file,sep=',',index_col=0)
    else:
        raise NameError("Invalid input format. Has to be either .tsv, .csv or .xlsx.")
#---------------------------------------------------------------------------------------------------#
def group_with_freq(df,col,group_unique=False,new_name=None):
    """ Groups identical values and calculates their frequency, returning an updated dataframe.

        Calculates the frequency of each value in the column named 'col' from
        dataframe 'df', adding it into a column named 'freq_'+col. It also assigns a group number
        to each unique value in the column 'group_'+col. If the parameter 'group_unique'
        is True, it groups sequences that appear only once in one group labelled as -1. The suffix
        of the new column is by default the name of 'col', but can be changed adding a string as
        'new_name' parameter.

        Parameters
        ----------
        df : pd.DataFrame
            Dataframe with the column of values to group
        col : string
            Name of the column holding the values to analyze.
        group_unique : bool, optional
            Wether to group samples that only appear once in a group labelled '-1' or mantain each
            unique element as a separate group. Default is False.
        new_name : string, optional
            Suffix of the neame of the new column, instead of using the name of the old column.
            Default is None.

        Returns
        -------
        pd.DataFrame
            Dataframe with additional columns for group number and frequency.
    """
    DF = df.copy()
    freq_col_name = 'freq_' + col if new_name is None else 'freq_' + new_name
    group_col_name = 'group_'+col if new_name is None else 'group_' + new_name

    DF[freq_col_name] = DF[col].map(DF[col].value_counts()).astype(pd.Int64Dtype())
    DF.sort_values(by=[freq_col_name,col],ascending=False,inplace=True)
    if group_unique:
        # Default dict for cluster numbers. Return -1 if unseen instance
        seq2idx = defaultdict(lambda : -1)
        n_rep = len(DF.loc[DF[freq_col_name]!=1])
        seqs = DF.loc[:,col].iloc[:n_rep].unique()
    else:
        seq2idx = {}
        seqs = DF[col].unique()

    for n,g in enumerate(seqs):
        seq2idx[g]= n
    # Array with cluster numbers for each sequence
    group = np.array([seq2idx[i[col]] if not isinstance(i[col], float) \
                      else pd.NA for _,i in DF.iterrows()])
    DF[group_col_name]=group # Add cluster number to df
    return DF
#---------------------------------------------------------------------------------------------------#
def generate_clone_sets(df,cols):
    """ Returns the clones as a list of sets of CDR3 chains.

        Forms the TCR clones from a dataset of T-cells containing CDR3 sequences from
        different loci. Specifically, takes the sequences in the loci given in the list
        'cols' from dataframe 'df' and creates python sets of sequences for each cell. The
        set will be invariant to the order of inclusion of the elements, ignoring the allele
        where the sequence comes from. Each set will be the clonal form of a cell. Loci containing
        NaNs as sequence are ignored.

        IMPORTANT ASSUMPTION: It is assumed that the probability of generating 2 identical
        sequences in different locus alleles for a given cell is zero.

        Parameters
        ----------
        df : pd.DataFrame
            T-cell dataset containing the CDR3 sequences on interest for each cell.
        cols : list of strings
            List with the names of the columns in 'df' where the sequences are. It has to
            have 4 loci, 2 sequences per locus corresponding to the 2 alleles.

        Raises
        ------
        ValueError
            If 'cols' has not exactly 4 entries.

        Returns
        -------
        list of sets
            List containing sets of sequences defining the clonal form of each cell.
    """
    if len(cols)<4:
        raise ValueError("Incomplete column list. It has to have 4 loci.")
    elif len(cols)>4:
        raise ValueError("Too many loci columns given. It has to have 4 loci.")
    seq_set = []
    for i in df.iterrows():
        set_as_list = []
        for s in i[1][cols].values:
            if not pd.isna(s):
                set_as_list.append(s)
        set_entry = set(set_as_list)
        seq_set.append(set(set_entry))
    return seq_set
#---------------------------------------------------------------------------------------------------#
def group_sets(set_list):
    """ Groups identical sets and returns a list with each set's group number.

        Iterates over a list of sets making one-by-one comparisons between elements
        to determine equality. The method uses 2 nested for loops, having a complexity of
        at most N^2, this due to the unhashable nature of python sets. A preliminar if
        statement determining if a set has or has not been evaluated before, reducing
        complexity and improving efficiently *slightly*. Identical sets are grouped
        together and assigned a group number, which is appended to a list which is returned.

        Parameters
        ----------
        set_list : list
            List containing sets to be compared.

        Returns
        -------
        list
            List containing group numbers. They are not necessarly consecutive.
    """
    l = len(set_list)
    group = np.zeros([l],dtype=int)
    for i in range(l):
        # Skips the rest of the loop for already assigned sets
        if group[i]!=0:
            continue
        t = [] # List of equality booleans
        for j in range(l):
            t.append(set_list[i]==set_list[j])
        idx = np.where(t)[0] # Array on indices where there is equality
        if set_list[i]: # For not empty sets
            group[idx]=i
        else: # For empty sets
            group[idx]=-1
    return group
#---------------------------------------------------------------------------------------------------#
def concat_seqs_in_set(set_list):
    """ Concatenates the elements of the sets of a list of sets.

        Iterates over a list of sets 'set_list', concatenating its elements
        provided that they are strings and returns the concatenated sequences as
        entries of a list.

        Parameters
        ----------
        seq_set: list
            List containing sets of sequences defining the clonal form of each cell.

        Raises
        ------
        ValueError
            If an element of any of the sets is not a string.

        Returns
        -------
        list
            List of concatenated strings of the same length as the input list.
            Can have empty string ('') entries.

    """
    clones = []
    for s in set_list:
        clone = ''
        for seq in s:
            if type(seq) is not str:
                raise ValueError("The elements of the set have to be strings")
            clone = clone+seq
        clones.append(clone)
    return clones
