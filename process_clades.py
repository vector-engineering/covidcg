#!/usr/bin/env python3
# coding: utf-8

'''Process clade information
Author: Albert Chen (Deverman Lab, Broad Institute)
'''

import pandas as pd
import numpy as np
import networkx as nx
import re

from collections import defaultdict

from fasta import read_fasta_file

def load_clades():
    # Load clades
    clades_df = pd.read_csv('clades.tsv', sep='\t')
    # print(clades_df)
    return clades_df


def build_clade_tree(clades_df):
    '''Build the clade tree
    '''

    clades = clades_df['Clade'].unique()
    # print('Clades: {}'.format(', '.join(clades)))

    # Group by name length
    clade_groups = defaultdict(list)
    for clade in clades:
        clade_groups[len(clade)].append(clade)
    clade_group_levels = sorted(clade_groups.keys())
    # print(clade_groups)

    # Build tree, level by level
    clade_tree = nx.DiGraph()
    # Add the root level first
    clade_tree.add_node('root')

    for i, k in enumerate(clade_group_levels):
        for clade in clade_groups[k]:
            # Find the ancestors in levels above
            # If no ancestor is found in the level directly above,
            # then keep going up until it's found
            ancestor = None
            cur_level = i - 1
            while ancestor is None and cur_level >= 0:
                ancestor = [
                    c for c in clade_groups[clade_group_levels[cur_level]] 
                    if re.search(c, clade)
                ]
                if ancestor:
                    ancestor = ancestor[0]
                else:
                    ancestor = None
                    cur_level -= 1

            # If no ancestors are found still, assign to the root node
            # This will be true of all the first-level nodes
            if ancestor is None:
                ancestor = 'root'
            
            # print(clade, ancestor)
            # Add node
            clade_tree.add_node(clade)
            # Add edge
            clade_tree.add_edge(clade, ancestor)

    # print(clade_tree.nodes)
    # print(clade_tree.edges)

    #print(clade_tree.edges('B2'))
    #print(clade_tree.edges('B'))

    return clade_tree


def build_clade_sequences(clades_df, clade_tree):
    '''Build sequences for each clade
    '''

    # Load the reference sequence
    with open('reference.fasta', 'r') as fp:
        lines = fp.readlines()
        ref = read_fasta_file(lines)
        ref_seq = list(ref.values())[0]
    # print(len(ref_seq))

    clade_seqs = {
        'root': ref_seq
    }
    # Track positions that are changing
    changing_positions = set()
    for clade in clade_tree.nodes: # These should be sorted by length
        # Skip the root node
        if clade == 'root':
            continue
        
        # Use the ancestor's sequence, if there is an ancestor.
        # If not, use the reference sequence
        ancestor = list(clade_tree.edges(clade))[0][1]
        seq = clade_seqs[ancestor]
        # print(clade, ancestor)

        # Get the mutations from the original clades_df
        muts = clades_df.loc[clades_df['Clade'] == clade, :].reset_index(drop=True)
        # print(muts)
        # Modify the sequence
        for i in range(len(muts)):
            pos = muts.at[i, 'Pos']
            ref = muts.at[i, 'Ref']
            alt = muts.at[i, 'Alt']

            # Add position to the master list
            changing_positions.add(pos)

            # Check to see if the clade reference matches the loaded reference
            # print(ref_seq[pos - 1], ref)

            # List operations as easier
            seq = list(seq)
            # If it's a SNP, then just replace the base
            if len(ref) == 1 and len(alt) == 1:
                seq[pos - 1] = alt
            # If it's an insertion
            elif len(ref) == 1 and len(alt) > 1:
                # Replace the first base
                seq[pos - 1] = alt[0]
                # Insert bases
                seq.insert(pos, alt[1:])
            # If it's a deletion
            elif len(ref) > 1 and len(alt) == 1:
                # Replace the first base
                seq[pos - 1] = alt
                # Pop bases
                for j in range(len(ref) - 1):
                    seq.pop(pos)
            # Back to a string
            seq = ''.join(seq)

        clade_seqs[clade] = seq

    # print(clade_seqs)
    # print(changing_positions)

    return clade_seqs, changing_positions


def main():
    clades_df = load_clades()
    clade_tree = build_clade_tree(clades_df)
    clade_seqs, changing_positions = build_clade_sequences(clades_df, clade_tree)


if __name__ == '__main__':
    main()