#!/usr/bin/env python3
# coding: utf-8

'''Merge interactions from pymol_interacting_nab.py into spike_structures.csv

Author: Albert Chen (Deverman Lab - Broad Institute)
'''

import pandas as pd

from pathlib import Path

project_root_path = Path(__file__).resolve().parent.parent
data_dir = (project_root_path / 'data').resolve() # Resolve any symlinks --> absolute path

df = pd.read_csv(data_dir / 'spike_structures.csv')
interact_df = pd.read_csv(data_dir / 'spike_interactions.csv')

df['interacting_residues'] = df['PDB ID'].map(
    pd.Series(
        interact_df['interacting_residues'].values,
        index=interact_df['pdb_id']
    )
)

df.to_csv(data_dir / 'spike_structures2.csv', index=False)
