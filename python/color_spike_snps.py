#!/usr/bin/env python3
# coding: utf-8

'''Color Spike AAs by mutation/SNP frequency

In PyMOL, first navigate to this python/ folder using cd
Then run this script with 'run color_spike_snps.py'

Author: Albert Chen (Deverman Lab - Broad Institute)
'''

import os
import sys

from pathlib import Path
from pymol import cmd, stored

sys.path.append(os.getcwd())

pdb_name = '6X2A'

project_root_path = Path(os.getcwd()).resolve().parent
data_dir = (project_root_path / 'data').resolve() # Resolve any symlinks --> absolute path

# Load spike color data
# We'll have to do it without pandas
# Get it into a dict of pos -> hex color
spike_pos_colors = {}

with (data_dir / 'spike_snp_freq_by_pos.csv').open('r') as fp:
    # Columns: pos, count, count_log, count_norm, count_log_norm, color_hex, color_log_hex
    i = 0
    for line in fp:
        i += 1
        # Skip header
        if i == 1:
            continue

        # Remove newline char
        line = line.rstrip()

        # Split by comma
        line_split = line.split(',')

        spike_pos_colors[line_split[0]] = line_split[6]

# print(spike_pos_colors)

def load_spike_structure():
    cmd.delete('all')

    cmd.set('fetch_path', 'pdbs')
    cmd.fetch(pdb_name, name=pdb_name)

    # Get rid of water/hydrogens in the structure
    cmd.remove('solvent')
    cmd.remove('hydrogens')

    # Clean up the rest... remove all non-protein atoms (PyMOL 2.1+ only)
    cmd.remove('(not polymer.protein)')

    # Hide cartoon and show surface
    cmd.hide('cartoon')
    cmd.show('surface')

def color_residues():
    for resi, color in spike_pos_colors.items():
        cmd.color(color, 'resi {}'.format(resi))

load_spike_structure()
color_residues()
