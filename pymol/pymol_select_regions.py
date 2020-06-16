#!/usr/bin/env python3
# coding: utf-8

'''Select regions on Spike

Author: Albert Chen (Deverman Lab - Broad Institute)
'''

import csv

from pathlib import Path
from pymol import cmd, stored

sys.path.append(os.getcwd())

project_root_path = Path(os.getcwd()).resolve().parent
data_dir = (project_root_path / 'data').resolve() # Resolve any symlinks --> absolute path


# Read in spike_structures.csv into a dict
with (data_dir / 'spike_structures.csv').open('r') as fp:
    reader = csv.reader(fp, delimiter=',', quotechar='"')
    i = 0
    col_names = []
    struct_rows = {}
    for row in reader:
        i += 1
        # Process header separately
        # Save the first column (PDB ID) as the index key
        if i == 1:
            col_names = row[1:]
            continue
        
        # Store values in dict
        struct_rows[row[0]] = {}
        for j, col in enumerate(col_names):
            struct_rows[row[0]][col] = row[j + 1]

# print(struct_rows)

sele = 'spike'

# For the ACE2 RBD, using contacts from 6M0J
cmd.select('rbd', '{} and (resi 417 or resi 446 or resi 449 or resi 453 or resi 455 or resi 456 or resi 475 or resi 486 or resi 487 or resi 489 or resi 493 or resi 493 or resi 496 or resi 498 or resi 500 or resi 501 or resi 502 or resi 505)'.format(sele))

# Get CB6 Nab contacts (7C01)
cmd.select('cb6', '{} and ({})'.format(
    sele, ' or '.join(['resi ' + x for x in struct_rows['7C01']['interacting_residues'].split(';')])
))

# Get B38 Nab contacts (7BZ5)
cmd.select('b38', '{} and ({})'.format(
   sele, ' or '.join(['resi ' + x for x in struct_rows['7BZ5']['interacting_residues'].split(';')])
))

# Get S309 Nab contacts (6WPT)
cmd.select('s309', '{} and ({})'.format(
    sele, ' or '.join(['resi ' + x for x in struct_rows['6WPT']['interacting_residues'].split(';')])
))

# Get 7BWJ contacts (7BWJ)
cmd.select('7bwj', '{} and ({})'.format(
    sele, ' or '.join(['resi ' + x for x in struct_rows['7BWJ']['interacting_residues'].split(';')])
))

# CR3022 contacts (6W41)
cmd.select('cr3022', '{} and ({})'.format(
    sele, ' or '.join(['resi ' + x for x in struct_rows['6W41']['interacting_residues'].split(';')])
))

# H11-D4 contacts (6Z43)
cmd.select('h11d4', '{} and ({})'.format(
    sele, ' or '.join(['resi ' + x for x in struct_rows['6Z43']['interacting_residues'].split(';')])
))