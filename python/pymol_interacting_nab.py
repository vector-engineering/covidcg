# coding: utf-8
#!/usr/bin/env python3

'''Get residues in Spike that interface with Nabs or other antibodies/nanobodies

In PyMOL, first navigate to this python/ folder using cd
Then run this script with 'run pymol_interacting_nab.py'

Author: Albert Chen (Deverman Lab - Broad Institute)
'''

import os
import sys

from pathlib import Path
from pymol import cmd, stored

sys.path.append(os.getcwd())

project_root_path = Path(os.getcwd()).resolve().parent
data_dir = (project_root_path / 'data').resolve() # Resolve any symlinks --> absolute path

cmd.set('fetch_path', 'pdbs')

# Load spike_structures.csv
structures = {} # structure -> dict
with (data_dir / 'spike_structures.csv').open('r') as fp:
    # PDB ID, PDB URL, Name, Resolution, Method, Binder Name, Monomer Chain, Interacting Chain, Notes
    i = 0
    for line in fp:
        i += 1
        # Skip header
        if i == 1:
            continue

        # Strip newline at end
        line = line.rstrip()
        # Split into chunks by comma
        line = line.split(',')
        
        # PDB ID
        pdb_id = line[0]
        structures[pdb_id] = {}
        
        # Fill in fields
        structures[pdb_id]['name'] = line[1]
        structures[pdb_id]['binder'] = line[5]
        structures[pdb_id]['monomer_chain'] = line[6]
        structures[pdb_id]['binder_chain'] = line[7].split(';')

print(structures)


# Interaction distance in angstroms
interaction_distance = 8.0

# dict of model -> list of residue indexes
interacting_aas = {}

for structure, vals in structures.items():
    # Skip if there's no binder
    if not vals['binder'].strip():
        continue

    print(structure, vals['name'], vals['binder'])

    # Clear workspace and get structure from PDB
    cmd.delete('all')
    cmd.fetch(structure, name=structure)

    # Get rid of water/hydrogens in the structure
    cmd.remove('solvent')
    cmd.remove('hydrogens')

    # Clean up the rest... remove all non-protein atoms (PyMOL 2.1+ only)
    cmd.remove('(not polymer.protein)')

    
    # Get all atoms in chain A within [distance] of chain H or chain L. 
    # Expand to entire residue, and then get just the alpha carbons
    cmd.select('interacting', '((byres (chain {}) within {} of ({})) and name ca)'.format(
        vals['monomer_chain'],
        interaction_distance,
        ' or '.join(['chain ' + x for x in vals['binder_chain']])
    ))
    # Store residue indexes
    stored.interacting = []
    cmd.iterate('interacting', 'stored.interacting.append(resi)')

    # Store in master dict
    interacting_aas[structure] = stored.interacting

# Clean up
cmd.delete('all')

# ----------------------------------------------------------------------------------
# Post-processing
# ----------------------------------------------------------------------------------

print(interacting_aas)

# Write to CSV
with (data_dir / 'spike_interactions.csv').open('w') as fp:
    fp.write('pdb_id,interacting_residues\n')
    for pdb_id, contacts in interacting_aas.items():
        fp.write(pdb_id + ',')
        fp.write(';'.join(contacts))
        fp.write('\n')
