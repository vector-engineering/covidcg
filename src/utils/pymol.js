import { reds } from '../constants/colors';

const numColors = reds.length;

export function mutationHeatmapToPymolScript({ pdbId, snvs }) {
  let script = `#!/usr/bin/env python3
# coding: utf-8

## AUTO-GENERATED BY COVIDCG.ORG

from pymol import cmd, stored

pdb_name = '${pdbId}'

cmd.delete('all')
cmd.fetch(pdb_name, name=pdb_name)

# Get rid of water/hydrogens in the structure
cmd.remove('solvent')
cmd.remove('hydrogens')
# Clean up the rest... remove all non-protein atoms (PyMOL 2.1+ only)
cmd.remove('(not polymer.protein)')

# Hide cartoon and show surface
cmd.hide('cartoon')
cmd.show('surface')

# Default color
cmd.color('white')

# RESI COLORING CODE
`;

  snvs.forEach((snv) => {
    const colorInd = Math.floor((snv.fraction - 0.001) * numColors);
    // PyMOL needs colors in "0xRRGGBB" instead of "#RRGGBB"
    script += `cmd.color('0x${reds[colorInd].substr(1)}', 'resi ${snv.pos}')\n`;
  });

  return script;
}
