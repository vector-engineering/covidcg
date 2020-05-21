#!/usr/bin/env python3
# coding: utf-8

'''Assign sequences to clades using pangolin
Author: Albert Chen (Deverman Lab, Broad Institute)
'''

import argparse
import os
import re
import subprocess

from pathlib import Path

# TODO: add argument for skipping existing lineage_report
# note, it's probably a good idea to re-run all when a new lineage report comes out
# just in case old sequences get re-assigned.
# TODO: add threads argument to pass onto pangolin

# We're using all cores. gotta go fast
print('Running with {} CPUs'.format(os.cpu_count()))

# Get all fasta files from data/
data_path = Path('data')
fasta_files = sorted(data_path.glob('*.fasta'))

# TODO: check for existing

print('Found {} fasta files: \n  - {}'.format(len(fasta_files), '\n  - '.join([f.name for f in fasta_files])))

# Run pangolin
for i, ff in enumerate(fasta_files):
    print('Processing {} / {}: {}'.format(i + 1, len(fasta_files), ff.name))
        
    # Output file is the same as input, but chop off "gisaid" from the beginning
    # And it'll start with "lineage_report"
    outfile = 'lineage_report_' + re.sub('gisaid_', '', ff.stem) + '.csv'

    subprocess.run([
        'pangolin', str(ff),
        '-o', 'processed_data',
        '--outfile', outfile,
        '--threads', str(os.cpu_count())
    ])

print('Done')
