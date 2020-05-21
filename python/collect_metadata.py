#!/usr/bin/env python3
# coding: utf-8

'''Collect GISAID metadata and acknowledgements into one file
Author: Albert Chen (Deverman Lab, Broad Institute)
'''

import pandas as pd

from pathlib import Path

# Resolve any symlinks --> absolute path
data_dir = Path('data').resolve()

# Get all files, we'll sort them in a sec
all_data_files = sorted(data_dir.iterdir())

# ------------------------
# COLLECT PATIENT METADATA
# ------------------------

patient_meta_files = [f for f in all_data_files if 'patient_meta' in f.name]

print('Collecting {} patient metadata files...'.format(len(patient_meta_files)), end='', flush=True)
patient_meta_df = pd.DataFrame()
for f in patient_meta_files:
    _df = pd.read_csv(f, sep='\t', skiprows=2)
    patient_meta_df = pd.concat([patient_meta_df, _df], ignore_index=True)

# Save dataframe
patient_meta_df.to_csv('processed_data/patient_meta.csv', index=False)
print('done', flush=True)

# --------------------------------------
# COLLECT SEQUENCING TECHNOLOGY METADATA
# --------------------------------------

seq_meta_files = [f for f in all_data_files if 'seq_meta' in f.name]

print('Collecting {} sequencing tech metadata files...'.format(len(seq_meta_files)), end='', flush=True)
seq_meta_df = pd.DataFrame()
for f in seq_meta_files:
    _df = pd.read_csv(f, sep='\t', skiprows=2)
    seq_meta_df = pd.concat([seq_meta_df, _df], ignore_index=True)

# Save dataframe
seq_meta_df.to_csv('processed_data/seq_meta.csv', index=False)
print('done', flush=True)

# ------------------------
# COLLECT ACKNOWLEDGEMENTS
# ------------------------

acknowledgement_files = [f for f in all_data_files if 'ack' in f.name]

print('Collecting {} acknowledgement files...'.format(len(acknowledgement_files)), end='', flush=True)
acknowledgement_df = pd.DataFrame()
for f in acknowledgement_files:
    _df = pd.read_excel(f, skiprows=[0, 1, 3])
    acknowledgement_df = pd.concat([acknowledgement_df, _df], ignore_index=True)

# Save dataframe
acknowledgement_df.to_csv('processed_data/acknowledgements.csv', index=False)
print('done', flush=True)

print('Saving unique acknowledgements...', end='', flush=True)
# Save unique acknolwedgements
unique_ack_df = acknowledgement_df.copy()
# Drop sequence information columns
unique_ack_df.drop(columns=['Accession ID', 'Virus name', 'Collection date'], inplace=True)
# Drop duplicates
unique_ack_df.drop_duplicates(inplace=True)
# Save to disk
unique_ack_df.to_csv('processed_data/unique_ack.csv', index=False)
print('done', flush=True)