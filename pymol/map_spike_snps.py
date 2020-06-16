#!/usr/bin/env python3
# coding: utf-8

'''Map Spike AA SNPs onto the trimer structure

Author: Albert Chen (Deverman Lab - Broad Institute)
'''

import numpy as np
import pandas as pd
import sys

from matplotlib import cm
from pathlib import Path

project_root_path = Path(__file__).resolve().parent.parent
data_dir = (project_root_path / 'data').resolve() # Resolve any symlinks --> absolute path

# First we need to aggregate all Spike AA SNPs by position
aa_snp_files = sorted((data_dir / 'aa_snp').glob('*.csv'))
# print(aa_snp_files)

aa_snp_df = pd.DataFrame()
for aa_snp_file in aa_snp_files:
    aa_snp_df = pd.concat([aa_snp_df, pd.read_csv(aa_snp_file)], ignore_index=True)

aa_snp_pos_df = (
    aa_snp_df
    .loc[aa_snp_df['gene'] == 'S', :]
    .groupby('pos')
    ['taxon']
    .count()
)

# Fill in missing positions (1-indexed) with 0
count_dict = {}
for i in range(1, 1275):
    count_dict[i] = 0

for pos, count in aa_snp_pos_df.iteritems():
    count_dict[pos + 1] = count

# Set D614G to 0. to make the colormap better
count_dict[614] = 0

# Re-cast to dataframe
aa_snp_pos_df = pd.DataFrame.from_dict(count_dict, orient='index', columns=['count'])

# Log-transform counts
aa_snp_pos_df['count_log'] = np.log10(aa_snp_pos_df['count'])

# Normalize counts by max
min_val = 0.1
max_val = 1.0
aa_snp_pos_df['count_norm'] = min_val + ((aa_snp_pos_df['count'] / aa_snp_pos_df['count'].max()) / (max_val - min_val))
aa_snp_pos_df['count_log_norm'] = min_val + ((aa_snp_pos_df['count_log'] / aa_snp_pos_df['count_log'].max()) / (max_val - min_val))

# Get colormap
cmap = cm.get_cmap('autumn')

# Map normalized counts to colormap
aa_snp_pos_df['color_rgb'] = aa_snp_pos_df['count_norm'].apply(cmap)
aa_snp_pos_df['color_log_rgb'] = aa_snp_pos_df['count_log_norm'].apply(cmap)

# Set 0 or neg inf counts to gray
for i, row in aa_snp_pos_df.iterrows():
    if row['count'] == 0 or row['count'] == 1:
        aa_snp_pos_df.loc[i, 'color_rgb'] = pd.Series([1]).apply(lambda x: (0.7, 0.7, 0.7, 1.0)).values
    if np.isinf(row['count_log']) or row['count_log'] == 0:
        aa_snp_pos_df.loc[i, 'color_log_rgb'] = pd.Series([1]).apply(lambda x: (0.7, 0.7, 0.7, 1.0)).values

# Convert to hex
def rgb_to_hex(rgb):
    out = 0
    # Red channel
    out += (int(rgb[0] * 255) * (256 ** 2))
    # Blue channel
    out += (int(rgb[1] * 255) * (256))
    # Green channel
    out += (int(rgb[2] * 255))
    return hex(out)
aa_snp_pos_df['color_hex'] = aa_snp_pos_df['color_rgb'].apply(rgb_to_hex)
aa_snp_pos_df['color_log_hex'] = aa_snp_pos_df['color_log_rgb'].apply(rgb_to_hex)

# Drop RGB columns
aa_snp_pos_df.drop(columns=['color_rgb', 'color_log_rgb'], inplace=True)

print(aa_snp_pos_df)

# Write to csv
aa_snp_pos_df.to_csv(data_dir / 'spike_snp_freq_by_pos.csv', index_label='pos')


