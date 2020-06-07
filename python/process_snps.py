#!/usr/bin/env python3
# coding: utf-8

'''Load SNP data, create SNP signatures

Author: Albert Chen - Deverman Lab, Broad Institute
'''

import numpy as np
import pandas as pd

from pathlib import Path

project_root_path = Path(__file__).resolve().parent.parent
data_dir = (project_root_path / 'data').resolve() # Resolve any symlinks --> absolute path

def load_dna_snps():
    # Load all DNA SNP files
    dna_snp_files = sorted((data_dir / 'dna_snp').glob('*.csv'))
    print('Loading {} DNA SNP files...'.format(len(dna_snp_files)), end='', flush=True)
    # Load into dataframe
    dna_snp_df = pd.DataFrame()
    for f in dna_snp_files:
        dna_snp_df = pd.concat([dna_snp_df, pd.read_csv(f)], ignore_index=True)
    # Extract the GISAID ID
    # dna_snp_df['gisaid_id'] = dna_snp_df['taxon'].str.split('|', expand=True).apply(lambda x: x[1])

    # Drop duplicate entries
    dna_snp_df.drop_duplicates(['taxon', 'pos'], inplace=True)

    dna_snp_df = dna_snp_df.reset_index(drop=True)
    print('done. Loaded {} DNA SNPs'.format(len(dna_snp_df)), flush=True)
    # dna_snp_df.to_csv('dna_snp.csv', index=False)

    return dna_snp_df


def load_aa_snps():
    # Load all AA SNP files
    aa_snp_files = sorted((data_dir / 'aa_snp').glob('*.csv'))
    print('Loading {} AA SNP files...'.format(len(aa_snp_files)), end='', flush=True)
    # Load into dataframe
    aa_snp_df = pd.DataFrame()
    for f in aa_snp_files:
        aa_snp_df = pd.concat([aa_snp_df, pd.read_csv(f)], ignore_index=True)
    # Extract the GISAID ID
    # aa_snp_df['gisaid_id'] = aa_snp_df['taxon'].str.split('|', expand=True).apply(lambda x: x[1])

    # Drop duplicate entries
    aa_snp_df.drop_duplicates(['taxon', 'gene', 'pos'], inplace=True)

    aa_snp_df = aa_snp_df.reset_index(drop=True)

    # Position 0-indexed -> 1-indexed
    aa_snp_df['pos'] = aa_snp_df['pos'] + 1

    print('done. Loaded {} AA SNPs'.format(len(aa_snp_df)), flush=True)

    return aa_snp_df


def generate_snp_signatures(dna_snp_df, aa_snp_df):
    # SNPs must occur at least this many times to pass filters
    count_threshold = 10
    # SNP signatures must be present in at least this many sequences to pass filters
    sig_count_threshold = 10

    # Collapse by taxon and count occurrences
    print('Collapsing by taxon and counting occurrences...', end='', flush=True)
    aa_snp_count_df = (
        aa_snp_df
        .groupby(['gene', 'pos', 'ref', 'alt'], as_index=False)
        .count()
        .rename(columns={'taxon': 'count'})
    )
    dna_snp_count_df = (
        dna_snp_df
        .groupby(['pos', 'ref', 'alt'], as_index=False)
        .count()
        .rename(columns={'taxon': 'count'})
    )
    print('done')

    # Filter out SNPs
    print('Filtering low-occurrence SNPs...', end='', flush=True)
    valid_aa_snps = (
        aa_snp_count_df
        .loc[aa_snp_count_df['count'] >= count_threshold, :]
        .reset_index(drop=True)
        .sort_values('count', ascending=False)
    )
    valid_dna_snps = (
        dna_snp_count_df
        .loc[dna_snp_count_df['count'] >= count_threshold, :]
        .reset_index(drop=True)
        .sort_values('count', ascending=False)
    )
    print('done')

    # Create unique SNP string
    print('Creating unique SNP strings...', end='', flush=True)
    valid_aa_snps['snp_str'] = (
        valid_aa_snps
        [['gene', 'pos', 'ref', 'alt']]
        .applymap(str)
        .agg('|'.join, axis=1)
    )
    valid_dna_snps['snp_str'] = (
        valid_dna_snps
        [['pos', 'ref', 'alt']]
        .applymap(str)
        .agg('|'.join, axis=1)
    )
    print('done')

    # Save to disk
    #valid_dna_snps.to_csv(data_dir / 'valid_dna_snps.csv', index=False)
    #valid_aa_snps.to_csv(data_dir / 'valid_aa_snps.csv', index=False)

    # Generate SNP strings for initial dataframes
    print('Generating SNP strings for initial dataframes...', end='', flush=True)
    aa_snp_df['snp_str'] = aa_snp_df['gene'].str.cat([
        aa_snp_df['pos'].astype(str), 
        aa_snp_df['ref'],
        aa_snp_df['alt']
    ], sep='|')
    dna_snp_df['snp_str'] = dna_snp_df['pos'].astype(str).str.cat([
        dna_snp_df['ref'],
        dna_snp_df['alt']
    ], sep='|')
    print('done')

    # Filter SNPs by valid SNPs
    print('Generating SNP signatures...', end='', flush=True)
    aa_snp_df = aa_snp_df.loc[aa_snp_df['snp_str'].isin(valid_aa_snps['snp_str']), :].reset_index(drop=True)
    dna_snp_df = dna_snp_df.loc[dna_snp_df['snp_str'].isin(valid_dna_snps['snp_str']), :].reset_index(drop=True)

    # Group by taxon and make a ';' delimited list of snp_strs
    aa_snp_group_df = aa_snp_df.groupby('taxon')['snp_str'].agg(';'.join).reset_index()
    dna_snp_group_df = dna_snp_df.groupby('taxon')['snp_str'].agg(';'.join).reset_index()

    # Count occurences of SNP signatures
    aa_snp_sig_count_df = aa_snp_group_df.groupby('snp_str').count().sort_values('taxon', ascending=False)
    dna_snp_sig_count_df = dna_snp_group_df.groupby('snp_str').count().sort_values('taxon', ascending=False)

    # Only take SNP signatures with >= N sequences
    valid_aa_snp_sigs = aa_snp_sig_count_df.index.values[aa_snp_sig_count_df['taxon'] >= sig_count_threshold]
    valid_dna_snp_sigs = dna_snp_sig_count_df.index.values[dna_snp_sig_count_df['taxon'] >= sig_count_threshold]

    # Expand string form for easier searching
    valid_aa_snp_sigs_explode = [set(s.split(';')) for s in valid_aa_snp_sigs]
    valid_dna_snp_sigs_explode = [set(s.split(';')) for s in valid_dna_snp_sigs]
    print('done')

    print('Assigning AA SNP signatures...', end='', flush=True)
    # Assign taxons that don't already fit perfectly into a group, to a group
    aa_snp_group_df['snp_sig'] = None
    dna_snp_group_df['snp_sig'] = None

    for i, row in aa_snp_group_df.iterrows():
        # Testing
        #if i > 1:
        #    break
        
        # Check to see if its snp_str is already valid
        if row['snp_str'] in valid_aa_snp_sigs:
            # print('Exact match')
            row['snp_sig'] = row['snp_str']
            continue
        
        # print('No match')
        
        # Otherwise compute the "hamming distance" between this snp_str and
        # all valid snp_strs
        scores = [len(a & set(row['snp_str'].split(';'))) for a in valid_aa_snp_sigs_explode]
        # And take the snp_sig with the largest score.
        # TODO: This isn't perfect. since argmax() takes the FIRST largest value if the largest
        # value occurs more than once. idk what the best way to do this is... if we see
        # a lot of split assignments then we should really merge the signatures that they're split
        # between.
        row['snp_sig'] = valid_aa_snp_sigs[np.argmax(scores)]
    print('done')
        
    print('Assigning DNA SNP signatures...', end='', flush=True)
    for i, row in dna_snp_group_df.iterrows():
        if row['snp_str'] in valid_dna_snp_sigs:
            row['snp_sig'] = row['snp_str']
            continue
        scores = [len(a & set(row['snp_str'].split(';'))) for a in valid_dna_snp_sigs_explode]
        row['snp_sig'] = valid_dna_snp_sigs[np.argmax(scores)]
    print('done')

    #print(aa_snp_group_df)
    #print(dna_snp_group_df)

    # Save to disk
    dna_snp_group_df.to_csv(data_dir / 'dna_snp_group.csv', index=False)
    aa_snp_group_df.to_csv(data_dir / 'aa_snp_group.csv', index=False)

    return dna_snp_group_df, aa_snp_group_df


def main():
    dna_snp_df = load_dna_snps()
    aa_snp_df = load_aa_snps()

    # print(dna_snp_df)

    dna_snp_group_df, aa_snp_group_df = generate_snp_signatures(dna_snp_df, aa_snp_df)

if __name__ == '__main__':
    main()
