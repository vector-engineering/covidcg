#!/usr/bin/env python3
# coding: utf-8

'''Process clade information
Author: Albert Chen (Deverman Lab, Broad Institute)
'''

import pandas as pd
import numpy as np
import re

from pathlib import Path

from fasta import read_fasta_file
from util import translate

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
    dna_snp_df['gisaid_id'] = dna_snp_df['taxon'].str.split('|', expand=True)[1]
    dna_snp_df = dna_snp_df.reset_index(drop=True)
    print('done. Loaded {} DNA SNPs'.format(len(dna_snp_df)), flush=True)
    # dna_snp_df.to_csv('dna_snp.csv', index=False)

    return dna_snp_df


def load_aa_snps():
    pass
    # Load all AA SNP files
    # aa_snp_files = sorted((data_dir / 'aa_snp').glob('*.csv'))
    # print('Loading {} AA SNP files...'.format(len(aa_snp_files)), end='', flush=True)
    # # Load into dataframe
    # aa_snp_df = pd.DataFrame()
    # for f in aa_snp_files:
    #     aa_snp_df = pd.concat([aa_snp_df, pd.read_csv(f)], ignore_index=True)
    # # Drop index column
    # aa_snp_df.drop(columns=['index'], inplace=True)
    # # Extract the GISAID ID
    # aa_snp_df['gisaid_id'] = aa_snp_df['taxon'].str.split('|', expand=True)[1]
    # print('done. Loaded {} DNA SNPs'.format(len(aa_snp_df)), flush=True)


def load_taxon_lineages():
    # Load lineage files
    lineage_files = sorted((data_dir / 'lineage_meta').glob('*.csv'))
    print('Loading {} lineage metadata files...'.format(len(lineage_files)), end='', flush=True)
    lineage_df = pd.DataFrame()
    for f in lineage_files:
        lineage_df = pd.concat([lineage_df, pd.read_csv(f)], ignore_index=True)
    unique_lineages = sorted(lineage_df['lineage'].unique())
    print('done. Loaded {} unique lineages'.format(len(unique_lineages)), flush=True)

    # Convert sample_date to datetime
    lineage_df['sample_date'] = pd.to_datetime(lineage_df['sample_date'], yearfirst=True)

    # Save combined lineages file
    all_lineages_path = data_dir / 'all_lineages.csv'
    lineage_df.to_csv(all_lineages_path, index=False)
    print('Saved all lineage assignments to {}'.format(all_lineages_path))

    return lineage_df


def get_consensus_dna_snps(dna_snp_df, lineage_df):

    # Fraction of taxons that need to have a SNP for it to be considered a consensus
    # SNP for a lineage
    # TODO: make this a CLI arg
    consensus_fraction = 0.9

    # Store lineage - SNP associations in here
    lineage_snp_df = pd.DataFrame()
    
    for lineage in unique_lineages:
        # Get all taxon GISAID IDs assigned to this lineage
        lin_df = lineage_df.loc[lineage_df['lineage'] == lineage, :].reset_index(drop=True)

        # Count SNPs per GISAID ID, for this lineage
        snp_group_df = (
            dna_snp_df.loc[dna_snp_df['gisaid_id'].isin(lin_df['gisaid_id']), :]
            .groupby(['pos', 'ref', 'alt'], as_index=False)['gisaid_id']
            .count()
            # Sort by frequency, descending order
            .sort_values('gisaid_id', ascending=False)
        )

        # SNPs from Clade A:
        consensus_snp_df = snp_group_df.loc[
            snp_group_df['gisaid_id'] > (len(lin_df) * consensus_fraction), :
        ].reset_index(drop=True)
        consensus_snp_df['lineage'] = lineage
        # Order by position
        consensus_snp_df = consensus_snp_df.sort_values('pos')

        # Append to master dataframe
        lineage_snp_df = pd.concat(
            [
                lineage_snp_df,
                consensus_snp_df[['lineage', 'pos', 'ref', 'alt']]
            ],
            ignore_index=True
        )

    # Save to disk
    print('Saving lineage - SNP definitions'.format(str(lineage_snp_out_path)))
    # print(lineage_snp_df)
    lineage_snp_df.to_csv(data_dir / 'lineage_snp.csv', index=False)
    lineage_snp_df.to_json(data_dir / 'lineage_snp.json', orient='records')
    

    return lineage_snp_df


def main():

    dna_snp_df = load_dna_snps()
    lineage_df = load_taxon_lineages()

    get_consensus_dna_snps(dna_snp_df, lineage_df)


if __name__ == '__main__':
    main()