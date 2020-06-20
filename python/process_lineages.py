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
from process_snps import load_dna_snps, load_aa_snps
from util import translate

project_root_path = Path(__file__).resolve().parent.parent
data_dir = (project_root_path / 'data').resolve() # Resolve any symlinks --> absolute path

def load_taxon_lineages():
    # Load lineage files
    meta_files = sorted((data_dir / 'patient_meta').glob('*.tsv'))
    print('Loading {} metadata files...'.format(len(meta_files)), end='', flush=True)
    lineage_df = pd.DataFrame()
    for f in meta_files:
        lineage_df = pd.concat([lineage_df, pd.read_csv(f, sep='\t', skiprows=2)], ignore_index=True)
    lineage_df = lineage_df.reset_index(drop=True)

    lineage_df.rename(columns={
        'Lineage': 'lineage',
        'Collection date': 'sample_date',
        'Accession ID': 'gisaid_id'
    }, inplace=True)

    # Filter out taxons with any lineage
    lineage_df = lineage_df.loc[lineage_df['lineage'].apply(lambda x: (x is not None) and (type(x) is str)), :].reset_index(drop=True)

    unique_lineages = sorted(lineage_df['lineage'].unique())
    print('done. Loaded {} unique lineages'.format(len(unique_lineages)), flush=True)

    # Convert sample_date to datetime
    lineage_df['sample_date'] = pd.to_datetime(lineage_df['sample_date'], yearfirst=True)

    # Save combined lineages file
    all_lineages_path = data_dir / 'all_lineages.csv'
    lineage_df.to_csv(all_lineages_path, index=False)
    # print('Saved all lineage assignments to {}'.format(all_lineages_path))

    return lineage_df


def get_consensus_dna_snps(dna_snp_df, lineage_df):

    # Fraction of taxons that need to have a SNP for it to be considered a consensus
    # SNP for a lineage
    # TODO: make this a CLI arg
    consensus_fraction = 0.9

    # Store lineage - SNP associations in here
    lineage_snp_df = pd.DataFrame()
    
    unique_lineages = sorted(lineage_df['lineage'].unique())
    for i, lineage in enumerate(unique_lineages):
        # print(i)
        # Get all taxon GISAID IDs assigned to this lineage
        lin_df = lineage_df.copy()
        # print('copied')
        lin_df = lin_df.loc[lineage_df['lineage'] == lineage, :].reset_index(drop=True)
        # print('subsetted')

        # Count SNPs per GISAID ID, for this lineage
        snp_group_df = (
            dna_snp_df.loc[dna_snp_df['gisaid_id'].isin(lin_df['gisaid_id']), :]
            .groupby(['pos', 'ref', 'alt'], as_index=False)['gisaid_id']
            .count()
            # Sort by frequency, descending order
            .sort_values('gisaid_id', ascending=False)
        )
        # Calculate consensus percentage
        snp_group_df['consensus'] = snp_group_df['gisaid_id'] / len(lin_df)

        # print('counted')

        # SNPs from Clade A:
        consensus_snp_df = snp_group_df.loc[
            snp_group_df['consensus'] >  consensus_fraction, :
        ].reset_index(drop=True)
        consensus_snp_df['lineage'] = lineage
        # Order by position
        consensus_snp_df = consensus_snp_df.sort_values('pos')

        # Append to master dataframe
        lineage_snp_df = pd.concat(
            [   
                lineage_snp_df,
                consensus_snp_df[['lineage', 'pos', 'ref', 'alt', 'consensus']]
            ],
            ignore_index=True
        )

    # print(lineage_snp_df.to_csv())

    # Save to disk
    print('Saving lineage - DNA SNP definitions')
    # print(lineage_snp_df)
    lineage_snp_df.to_csv(data_dir / 'lineage_dna_snp.csv', index=False)
    lineage_snp_df.to_json(data_dir / 'lineage_dna_snp.json', orient='records')
    
    return lineage_snp_df


def get_consensus_aa_snps(aa_snp_df, lineage_df):

    # Fraction of taxons that need to have a SNP for it to be considered a consensus
    # SNP for a lineage
    # TODO: make this a CLI arg
    consensus_fraction = 0.9

    # Store lineage - SNP associations in here
    lineage_snp_df = pd.DataFrame()
    
    unique_lineages = sorted(lineage_df['lineage'].unique())
    for lineage in unique_lineages:
        # Get all taxon GISAID IDs assigned to this lineage
        lin_df = lineage_df.copy().loc[lineage_df['lineage'] == lineage, :].reset_index(drop=True)

        # Count SNPs per GISAID ID, for this lineage
        snp_group_df = (
            aa_snp_df.loc[aa_snp_df['gisaid_id'].isin(lin_df['gisaid_id']), :]
            .groupby(['gene', 'pos', 'ref', 'alt'], as_index=False)['gisaid_id']
            .count()
            # Sort by frequency, descending order
            .sort_values('gisaid_id', ascending=False)
        )
        # Calculate consensus percentage
        snp_group_df['consensus'] = snp_group_df['gisaid_id'] / len(lin_df)

        # SNPs from Clade A:
        consensus_snp_df = snp_group_df.loc[
            snp_group_df['consensus'] > consensus_fraction, :
        ].reset_index(drop=True)
        consensus_snp_df['lineage'] = lineage
        # Order by position
        consensus_snp_df = consensus_snp_df.sort_values('pos')

        # Append to master dataframe
        lineage_snp_df = pd.concat(
            [
                lineage_snp_df,
                consensus_snp_df[['lineage', 'gene', 'pos', 'ref', 'alt', 'consensus']]
            ],
            ignore_index=True
        )

    # Save to disk
    print('Saving lineage - AA SNP definitions')
    # print(lineage_snp_df)
    lineage_snp_df.to_csv(data_dir / 'lineage_aa_snp.csv', index=False)
    lineage_snp_df.to_json(data_dir / 'lineage_aa_snp.json', orient='records')

    return lineage_snp_df


def main():

    dna_snp_df = load_dna_snps()
    aa_snp_df = load_aa_snps()
    lineage_df = load_taxon_lineages()

    # print('getting NT consensus')
    get_consensus_dna_snps(dna_snp_df, lineage_df)
    
    # print('getting AA consensus')
    get_consensus_aa_snps(aa_snp_df, lineage_df)


if __name__ == '__main__':
    main()