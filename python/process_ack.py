#!/usr/bin/env python3
# coding: utf-8

'''Process acknowledgement files

Author: Albert Chen (Deverman Lab - Broad Institute)
'''

import numpy as np
import pandas as pd

from pathlib import Path

project_root_path = Path(__file__).resolve().parent.parent
data_dir = (project_root_path / 'data').resolve() # Resolve any symlinks --> absolute path

def process_ack():
    '''COLLECT ACKNOWLEDGEMENTS	
    '''

    ack_files = sorted((data_dir / 'acknowledgements').glob('*.xls'))
    print('Collecting {} acknowledgement files...'.format(len(ack_files)), end='', flush=True)	
    ack_df = pd.DataFrame()	
    for f in ack_files:	
        _df = pd.read_excel(f, skiprows=[0, 1, 3])	
        ack_df = pd.concat([ack_df, _df], ignore_index=True)	

    # Rename columns, and keep a subset
    ack_df = (
        ack_df
        .rename(columns={
            'Accession ID': 'gisaid_id', 
            'Originating lab': 'originating_lab', 
            'Submitting lab': 'submitting_lab',
            'Authors': 'authors'
        })
        [['gisaid_id', 'originating_lab', 'submitting_lab', 'authors']]
    )
    print('done')

    # Get unique labs, by originating lab, submitting lab, and authors
    # Where the dataframe index becomes the acknowledgement ID
    print('Getting unique acknowledgement entries...', end='', flush=True)
    unique_ack_df = (
        ack_df
        .drop(columns=['gisaid_id'])
        .drop_duplicates(['originating_lab', 'submitting_lab', 'authors'])
        .reset_index(drop=True)
        # Create index column
        .reset_index()
        # Flip index and text columns, for fast lookups
        .set_index(['originating_lab', 'submitting_lab', 'authors'])
    )
    print('done')
    # print(unique_ack_df)

    # Map acknowledgement IDs back to GISAID ids
    print('Mapping acknowledgement IDs to GISAID IDs...', end='', flush=True)
    ack_df['ack_id'] = -1
    for i, row in ack_df.iterrows():
        ack_df.loc[i, 'ack_id'] = (
            unique_ack_df.loc[
                row['originating_lab'], 
                row['submitting_lab'], 
                row['authors']
            ]['index']
        )

    # Turn into a series
    ack_df = pd.Series(ack_df['ack_id'].values, index=ack_df['gisaid_id'].values)

    print('done')

    # print(ack_df)

    print('Saving acknowledgement files...', end='', flush=True)
    
    ack_df.to_csv(data_dir / 'taxon_acknowledgements.csv', header=['ack_id'], index_label='gisaid_id')	
    ack_df.to_json(data_dir / 'taxon_acknowledgements.json')	

    # Drop MultiIndex to columns, make index the real index again
    unique_ack_df = unique_ack_df.reset_index().set_index('index')
    unique_ack_df.to_csv(data_dir / 'acknowledgement_map.csv', index=True, index_label='ack_id')
    unique_ack_df.to_json(data_dir / 'acknowledgement_map.json', orient='index')	

    print('done', flush=True)	

if __name__ == '__main__':
    process_ack()
