#!/usr/bin/env python3
# coding: utf-8

'''Merge location and lineage data
'''

import json
import pandas as pd
import numpy as np

from pathlib import Path

from fasta import read_fasta_file
from process_geo import load_geo_data
from process_lineages import load_taxon_lineages

project_root_path = Path(__file__).resolve().parent.parent
data_dir = (project_root_path / 'data').resolve() # Resolve any symlinks --> absolute path

def main():
    location_df, unique_location_df = load_geo_data()

    # Drop duplicate GISAID entries
    location_df.drop_duplicates('gisaid_id', inplace=True)

    # Load in lineage data
    # lineage_snp_df = pd.read_csv(data_dir / 'lineage_snp.csv')
    taxon_lineage_df = load_taxon_lineages()

    # Join the location_id onto taxon_lineage_df using the gisaid_id as a key
    taxon_lineage_df['location_id'] = taxon_lineage_df['gisaid_id'].map(pd.Series(
        location_df['location_id'].values, index=location_df['gisaid_id'].values
    ))

    # Drop null location IDs (these were probably taxons that were filtered out
    # during preprocessing)
    taxon_lineage_df = taxon_lineage_df.loc[~pd.isnull(taxon_lineage_df['location_id']), :].reset_index(drop=True)

    # Convert location IDs to integers
    taxon_lineage_df['location_id'] = taxon_lineage_df['location_id'].astype(int)    

    # Group by location_id, lineage, and sample_date
    case_data_df = (
        taxon_lineage_df
        .groupby(['location_id', 'lineage', 'sample_date'], as_index=False)
        ['gisaid_id']
        .count()
        .rename(columns={
            'location_id': 'loc_id',
            'gisaid_id': 'cases',
            'sample_date': 'date'
        })
    )
    # print(case_data_df)


    # Save to disk
    print('Saving case data')
    case_data_df.to_csv(data_dir / 'case_data.csv', index=False)
    case_data_df.to_json(data_dir / 'case_data.json', orient='records')

    # Write the reference fasta file to json
    # Load the reference sequence
    ref_fasta_path = (data_dir / 'reference.fasta')
    with ref_fasta_path.open('r') as fp:
        lines = fp.readlines()
        ref = read_fasta_file(lines)
        _ref_seq = list(ref.values())[0]

    ref_json_path = data_dir / 'reference.json'
    with ref_json_path.open('w') as fp:
        fp.write(json.dumps({'ref_seq': _ref_seq}))


if __name__ == '__main__':
    main()