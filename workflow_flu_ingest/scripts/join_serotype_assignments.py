#!/usr/bin/env python3
# coding: utf-8

"""Join serotype assignments to metadata

Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""

import argparse
import json
import pandas as pd

def main():

    parser = argparse.ArgumentParser()

    parser.add_argument('--metadata-virus', type=str, required=True, help='Metadata virus CSV file')
    parser.add_argument('--serotype-assignments', type=str, required=True, help='Serotype assignments CSV file')
    parser.add_argument('--out-metadata-virus', type=str, required=True, help='Output metadata virus CSV file')

    args = parser.parse_args()

    metadata_virus = pd.read_csv(args.metadata_virus, index_col='virus_name')
    serotype_assignments = pd.read_csv(args.serotype_assignments, index_col='Accession ID')

    # Read Accession ID JSONs
    metadata_virus.loc[:, 'accession_ids'] = metadata_virus['accession_ids'].apply(json.loads)

    # Join virus names
    serotype_assignments = serotype_assignments.merge(
        metadata_virus[['accession_ids']].explode('accession_ids').reset_index(),
        how='inner', left_index=True, right_on='accession_ids', copy=False, sort=False
    )

    # Join assignments onto virus DF
    metadata_virus = metadata_virus.join(
        (
            serotype_assignments
            .drop(columns=['accession_ids'])
            .rename(columns={'genus': 'assign_genus', 'serotype': 'assign_serotype'})
            .set_index('virus_name')
        ),
        how='left'
    )
    metadata_virus.loc[:, 'serotype'] = metadata_virus['serotype'].combine_first(metadata_virus['assign_serotype'])
    metadata_virus.drop(columns=['assign_genus', 'assign_serotype'], inplace=True)
    
    # Remove viruses without serotype assignment
    metadata_virus.drop(metadata_virus.index[metadata_virus['serotype'].isna()], inplace=True)

    # Save to disk
    # First reserialize accession IDs
    metadata_virus.loc[:, 'accession_ids'] = metadata_virus['accession_ids'].apply(json.dumps)
    metadata_virus.to_csv(args.out_metadata_virus)


if __name__ == '__main__':
    main()
