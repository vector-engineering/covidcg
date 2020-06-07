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
from process_snps import load_dna_snps, load_aa_snps, generate_snp_signatures
from reference import ref_seq, genes, gene_aa

project_root_path = Path(__file__).resolve().parent.parent
data_dir = (project_root_path / 'data').resolve() # Resolve any symlinks --> absolute path

def main():
    
    # ------------------------
    # Start with location data
    # ------------------------

    location_df, unique_location_df = load_geo_data()

    # Drop duplicate GISAID entries
    location_df.drop_duplicates('gisaid_id', inplace=True)

    # --------------------
    # Load in lineage data
    # --------------------

    # lineage_snp_df = pd.read_csv(data_dir / 'lineage_snp.csv')
    taxon_df = load_taxon_lineages()

    # Join the location_id onto taxon_df using the gisaid_id as a key
    taxon_df['location_id'] = taxon_df['gisaid_id'].map(pd.Series(
        location_df['location_id'].values, index=location_df['gisaid_id'].values
    ))

    # Drop null location IDs (these were probably taxons that were filtered out
    # during preprocessing)
    taxon_df = taxon_df.loc[~pd.isnull(taxon_df['location_id']), :].reset_index(drop=True)

    # Convert location IDs to integers
    taxon_df['location_id'] = taxon_df['location_id'].astype(int)    

    # ----------------
    # Load in SNP data
    # ----------------

    dna_snp_df = load_dna_snps()
    aa_snp_df = load_aa_snps()

    # Group SNPs by taxon, find and assign SNP signatures
    dna_snp_group_df, aa_snp_group_df = generate_snp_signatures(dna_snp_df, aa_snp_df)

    #print(dna_snp_group_df)
    #print(aa_snp_group_df)

    # Extract the GISAID id from the taxon column
    dna_snp_group_df['gisaid_id'] = dna_snp_group_df['taxon'].str.split('|').apply(lambda x: x[1])
    aa_snp_group_df['gisaid_id'] = aa_snp_group_df['taxon'].str.split('|').apply(lambda x: x[1])

    # Map SNPs to integer IDs
    print('Mapping SNPs to integers...', end='', flush=True)
    dna_snp_map = pd.Series(
        np.unique(np.concatenate(dna_snp_group_df['snp_str'].str.split(';').values).ravel())
    )
    aa_snp_map = pd.Series(
        np.unique(np.concatenate(aa_snp_group_df['snp_str'].str.split(';').values).ravel())
    )
    # Flip index and values
    dna_snp_map = pd.Series(dna_snp_map.index.values, index=dna_snp_map.values)
    aa_snp_map = pd.Series(aa_snp_map.index.values, index=aa_snp_map.values)

    # Save maps
    dna_snp_map.to_csv(data_dir / 'dna_snp_map.csv', index_label='snp', header=['id'])
    dna_snp_map.to_json(data_dir / 'dna_snp_map.json', orient='index')
    aa_snp_map.to_csv(data_dir / 'aa_snp_map.csv', index_label='snp', header=['id'])
    aa_snp_map.to_json(data_dir / 'aa_snp_map.json', orient='index')

    # Convert SNP strings to integer lists
    dna_snp_group_df['snp_str'] = dna_snp_group_df['snp_str'].str.split(';').apply(
        lambda x: ';'.join([str(dna_snp_map[a]) for a in x] if x else None
    ))
    dna_snp_group_df['snp_sig'] = dna_snp_group_df['snp_sig'].str.split(';').apply(
        lambda x: ';'.join([str(dna_snp_map[a]) for a in x] if x else None
    ))
    aa_snp_group_df['snp_str'] = aa_snp_group_df['snp_str'].str.split(';').apply(
        lambda x: ';'.join([str(aa_snp_map[a]) for a in x] if x else None
    ))
    aa_snp_group_df['snp_sig'] = aa_snp_group_df['snp_sig'].str.split(';').apply(
        lambda x: ';'.join([str(aa_snp_map[a]) for a in x] if x else None
    ))
    print('done. Saved SNP -> integer maps')

    # Join to taxon_df
    print('Joining to case dataframe...', end='', flush=True)
    taxon_df['dna_snp_str'] = taxon_df['gisaid_id'].map(pd.Series(
        dna_snp_group_df['snp_str'].values, index=dna_snp_group_df['gisaid_id'].values
    ))
    taxon_df['dna_snp_sig'] = taxon_df['gisaid_id'].map(pd.Series(
        dna_snp_group_df['snp_sig'].values, index=dna_snp_group_df['gisaid_id'].values
    ))
    taxon_df['aa_snp_str'] = taxon_df['gisaid_id'].map(pd.Series(
        aa_snp_group_df['snp_str'].values, index=aa_snp_group_df['gisaid_id'].values
    ))
    taxon_df['aa_snp_sig'] = taxon_df['gisaid_id'].map(pd.Series(
        aa_snp_group_df['snp_sig'].values, index=aa_snp_group_df['gisaid_id'].values
    ))
    print('done')

    # Convert from string back to a list of ints
    print('Converting SNP strings back into lists of integers...', end='', flush=True)
    taxon_df['dna_snp_str'] = taxon_df['dna_snp_str'].apply(
        lambda x: [int(snp) for snp in x.split(';')] if type(x) is str else None
    )
    taxon_df['aa_snp_str'] = taxon_df['aa_snp_str'].apply(
        lambda x: [int(snp) for snp in x.split(';')] if type(x) is str else None
    )
    taxon_df['dna_snp_sig'] = taxon_df['dna_snp_sig'].apply(
        lambda x: [int(snp) for snp in x.split(';')] if type(x) is str else None
    )
    taxon_df['aa_snp_sig'] = taxon_df['aa_snp_sig'].apply(
        lambda x: [int(snp) for snp in x.split(';')] if type(x) is str else None
    )
    print('done')
    
    print(taxon_df)

    # ---------------
    # Post-processing
    # ---------------
    print('Post-processing...', end='', flush=True)

    # Drop extra columns
    taxon_df.drop(columns=['name', 'taxon', 'SH-alrt', 'UFbootstrap', 'lineages_version', 'status', 'note'], inplace=True)

    print('done')

    # ------------
    # Save to disk
    # ------------

    print('Saving case data...', end='', flush=True)
    taxon_df.to_csv(data_dir / 'case_data2.csv', index=False)
    taxon_df.to_json(data_dir / 'case_data2.json', orient='records')
    

    # Write the reference fasta file to json

    ref_json_path = data_dir / 'reference.json'

    ref_obj = {
        'ref_seq': ref_seq,
        'gene_aa': gene_aa
    }

    with ref_json_path.open('w') as fp:
        fp.write(json.dumps(ref_obj))

    print('done')


if __name__ == '__main__':
    main()