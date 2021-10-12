# pip install google-cloud-bigtable

import argparse
import json
import os
import pandas as pd
from google.cloud import bigtable
from sklearn.preprocessing import MultiLabelBinarizer

# TODO: Move project_id, instance_id, and table_id into environment variables

# 1. Load metadata definitions
#    - Load example_data_genbank/metadata_map.json (JSON file)
#    - Load example_data_genbank/location_map.json (loaded via. pandas)
# 2. Load metadata dataframe (example_data_genbank/case_data.json)
# 3. Join metadata definitions onto the dataframe
# 4. Plan out/make column families
# 5. Push to bigtable

def main():

    project_id="test-project-1-2-3"
    instance_id="covid-cg"
    table_id="gisaid"

    # Create a Cloud Bigtable client.
    client = bigtable.Client(project=project_id)

    # Connect to an existing Cloud Bigtable instance.
    instance = client.instance(instance_id)

    # Open an existing table.
    table = instance.table(table_id)

    case_data = pd.read_json('example_data_genbank/case_data.json')
    location_map = pd.read_json('example_data_genbank/location_map.json')
    with open('example_data_genbank/metadata_map.json', "r") as fp:
        metadata_map = json.loads(fp.read())

    df = case_data.join(location_map, on="location_id")

    print(df.columns)
    print(metadata_map.keys())
    """
    'Accession ID', 

    # DATES
    'collection_date',
    'submission_date'

    # METADATA (integers, map is in metadata_map)
    'database', 'strain', 'host', 'isolation_source', 'biosample_accession',  
    'authors', 'publications', 
    
    # METADATA - string format
    'lineage', 

    # MUTATION IDs
    'dna_snp_str', 'gene_aa_snp_str',
    'protein_aa_snp_str', 
    
    # LOCATION INFORMATION
    'location_id', 'region', 'country', 'division',  
    'location'
    """

    # Create unique row key
    # [location_id]:[date]:[Accession ID]

    """
    Collection date is in ISO format (YYYY-MM-DD)
    """

    START_DATE = pd.to_datetime('2019-01-01')
    days_from_start = df['collection_date'].apply(lambda x: (pd.to_datetime(x) - START_DATE).days)

    bigtable_key = df['location_id'].astype(str) + ':' + days_from_start.astype(str) + ':' + df['Accession ID']

    # Helper function to map over lists and map to a dict.
    # Assumes the dictionary keys are ints.
    def geneHelper(l, d):
        for i in range(0, len(l)):
            l[i] = d[l[i]]
        return l
    # Map from metadata_map to df.
    print(df['dna_snp_str'].head())
    # Iterate over columns.
    for c in df.columns:
        # Map columns to the metadata map where the key is the int and the value
        # is the value e.g. Homo sapiens, etc.
        try:
            d = metadata_map[c.replace('_str', '')]
            # Reverse item and key order so that ints can be used as keys.
            #mappingData = {v: k for k, v in d.items()}
            d = {int(k): v for k, v in d.items()}
            df[c] = df[c].map(d)
        # Map columns where the columns contain lists. This is for all of the
        # sequences where 0, 1, or more can be present.
        except:
            try:
                d = metadata_map[c.replace('_str', '')]
                d = {v: k for k, v in d.items()}
                df[c] = df[c].apply(lambda x: geneHelper(x, d))
            except:
                pass
    print(df['dna_snp_str'].head())
    # TODO: config_genbank.yaml as map for ints in df

    # Create core metadata family.
    # Create sparse matrix for from location and location id.
    # TODO: Map for locations.
    coreMetadata = pd.get_dummies(df[['location', 'location_id']])
    coreMetadata[['collection_data', 'submission_date', 'Accession ID']] = df[[
        'collection_date', 'submission_date', 'Accession ID']]
    # Create other metadata family.
    otherMetadata = df[['database', 'strain', 'host',
                        'isolation_source', 'biosample_accession']]
    # Mutations family.
    mlb = MultiLabelBinarizer()
    # Binarize all mutations, and add _type suffix to column name to prevent collisions.
    dnaMuts = pd.DataFrame(mlb.fit_transform(
        df['dna_snp_str']),
        columns=mlb.classes_, index=df.index).add_suffix('_dna')
    geneMuts = pd.DataFrame(mlb.fit_transform(
        df['gene_aa_snp_str']),
        columns=mlb.classes_, index=df.index).add_suffix('_gene')
    proteinMuts = pd.DataFrame(mlb.fit_transform(
        df['protein_aa_snp_str']),
        columns=mlb.classes_, index=df.index).add_suffix('_protein')
    print(dnaMuts)
    mutationFam = pd.concat([dnaMuts, geneMuts, proteinMuts])
    print(mutationFam.shape)

    """
    `metadata_cols` field in config_genbank.yaml has a list of metadata columns that will need to be "unmapped" from integers back to strings. The int->string map itself is in `metadata_map.json` (already loaded in here)

    Map mutations back from integers -> strings. This map should also be in `metadata_map.json` file. Three levels of mutations, 'dna', 'gene_aa', and 'protein_aa'. These are currently integers in the main DF, in columns: 'dna_snp_str', 'gene_aa_snp_str',
    'protein_aa_snp_str'
    - use pandas dummy variables function
    
    
    Separate into column families:
        - core_metadata
            - location
            - collection date
            - submission date
            - Accession ID
            - location_id
        - other metadata (less frequently accessed)
            - 'database', 'strain', 'host', 'isolation_source', 'biosample_accession',  
    'authors', 'publications', 
        - DNA mutations
        - GENE AA mutations
        - PROTEIN AA mutations
        - Location
            - region (continent)
            - country
            - division (state, province)
            - location (county, city, etc)
        - sequence
            - whole genome sequence
            - NOT NOW
    
    """

    # row_key = 'r1'
    # row = table.read_row(row_key.encode('utf-8'))

    # column_family_id = 'metadata'
    # column_id = 'c1'.encode('utf-8')
    # value = row.cells[column_family_id][column_id][0].value.decode('utf-8')

    # print('Row key: {}\nData: {}'.format(row_key, value))


if __name__ == '__main__':
    # parser = argparse.ArgumentParser(
    #     description=__doc__,
    #     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # parser.add_argument('project_id', help='Your Cloud Platform project ID.')
    # parser.add_argument(
    #     'instance_id', help='ID of the Cloud Bigtable instance to connect to.')
    # parser.add_argument(
    #     '--table',
    #     help='Existing table used in the quickstart.',
    #     default='my-table')

    # args = parser.parse_args()
    # main(args.project_id, args.instance_id, args.table)
    main()
