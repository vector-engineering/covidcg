import json
import numpy as np
import os
import pandas as pd
import time

from cassandra import ConsistencyLevel
from cassandra.cluster import Cluster
from cassandra.query import BatchStatement

from pathlib import Path
from yaml import load, dump
try:
    from yaml import CLoader as Loader, CDumper as Dumper
except ImportError:
    from yaml import Loader, Dumper

from cg_server.color import get_categorical_colormap
from cg_server.config import config

# root/services/server/cg_server/db_seed/seed.py
# root/cassandra_test/cg_server/seed.py
project_root = Path(__file__).parent.parent.parent
data_path = Path(os.getenv("DATA_PATH", project_root / config["data_folder"]))
static_data_path = Path(
    os.getenv("STATIC_DATA_PATH", project_root / config["static_data_folder"])
)

def seed(session):
    print('SEEDING DATABASE')

    # TODO:
    # Google cloud deployment, test performance
    # Genome table?

    session.execute("""
    CREATE KEYSPACE IF NOT EXISTS cg_gisaid 
    WITH REPLICATION = {
        'class': 'SimpleStrategy',
        'replication_factor': 1
    };
    """)

    session.execute("USE cg_gisaid;")

    df = pd.read_json(data_path / 'case_data.json')
    # print(df.head())

    with (data_path / 'metadata_map.json').open('r') as fp:
        metadata_map = json.loads(fp.read())
    
    location_map = pd.read_json(data_path / 'location_map.json')

    # WRITE BLOBS
    start = time.time()
    print('Writing blobs')

    session.execute("DROP TABLE IF EXISTS blobs;")
    session.execute("""
    CREATE TABLE blobs (
        bkey text PRIMARY KEY,
        bvalue blob
    );""")

    # Write metadata map
    session.execute("INSERT INTO blobs (bkey, bvalue) VALUES (%s, %s)", ('metadata_map', bytearray(json.dumps(metadata_map), 'utf-8')))

    # Write location map
    with (data_path / 'location_map.json').open('r') as fp:
        session.execute("INSERT INTO blobs (bkey, bvalue) VALUES (%s, %s)", ('location_map', bytearray(fp.read(), 'utf-8')))

    # Write geo select tree
    with (data_path / 'geo_select_tree.json').open('r') as fp:
        session.execute("INSERT INTO blobs (bkey, bvalue) VALUES (%s, %s)", ('geo_select_tree', bytearray(fp.read(), 'utf-8')))

    # Groups
    with (data_path / "group_consensus_snps.json").open("r") as fp:
        group_consensus_snps = json.loads(fp.read())

    groups = []
    for grouping in group_consensus_snps.keys():
        for group_name in list(group_consensus_snps[grouping].keys()):
            color = get_categorical_colormap(group_name)
            groups.append({'group': grouping, 'name': group_name, 'color': color})
    session.execute("INSERT INTO blobs (bkey, bvalue) VALUES (%s, %s)", ('groups', bytearray(json.dumps(groups), 'utf-8')))
    
    print('Finished writing blobs in {:.3f} s'.format(time.time()-start))

    dna_snp_map = {v: k for k, v in metadata_map['dna_snp'].items()}
    gene_aa_snp_map = {v: k for k, v in metadata_map['gene_aa_snp'].items()}
    protein_aa_snp_map = {v: k for k, v in metadata_map['protein_aa_snp'].items()}

    collection_dates = pd.to_datetime(df['collection_date']).dt.to_pydatetime()
    # print(df['collection_date'].dt.to_pydatetime())
    df['dna_snp_str'] = df['dna_snp_str'].apply(lambda x: [dna_snp_map[_x] for _x in x])
    df['gene_aa_snp_str'] = df['gene_aa_snp_str'].apply(lambda x: [gene_aa_snp_map[_x] for _x in x])
    df['protein_aa_snp_str'] = df['protein_aa_snp_str'].apply(lambda x: [protein_aa_snp_map[_x] for _x in x])
    df = df.rename(columns={
        'dna_snp_str': 'dna_snp',
        'gene_aa_snp_str': 'gene_aa_snp',
        'protein_aa_snp_str': 'protein_aa_snp'
    })
    # print(df.head())

    # GROUP MUTATION FREQUENCIES
    start = time.time()
    print('Writing group mutation frequencies')

    with (data_path / "group_snv_frequencies.json").open("r") as fp:
            group_snv_frequencies = json.loads(fp.read())

    session.execute("DROP TABLE IF EXISTS group_mutations;")
    session.execute("""
    CREATE TABLE group_mutations (
        mutation_type text,
        group text,
        mutation text,
        count int,
        fraction float,
        PRIMARY KEY ((mutation_type, group), mutation)
    );
    """)

    insert_mut_freq = session.prepare("""
    INSERT INTO group_mutations (mutation_type, group, mutation, count, fraction) VALUES (?,?,?,?,?);
    """)

    snp_fields = ["dna", "gene_aa", "protein_aa"]
    for grouping in group_snv_frequencies.keys():
        for snp_field in snp_fields:

            if snp_field == 'dna':
                snp_map = dna_snp_map
            elif snp_field == 'gene_aa':
                snp_map = gene_aa_snp_map
            elif snp_field == 'protein_aa':
                snp_map = protein_aa_snp_map
            
            for item in group_snv_frequencies[grouping][snp_field]:
                session.execute(insert_mut_freq, (snp_field, item['group'], snp_map[item['snv_id']], item['count'], item['fraction']))

    print('Finished writing group mutation frequencies in {:.3f} s'.format(time.time()-start))


    insert_batch_size = 10
    n_batches = int(np.ceil(len(df) / insert_batch_size))

    metadata_cols = list(config['metadata_cols'].keys())
    metadata_col_defs = ',\n'.join(['{} int'.format(col) for col in metadata_cols])

    mutation_types = ['dna_snp', 'gene_aa_snp', 'protein_aa_snp']

    for mutation_type in mutation_types:
        print('Inserting {}'.format(mutation_type))

        session.execute("DROP TABLE IF EXISTS {};".format(mutation_type))

        """
        Primary key
        First field is the partition key, rest are clustering keys
        Constraints:
        - Primary key has to be unique
        - Partition key *needs* to be provided
        - Partition key can be a compound, i.e, (key1, key2), but then both have to 
        be provided when doing WHERE filtering
        - Partition key can only be filtered on with = or IN operators (no inequalities)
        - Clustering keys can be ommitted, but cannot do a WHERE on a clustering key
        that is downstream of another clustering key. i.e., if (key1, key2, key3), cannot
        do WHERE key1 = a AND key3 = c. Must do: key1 = a AND key2 = b AND key3 = c
        """

        session.execute("""
        CREATE TABLE IF NOT EXISTS {mutation_type} (
            sequence_id int,
            location_id int,
            collection_date timestamp,
            {metadata_col_defs},
            mutations set<text>,
            PRIMARY KEY (location_id, collection_date, sequence_id)
        )
        """.format(
            mutation_type=mutation_type,
            metadata_col_defs=metadata_col_defs
        ))

        insert_mut = session.prepare("""
        INSERT INTO {mutation_type} (collection_date, location_id, sequence_id, mutations, {metadata_cols}) 
        VALUES (?,?,?,?,{qmarks})
        """.format(
            mutation_type=mutation_type,
            metadata_cols=','.join(metadata_cols),
            qmarks=','.join(['?']*len(metadata_cols))
        ))

        start = time.time()
        for b in range(n_batches):
            batch = BatchStatement(consistency_level=ConsistencyLevel.QUORUM)

            # print(b, time.time()-start)
            _df = df.iloc[b*insert_batch_size:(b+1)*insert_batch_size]
            
            for i, row in _df.iterrows():
                metadata_vals = [row[col] for col in metadata_cols]
                batch.add(insert_mut, (collection_dates[i], row['location_id'], i, set(row[mutation_type]), *metadata_vals))
            
            session.execute(batch)

        print('Inserted {} sequences in {:.3f} s'.format(len(df), time.time()-start))

        # CREATE INDICES
        for col in metadata_cols:
            index_name = mutation_type + '_' + col + '_idx'
            session.execute("DROP INDEX IF EXISTS {index_name}".format(index_name=index_name))
            session.execute("CREATE INDEX {index_name} ON {mutation_type}({metadata_col});".format(
                index_name=index_name,
                mutation_type=mutation_type,
                metadata_col=col
            ))

    # GROUP COLS
    group_cols = list(config['group_cols'].keys())
    for group_col in group_cols:
        print('Inserting {}'.format(group_col))

        session.execute("DROP TABLE IF EXISTS {group_col};".format(group_col=group_col))
        session.execute("""
        CREATE TABLE IF NOT EXISTS {group_col} (
            sequence_id int,
            location_id int,
            collection_date timestamp,
            {metadata_col_defs},
            group text,
            PRIMARY KEY (location_id, collection_date, sequence_id)
        )
        """.format(
            group_col=group_col,
            metadata_col_defs=metadata_col_defs
        ))

        insert_group = session.prepare("""
        INSERT INTO {group_col} (collection_date, location_id, sequence_id, group, {metadata_cols}) 
        VALUES (?,?,?,?,{qmarks})
        """.format(
            group_col=group_col,
            metadata_cols=','.join(metadata_cols),
            qmarks=','.join(['?']*len(metadata_cols))
        ))

        start = time.time()
        for b in range(n_batches):
            batch = BatchStatement(consistency_level=ConsistencyLevel.QUORUM)

            # print(b, time.time()-start)
            _df = df.iloc[b*insert_batch_size:(b+1)*insert_batch_size]
            
            for i, row in _df.iterrows():
                metadata_vals = [row[col] for col in metadata_cols]
                batch.add(insert_group, (collection_dates[i], row['location_id'], i, row[group_col], *metadata_vals))
            
            session.execute(batch)

        print('Inserted {} sequences in {:.3f} s'.format(len(df), time.time()-start))
