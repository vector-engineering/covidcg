import pandas as pd
import numpy as np
import json
import time

from cassandra import ConsistencyLevel
from cassandra.cluster import Cluster
from cassandra.query import BatchStatement

cluster = Cluster()
session = cluster.connect()

# session.execute("""
# CREATE KEYSPACE IF NOT EXISTS cg_gisaid 
# WITH REPLICATION = {
#     'class': 'SimpleStrategy',
#     'replication_factor': 1
# };
# """)

session.execute("USE cg_gisaid;")

session.execute("DROP TABLE IF EXISTS dna_mut;")

# Primary key
# First field is the partition key, rest are clustering keys
# Constraints:
# - Primary key has to be unique
# - Partition key *needs* to be provided
# - Partition key can be a compound, i.e, (key1, key2), but then both have to 
#   be provided when doing WHERE filtering
# - Partition key can only be filtered on with = or IN operators (no inequalities)
# - Clustering keys can be ommitted, but cannot do a WHERE on a clustering key
#   that is downstream of another clustering key. i.e., if (key1, key2, key3), cannot
#   do WHERE key1 = a AND key3 = c. Must do: key1 = a AND key2 = b AND key3 = c
session.execute("""
CREATE TABLE IF NOT EXISTS dna_mut (
    sequence_id int,
    location_id int,
    collection_date timestamp,
    mutations set<text>,
    PRIMARY KEY (location_id, collection_date, sequence_id)
)
""")

df = pd.read_json('../example_data_genbank/case_data.json')
# print(df.head())

with open('../example_data_genbank/metadata_map.json') as fp:
    metadata_map = json.loads(fp.read())

dna_snp_map = {v: k for k, v in metadata_map['dna_snp'].items()}
gene_aa_snp_map = {v: k for k, v in metadata_map['gene_aa_snp'].items()}
protein_aa_snp_map = {v: k for k, v in metadata_map['protein_aa_snp'].items()}

collection_dates = pd.to_datetime(df['collection_date']).dt.to_pydatetime()
# print(df['collection_date'].dt.to_pydatetime())
df['dna_snp_str'] = df['dna_snp_str'].apply(lambda x: [dna_snp_map[_x] for _x in x])
df['gene_aa_snp_str'] = df['gene_aa_snp_str'].apply(lambda x: [gene_aa_snp_map[_x] for _x in x])
df['protein_aa_snp_str'] = df['protein_aa_snp_str'].apply(lambda x: [protein_aa_snp_map[_x] for _x in x])

print(df.head())

insert_batch_size = 10
n_batches = int(np.ceil(len(df) / insert_batch_size))

insert_mut = session.prepare("INSERT INTO dna_mut (collection_date, location_id, sequence_id, mutations) VALUES (?,?,?,?)")

start = time.time()
for b in range(n_batches):
    batch = BatchStatement(consistency_level=ConsistencyLevel.QUORUM)

    # print(b, time.time()-start)
    _df = df.iloc[b*insert_batch_size:(b+1)*insert_batch_size]
    
    for i, row in _df.iterrows():
        batch.add(insert_mut, (collection_dates[i], row['location_id'], i, set(row['dna_snp_str'])))
     
    session.execute(batch)

print('Inserted {} rows in {:.3f} s'.format(len(df), time.time()-start))
