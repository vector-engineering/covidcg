import pandas as pd
import numpy as np
import json
import time
import datetime

from cassandra import query
from cassandra.cluster import Cluster

from collections import defaultdict


cluster = Cluster()
session = cluster.connect()

session.execute("USE cg_gisaid;")

start = time.time()
rows = session.execute("""
SELECT * FROM dna_mut WHERE location_id IN %s AND collection_date >= %s
""", (query.ValueSequence(list(range(800))), datetime.datetime.fromisoformat('2020-01-01')))
rows = list(rows)

print('fetched {} rows in {:.3f} ms'.format(len(rows), (time.time()-start)*1000))

start = time.time()

res = defaultdict(lambda: 0)
for i, row in enumerate(rows):
    # if i > 100:
    #     break
    # print(tuple(row.mutations))

    if row.location_id % 2:
        location = 'A'
    else:
        location = 'B'

    if row.mutations is not None:
        res[(location, row.collection_date, tuple(row.mutations))] += 1

print('aggregated {} rows in {:.3f} ms'.format(len(rows), (time.time()-start)*1000))

