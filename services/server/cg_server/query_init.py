# coding: utf-8

"""Send initial data package to front-end
Consists of SNV definitions, metadata definitions, etc

Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""

import pandas as pd

from cg_server.config import config


def query_metadata_map(cur):

    table_queries = []
    for field in config["metadata_cols"].keys():
        table_queries.append(
            """
            SELECT '{field}' as "field", "id"::text, "value"
	        FROM "metadata_{field}"
            """.format(
                field=field
            )
        )
    for snp_type in ["dna", "gene_aa", "protein_aa"]:
        table_queries.append(
            """
            SELECT '{snp_type}_snp' as "field", "snp_str" as "id", "id"::text as "value"
	        FROM "{snp_type}_snp"
            """.format(
                snp_type=snp_type
            )
        )
    table_queries = " UNION ALL ".join(table_queries)

    cur.execute(
        """
        SELECT "field", json_object_agg("id", "value") as "map"
        FROM (
            {table_queries}
        ) a
        GROUP BY "field"
        """.format(
            table_queries=table_queries
        )
    )

    metadata_map = dict(cur.fetchall())
    # print(list(metadata_map.keys()))

    return metadata_map


def query_stats(cur):
    cur.execute(
        """
        SELECT "value"
        FROM "stats"        
        """
    )
    return cur.fetchone()[0]


def query_geo_select_tree(cur):
    cur.execute(
        """
        SELECT "value"
        FROM "geo_select_tree"        
        """
    )
    return cur.fetchone()[0]


def query_country_score(cur):
    cur.execute(
        """
        SELECT "value"
        FROM "country_score"
        """
    )
    return cur.fetchone()[0]


def query_init(conn):

    with conn.cursor() as cur:
        metadata_map = query_metadata_map(cur)
        stats = query_stats(cur)
        geo_select_tree = query_geo_select_tree(cur)
        country_score = query_country_score(cur)

    return {
        "metadata_map": metadata_map,
        **stats,
        "geo_select_tree": geo_select_tree,
        "country_score": country_score,
    }

