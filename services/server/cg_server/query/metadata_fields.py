# coding: utf-8

"""Get metadata values from keys

Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""

import psycopg2
from psycopg2 import sql

from cg_server.config import config


def query_metadata_fields(conn, conn_pool, req):

    table_queries = []
    placeholder_map = {}

    for field in req.keys():
        if field not in config["metadata_cols"].keys():
            continue

        vals = req[field]
        if not vals or type(vals) is not list:
            continue

        placeholder_map[field] = tuple(vals)

        table_queries.append(
            sql.SQL(
                """
            SELECT
                {field} as "field",
                "id"::text,
                "value"
	        FROM {field_table}
            WHERE "id" IN {vals}
            """
            ).format(
                # This should be safe since I check above that the field
                # is defined inside the config file
                field=sql.Literal(field),
                field_table=sql.Identifier("metadata_" + field),
                vals=sql.Placeholder(field),
            )
        )

    table_queries = sql.SQL(" UNION ALL ").join(table_queries)

    if not table_queries.as_string(conn):
        conn_pool.putconn(conn)
        return {}

    query = sql.SQL(
        """
        SELECT "field", json_object_agg("id", "value") as "map"
        FROM (
            {table_queries}
        ) a
        GROUP BY "field"
        """
    ).format(table_queries=table_queries)

    res = {}
    with conn.cursor() as cur:
        cur.execute(
            query,
            # Fill in placeholders
            placeholder_map,
        )
        res = dict(cur.fetchall())

    conn_pool.putconn(conn)
    return res
