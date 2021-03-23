# coding: utf-8

"""Get metadata values from keys

Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""

import psycopg2
from psycopg2 import sql

from cg_server.config import config


def get_metadata_fields(conn, req):

    table_queries = []
    for field in req.keys():
        if field not in config["metadata_cols"].keys():
            continue

        vals = req[field]
        if not vals:
            continue

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

    table_queries = sql.SQL("UNION ALL").join(table_queries)
    query = sql.SQL(
        """
        SELECT "field", json_object_agg("id", "value") as "map"
        FROM (
            {table_queries}
        ) a
        GROUP BY "field"
        """
    ).format(table_queries=table_queries)
    # print(query.as_string(conn))
    # print(dict(zip(req.keys(), [tuple(vals) for vals in req.values()])))

    with conn.cursor() as cur:
        cur.execute(
            query,
            # Fill in placeholders
            dict(zip(req.keys(), [tuple(vals) for vals in req.values()])),
        )
        return dict(cur.fetchall())
