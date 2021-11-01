# coding: utf-8

"""Get metadata values from keys

Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""

import pandas as pd

from cg_server.config import config
from cg_server.query import build_sequence_location_where_filter
from psycopg2 import sql


def query_metadata(conn, req):
    with conn.cursor() as cur:

        sequence_where_filter = build_sequence_location_where_filter(req)

        # Iterate over each metadata column, and aggregate counts
        # per metadata value
        # Also, grab the metadata value string for each corresponding
        # metadata value ID
        # Union these all together so we only make one trip
        metadata_queries = []
        for field in config["metadata_cols"].keys():
            metadata_queries.append(
                sql.SQL(
                    """
                SELECT
                    counts."field",
                    counts."val_id",
                    counts."count",
                    def."value" AS "val_str"
                FROM (
                    SELECT 
                        {metadata_field_literal} AS "field", 
                        {metadata_field_ident} AS "val_id", 
                        count({metadata_field_ident}) AS "count"
                    FROM "metadata"
                    WHERE {sequence_where_filter}
                    GROUP BY {metadata_field_ident}
                ) counts
                INNER JOIN {metadata_def_table} def ON counts."val_id" = def."id"
                """
                ).format(
                    metadata_field_literal=sql.Literal(field),
                    metadata_field_ident=sql.Identifier(field),
                    metadata_def_table=sql.Identifier("metadata_" + field),
                    sequence_where_filter=sequence_where_filter,
                )
            )
            # print(metadata_queries[-1].as_string(conn))

        metadata_queries = sql.SQL(
            """
            UNION ALL
            """
        ).join(metadata_queries)

        cur.execute(metadata_queries)

        res = pd.DataFrame.from_records(
            cur.fetchall(), columns=["field", "val_id", "count", "val_str"]
        ).to_json(orient="records")

    return res
