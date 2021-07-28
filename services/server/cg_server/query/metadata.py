# coding: utf-8

"""Get metadata values from keys

Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""

import pandas as pd
import psycopg2

from cg_server.config import config
from cg_server.query.selection import create_sequence_temp_table
from psycopg2 import sql


def query_metadata(conn, req):

    with conn.cursor() as cur:
        temp_table_name = create_sequence_temp_table(cur, req)

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
                    counts.*,
                    def."value" as "val_str"
                FROM (
                    SELECT {metadata_field_literal} as "field", {metadata_field_ident} as "val_id", count({metadata_field_ident})
                    FROM {temp_table_name}
                    GROUP BY {metadata_field_ident}
                ) counts
                INNER JOIN {metadata_def_table} def ON counts."val_id" = def."id"
                """
                ).format(
                    metadata_field_literal=sql.Literal(field),
                    metadata_field_ident=sql.Identifier(field),
                    temp_table_name=sql.Identifier(temp_table_name),
                    metadata_def_table=sql.Identifier("metadata_" + field),
                )
            )

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
