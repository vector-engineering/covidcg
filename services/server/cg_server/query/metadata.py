# coding: utf-8

"""Get metadata values from keys

Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""

import pandas as pd
import psycopg2
import uuid

from cg_server.config import config
from cg_server.query.selection import build_sequence_query
from psycopg2 import sql


def query_metadata(conn, req):

    location_ids = req.get("location_ids", None)
    start_date = pd.to_datetime(req.get("start_date", None))
    end_date = pd.to_datetime(req.get("end_date", None))
    selected_metadata_fields = req.get("selected_metadata_fields", None)

    with conn.cursor() as cur:
        # First store the sequence query into a temp table for iterative
        # grouping on
        sequence_query = build_sequence_query(
            cur, conn, location_ids, start_date, end_date, selected_metadata_fields
        )
        temp_table_name = "sequence_selection_" + uuid.uuid4().hex
        cur.execute(
            sql.SQL(
                """
            CREATE TEMP TABLE {temp_table_name}
            ON COMMIT DROP
            AS ({sequence_query})
            """
            ).format(
                temp_table_name=sql.Identifier(temp_table_name),
                sequence_query=sequence_query,
            )
        )

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
