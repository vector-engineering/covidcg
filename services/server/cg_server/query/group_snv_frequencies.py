# coding: utf-8

"""Get SNV frequencies for a group (lineage/clade)

Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""

import pandas as pd
import psycopg2

from psycopg2 import sql


def query_group_snv_frequencies(conn, req):
    group = req["group"]
    snv_type = req["snv_type"]
    consensus_threshold = req["consensus_threshold"]

    with conn.cursor() as cur:
        table_name = "{group}_frequency_{snv_type}_snp".format(
            group=group, snv_type=snv_type
        )
        snv_table_name = "{snv_type}_snp".format(snv_type=snv_type)
        snv_cols = ["pos", "ref", "alt", "snv_name"]
        if snv_type == "gene_aa":
            snv_cols = ["gene",] + snv_cols
        elif snv_type == "protein_aa":
            snv_cols = ["protein",] + snv_cols

        snv_cols_expr = sql.SQL(",\n").join(
            [sql.SQL("snv.{col}").format(col=sql.Identifier(col)) for col in snv_cols]
        )

        cur.execute(
            sql.SQL(
                """
            SELECT
                f."name",
                f."count",
                f."fraction",
                {snv_cols_expr}
            FROM {table_name} f
            JOIN {snv_table_name} snv ON f."snv_id" = snv."id"
            WHERE "fraction" >= %(consensus_threshold)s
            """
            ).format(
                snv_cols_expr=snv_cols_expr,
                table_name=sql.Identifier(table_name),
                snv_table_name=sql.Identifier(snv_table_name),
            ),
            {"consensus_threshold": consensus_threshold},
        )

        res = pd.DataFrame.from_records(
            cur.fetchall(), columns=["name", "count", "fraction",] + snv_cols
        )

    # print(res)

    return res.to_json(orient="records")
