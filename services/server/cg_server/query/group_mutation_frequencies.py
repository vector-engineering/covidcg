# coding: utf-8

"""Get mutation frequencies for a group (lineage/clade)

Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""

import pandas as pd
from psycopg2 import sql


def query_group_mutation_frequencies(conn, req):
    group = req["group"]
    mutation_type = req["mutation_type"]
    consensus_threshold = req["consensus_threshold"]

    with conn.cursor() as cur:
        table_name = "{group}_frequency_{mutation_type}_mutation".format(
            group=group, mutation_type=mutation_type
        )
        mutation_table_name = "{mutation_type}_mutation".format(mutation_type=mutation_type)
        mutation_cols = ["pos", "ref", "alt", "mutation_name"]
        if mutation_type == "gene_aa":
            mutation_cols = ["gene",] + mutation_cols
        elif mutation_type == "protein_aa":
            mutation_cols = ["protein",] + mutation_cols

        mutation_cols_expr = sql.SQL(",\n").join(
            [
                sql.SQL("mut.{col}").format(col=sql.Identifier(col))
                for col in mutation_cols
            ]
        )

        cur.execute(
            sql.SQL(
                """
            SELECT
                f."name",
                f."count",
                f."fraction",
                f."mutation_id",
                {mutation_cols_expr}
            FROM {table_name} f
            JOIN {mutation_table_name} mut ON f."mutation_id" = mut."id"
            WHERE "fraction" >= %(consensus_threshold)s
            """
            ).format(
                mutation_cols_expr=mutation_cols_expr,
                table_name=sql.Identifier(table_name),
                mutation_table_name=sql.Identifier(mutation_table_name),
            ),
            {"consensus_threshold": consensus_threshold},
        )

        res = pd.DataFrame.from_records(
            cur.fetchall(),
            columns=["name", "count", "fraction", "mutation_id",] + mutation_cols,
        )

    # print(res)

    return res.to_json(orient="records")
