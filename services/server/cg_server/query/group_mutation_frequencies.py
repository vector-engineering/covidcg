# coding: utf-8

"""Get mutation frequencies for a group (lineage/clade)

Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""

import pandas as pd
from psycopg2 import sql

from cg_server.constants import constants
from cg_server.query.selection import build_sequence_location_where_filter


def query_group_mutation_frequencies(conn, req):
    group = req["group"]
    mutation_type = req["mutation_type"]
    consensus_threshold = req["consensus_threshold"]

    with conn.cursor() as cur:
        table_name = "{group}_frequency_{mutation_type}_mutation".format(
            group=group, mutation_type=mutation_type
        )
        mutation_table_name = "{mutation_type}_mutation".format(
            mutation_type=mutation_type
        )
        mutation_cols = ["pos", "ref", "alt", "mutation_name", "mutation_str"]
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
            columns=["name", "count", "fraction", "mutation_id"] + mutation_cols,
        )

    # print(res)

    return res.to_json(orient="records")


def query_group_mutation_frequencies_dynamic(conn, req):
    group_key = req.get("group_key", None)
    mutation_type = req.get("mutation_type", "dna")
    consensus_threshold = req.get("consensus_threshold", 0.0)

    mutation_cols = ["pos", "ref", "alt", "mutation_name", "mutation_str"]
    if mutation_type == "dna":
        mutation_table = "dna_mutation"
    elif mutation_type == "gene_aa":
        mutation_table = "gene_aa_mutation"
        mutation_cols = ["gene",] + mutation_cols
    elif mutation_type == "protein_aa":
        mutation_table = "protein_aa_mutation"
        mutation_cols = ["protein",] + mutation_cols

    sequence_where_filter = build_sequence_location_where_filter(req)
    sequence_mutation_table = "sequence_" + mutation_table

    mutation_cols_expr = sql.SQL(",\n").join(
        [sql.SQL("mut.{col}").format(col=sql.Identifier(col)) for col in mutation_cols]
    )

    with conn.cursor() as cur:
        cur.execute(
            sql.SQL(
                """
            WITH "group_counts" AS (
                SELECT {group_col}, COUNT("sequence_id")
                FROM {sequence_mutation_table}
                WHERE {sequence_where_filter}
                GROUP BY {group_col}
            ),
            "group_muts" AS (
                SELECT
                    {group_col}, "mutation", COUNT("sequence_id")
                FROM (
                    SELECT "sequence_id", {group_col}, UNNEST("mutations") as "mutation"
                    FROM {sequence_mutation_table}
                    WHERE {sequence_where_filter}
                ) "group_muts"
                GROUP BY {group_col}, "mutation"
            )
            SELECT
                "group_muts".{group_col} AS "name", 
                "group_muts"."count" AS "count", 
                ("group_muts"."count"::REAL / "group_counts"."count"::REAL) AS "fraction",
                "mutation" AS "mutation_id",
                {mutation_cols_expr}
            FROM "group_muts"
            JOIN "group_counts" ON "group_counts".{group_col} = "group_muts".{group_col}
            JOIN {mutation_table} "mut" ON "group_muts"."mutation" = "mut"."id"
            WHERE ("group_muts"."count"::REAL / "group_counts"."count"::REAL) >= %(consensus_threshold)s
            ORDER BY "name" ASC, "count" DESC
        """
            ).format(
                group_col=sql.Identifier(group_key),
                sequence_mutation_table=sql.Identifier(sequence_mutation_table),
                mutation_table=sql.Identifier(mutation_table),
                sequence_where_filter=sequence_where_filter,
                mutation_cols_expr=mutation_cols_expr,
            ),
            {"consensus_threshold": consensus_threshold},
        )

        res = pd.DataFrame.from_records(
            cur.fetchall(),
            columns=["name", "count", "fraction", "mutation_id"] + mutation_cols,
        )

    return res.to_json(orient="records")

