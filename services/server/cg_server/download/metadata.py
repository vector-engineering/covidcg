# coding: utf-8

"""Download metadata and SNVs from selected sequences

Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""

import pandas as pd

from flask import make_response
from psycopg2 import sql

from cg_server.config import config
from cg_server.constants import constants
from cg_server.query import build_sequence_location_where_filter


def download_metadata(conn, req):

    with conn.cursor() as cur:
        sequence_where_filter = build_sequence_location_where_filter(req)

        # Fields that the user wants
        selected_fields = req.get("selected_fields", [])
        snv_format = req.get("snv_format", constants["SNV_FORMAT"]["POS_REF_ALT"])

        sequence_cols = [
            "Accession ID",
            "collection_date",
            "submission_date",
        ]
        sequence_cols_expr = [
            sql.SQL("m.{}").format(sql.Identifier(col)) for col in sequence_cols
        ]
        joins = []

        for grouping in config["group_cols"].keys():
            if grouping not in selected_fields:
                continue

            sequence_cols.append(grouping)
            sequence_cols_expr.append(sql.SQL("m.{}").format(sql.Identifier(grouping)))

        for field in list(config["metadata_cols"].keys()) + list(
            constants["GEO_LEVELS"].values()
        ):
            if field not in selected_fields:
                continue

            sequence_cols.append(field)
            sequence_cols_expr.append(
                sql.SQL(
                    """
                    {metadata_table_name}."value" as {field}
                    """
                ).format(
                    metadata_table_name=sql.Identifier("metadata_" + field),
                    field=sql.Identifier(field),
                )
            )
            joins.append(
                sql.SQL(
                    """
                    INNER JOIN {metadata_table_name} {metadata_table_name}
                        ON m.{field} = {metadata_table_name}."id"
                    """
                ).format(
                    metadata_table_name=sql.Identifier("metadata_" + field),
                    field=sql.Identifier(field),
                )
            )

        snv_id_to_pos_maps = {}
        snv_id_to_name_maps = {}

        for snv_field in ["dna", "gene_aa", "protein_aa"]:
            if snv_field not in selected_fields:
                continue

            snv_table = "sequence_{}_snp".format(snv_field)

            sequence_cols.append(snv_field + "_snp")
            sequence_cols_expr.append(
                sql.SQL("""{snv_table_short}."mutations" AS {snp_col} """).format(
                    snv_table_short=sql.Identifier(snv_field),
                    snp_col=sql.Identifier(snv_field + "_mutations"),
                )
            )

            joins.append(
                sql.SQL(
                    """
                INNER JOIN (
                    SELECT "sequence_id", "mutations"
                    FROM {snv_table}
                ) {snv_table_short} ON {snv_table_short}."sequence_id" = m."sequence_id"
                """
                ).format(
                    snv_table_short=sql.Identifier(snv_field),
                    snv_table=sql.Identifier(snv_table),
                )
            )

            # Get ID -> name map
            snv_agg_field = "snp_str"
            if snv_format == constants["SNV_FORMAT"]["POS_REF_ALT"]:
                snv_agg_field = "snp_str"
            elif snv_format == constants["SNV_FORMAT"]["REF_POS_ALT"]:
                snv_agg_field = "snv_name"

            cur.execute(
                sql.SQL(
                    """
                SELECT "id", "pos", {snv_agg_field}
                FROM {snp_table}
                """
                ).format(
                    snv_agg_field=sql.Identifier(snv_agg_field),
                    snp_table=sql.Identifier(snv_field + "_snp"),
                )
            )
            snv_props = pd.DataFrame.from_records(
                cur.fetchall(), columns=["id", "pos", "name"]
            )
            snv_id_to_pos_maps[snv_field] = dict(
                zip(snv_props["id"].values, snv_props["pos"].values)
            )
            snv_id_to_name_maps[snv_field] = dict(
                zip(snv_props["id"].values, snv_props["name"].values)
            )

        query = sql.SQL(
            """
            SELECT
                {sequence_cols_expr}
            FROM metadata m
            {joins}
            WHERE {sequence_where_filter}
            """
        ).format(
            sequence_cols_expr=sql.SQL(",").join(sequence_cols_expr),
            joins=sql.SQL("\n").join(joins),
            sequence_where_filter=sequence_where_filter,
        )
        # print(query.as_string(conn))

        cur.execute(query)

        res_df = pd.DataFrame.from_records(cur.fetchall(), columns=sequence_cols,)

        # Replace mutation IDs with names
        for snv_field in ["dna", "gene_aa", "protein_aa"]:
            if snv_field not in selected_fields:
                continue

            snv_field_col = snv_field + "_snp"

            pos_map = snv_id_to_pos_maps[snv_field]
            name_map = snv_id_to_name_maps[snv_field]

            # Sort mutations by position
            res_df.loc[:, snv_field_col] = res_df[snv_field_col].apply(
                lambda x: sorted(x, key=lambda _x: pos_map[_x])
            )

            # Serialize list of mutation IDs
            res_df.loc[:, snv_field_col] = res_df[snv_field_col].apply(
                lambda x: ";".join([name_map[_x] for _x in x])
            )

    return make_response(res_df.to_csv(index=False), 200, {"Content-Type": "text/csv"})
