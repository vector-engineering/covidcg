# coding: utf-8

"""Download metadata and mutations from selected sequences

Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""

import pandas as pd

from flask import make_response
from psycopg2 import sql

from cg_server.config import config
from cg_server.constants import constants
from cg_server.query import build_sequence_location_where_filter, get_loc_level_ids


def download_metadata(conn, req):
    selected_reference = req.get("selected_reference", None)
    if not selected_reference:
        raise Exception("No reference specified")

    with conn.cursor() as cur:
        sequence_where_filter = build_sequence_location_where_filter(
            # req.get("group_key", None),
            None,  # No group key to prevent putting reference WHERE condition
            #        on the metadata table itself
            #        Add the reference conditions manually instead via the join condition
            get_loc_level_ids(req),
            req.get("start_date", None),
            req.get("end_date", None),
            req.get("subm_start_date", None),
            req.get("subm_end_date", None),
            req.get("selected_metadata_fields", None),
            req.get("selected_group_fields", None),
            req.get("sequence_length", None),
            req.get("percent_ambiguous", None),
            selected_reference,
        )

        # Fields that the user wants
        selected_fields = req.get("selected_fields", [])
        mutation_format = req.get(
            "mutation_format", constants["MUTATION_FORMAT"]["POS_REF_ALT"]
        )

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
                    LEFT JOIN {metadata_table_name} {metadata_table_name}
                        ON m.{field} = {metadata_table_name}."id"
                    """
                ).format(
                    metadata_table_name=sql.Identifier("metadata_" + field),
                    field=sql.Identifier(field),
                )
            )

        # Reference column
        sequence_cols.append("reference")
        sequence_cols_expr.append(
            sql.SQL("{} AS {}").format(
                sql.Literal(selected_reference), sql.Identifier("reference")
            )
        )

        mutation_id_to_pos_map = {}
        mutation_id_to_name_map = {}

        for mutation_field in ["dna", "gene_aa", "protein_aa"]:
            if mutation_field not in selected_fields:
                continue

            mutation_table = "sequence_{}_mutation".format(mutation_field)

            sequence_cols.append(mutation_field + "_mutation")
            sequence_cols_expr.append(
                sql.SQL(
                    """{mutation_table_short}."mutations" AS {mutation_col} """
                ).format(
                    mutation_table_short=sql.Identifier(mutation_field),
                    mutation_col=sql.Identifier(mutation_field + "_mutations"),
                )
            )

            joins.append(
                sql.SQL(
                    """
                INNER JOIN (
                    SELECT "sequence_id", "reference", "mutations"
                    FROM {mutation_table}
                ) {mutation_table_short} ON 
                    {mutation_table_short}."sequence_id" = m."sequence_id" AND
                    {mutation_table_short}."reference" = {reference_name}
                """
                ).format(
                    mutation_table_short=sql.Identifier(mutation_field),
                    mutation_table=sql.Identifier(mutation_table),
                    reference_name=sql.Literal(selected_reference),
                )
            )

            # Get ID -> name map
            mutation_agg_field = "mutation_str"
            if mutation_format == constants["MUTATION_FORMAT"]["POS_REF_ALT"]:
                mutation_agg_field = "mutation_str"
            elif mutation_format == constants["MUTATION_FORMAT"]["REF_POS_ALT"]:
                mutation_agg_field = "mutation_name"

            cur.execute(
                sql.SQL(
                    """
                SELECT "id", "pos", {mutation_agg_field}
                FROM {mutation_table}
                """
                ).format(
                    mutation_agg_field=sql.Identifier(mutation_agg_field),
                    mutation_table=sql.Identifier(mutation_field + "_mutation"),
                )
            )
            mutation_props = pd.DataFrame.from_records(
                cur.fetchall(), columns=["id", "pos", "name"]
            )
            mutation_id_to_pos_map[mutation_field] = dict(
                zip(mutation_props["id"].values, mutation_props["pos"].values)
            )
            mutation_id_to_name_map[mutation_field] = dict(
                zip(mutation_props["id"].values, mutation_props["name"].values)
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

        res_df = pd.DataFrame.from_records(
            cur.fetchall(),
            columns=sequence_cols,
        )

        # Replace mutation IDs with names
        for mutation_field in ["dna", "gene_aa", "protein_aa"]:
            if mutation_field not in selected_fields:
                continue

            mutation_field_col = mutation_field + "_mutation"

            pos_map = mutation_id_to_pos_map[mutation_field]
            name_map = mutation_id_to_name_map[mutation_field]

            # Sort mutations by position
            res_df.loc[:, mutation_field_col] = res_df[mutation_field_col].apply(
                lambda x: sorted(x, key=lambda _x: pos_map[_x])
            )

            # Serialize list of mutation IDs
            res_df.loc[:, mutation_field_col] = res_df[mutation_field_col].apply(
                lambda x: ";".join([name_map[_x] for _x in x])
            )

    return make_response(res_df.to_csv(index=False), 200, {"Content-Type": "text/csv"})
