# coding: utf-8

"""Get sequences/SNVs from database based on user selections

Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""

import pandas as pd
import psycopg2
import uuid

from psycopg2 import sql

from cg_server.config import config
from cg_server.constants import constants


def build_coordinate_filters(conn, dna_or_aa, coordinate_mode, coordinate_ranges):
    snv_cols = ["snp_str", "snv_name", "color", "pos", "ref", "alt"]
    snv_filter = []
    if dna_or_aa == constants["DNA_OR_AA"]["DNA"]:
        snv_table = "dna_snp"
    elif dna_or_aa == constants["DNA_OR_AA"]["AA"]:
        snv_cols.append("nt_pos")
        if coordinate_mode == constants["COORDINATE_MODES"]["COORD_GENE"]:
            snv_table = "gene_aa_snp"
            snv_cols.append("gene")
            snv_filter.append(
                sql.SQL('snp_data."gene" = {gene}').format(
                    gene=sql.Placeholder("selected_gene")
                )
            )
        elif coordinate_mode == constants["COORDINATE_MODES"]["COORD_PROTEIN"]:
            snv_table = "protein_aa_snp"
            snv_cols.append("protein")
            snv_filter.append(
                sql.SQL('snp_data."protein" = {protein}').format(
                    protein=sql.Placeholder("selected_protein")
                )
            )

    snv_filter = sql.SQL(" AND ").join(snv_filter)

    pos_filter = []
    pos_filter_injections = {}
    for i in range(len(coordinate_ranges)):
        pos_column = "pos" if dna_or_aa == constants["DNA_OR_AA"]["DNA"] else "nt_pos"
        pos_filter.append(
            sql.SQL(
                """
                (snp_data.{pos_column} >= {start} AND snp_data.{pos_column} <= {end})
                """
            ).format(
                pos_column=sql.Identifier(pos_column),
                start=sql.Placeholder("range_{}_start".format(i)),
                end=sql.Placeholder("range_{}_end".format(i)),
            )
        )
        # Build injection map
        pos_filter_injections["range_{}_start".format(i)] = coordinate_ranges[i][0]
        pos_filter_injections["range_{}_end".format(i)] = coordinate_ranges[i][1]

    pos_filter = sql.SQL(" OR ").join(pos_filter)

    # Compose final WHERE expression
    snv_filter = [snv_filter]
    # Only combine the snv_filter and pos_filter if the snv_filter exists
    if snv_filter[0].as_string(conn):
        snv_filter.append(sql.SQL(" AND "))

    snv_filter.append(pos_filter)
    snv_filter = sql.Composed(snv_filter)

    return snv_cols, snv_filter, snv_table, pos_filter_injections


def select_sequences(conn, cur, req):
    location_ids = req.get("location_ids", None)
    all_location_ids = sum(location_ids.values(), [])

    start_date = pd.to_datetime(req.get("start_date", None))
    end_date = pd.to_datetime(req.get("end_date", None))

    # Metadata filters will come in the form of a JSON object
    # of { metadata_field: [metadata_values] }
    selected_metadata_fields = req.get("selected_metadata_fields", None)
    metadata_filters = []
    # Build dictionary of metadata value tuples to inject
    metadata_vals = {}
    for md_key, md_vals in selected_metadata_fields.items():
        metadata_filters.append(
            sql.SQL("{field} IN {vals}").format(
                field=sql.Identifier(md_key), vals=sql.Placeholder(md_key)
            )
        )
        metadata_vals[md_key] = tuple(",".join([str(val) for val in md_vals]))

    metadata_filters = sql.SQL(" AND ").join(metadata_filters)
    if metadata_filters.as_string(conn):
        metadata_filters = sql.Composed([metadata_filters, sql.SQL(" AND ")])

    temp_table_name = "query_" + uuid.uuid4().hex

    cur.execute(
        sql.SQL(
            """
            CREATE TEMP TABLE {temp_table_name}
            ON COMMIT DROP
            AS (
                SELECT m.*
                FROM "metadata" m
                WHERE
                    {metadata_filters}
                    "collection_date" >= %(start_date)s AND
                    "collection_date" <= %(end_date)s AND
                    "location_id" IN %(location_ids)s
            );
            """
        ).format(
            temp_table_name=sql.Identifier(temp_table_name),
            metadata_filters=metadata_filters,
        ),
        {
            "start_date": start_date,
            "end_date": end_date,
            "location_ids": tuple(all_location_ids),
            # Inject metadata value filters
            **metadata_vals,
        },
    )

    return temp_table_name


def query_sequences(conn, req):

    with conn.cursor() as cur:

        temp_table_name = select_sequences(conn, cur, req)

        group_key = req.get("group_key", None)
        dna_or_aa = req.get("dna_or_aa", None)
        coordinate_mode = req.get("coordinate_mode", None)
        coordinate_ranges = req.get("coordinate_ranges", None)
        selected_gene = req.get("selected_gene", None)
        selected_protein = req.get("selected_protein", None)

        (
            snv_cols,
            snv_filter,
            snv_table,
            pos_filter_injections,
        ) = build_coordinate_filters(
            conn, dna_or_aa, coordinate_mode, coordinate_ranges,
        )
        sequence_snv_table = "sequence_" + snv_table

        snv_cols_selection = sql.SQL(",\n").join(
            [sql.SQL("snp.{}").format(sql.Identifier(col)) for col in snv_cols]
        )

        cur.execute(
            sql.SQL(
                """
                SELECT
                    seq."id" as "sequence_id", 
                    COALESCE(snp."snp_id", -1),
                    seq."Accession ID",
                    seq."location_id",
                    seq."collection_date",
                    {snv_cols_selection}
                FROM {temp_table_name} seq
                FULL JOIN (
                    SELECT
                        seq."id" as "sequence_id",
                        snp."snp_id",
                        snp_data.*
                    FROM {temp_table_name} seq
                    INNER JOIN {sequence_snv_table} snp ON seq."id" = snp."sequence_id"
                    INNER JOIN {snv_table} snp_data ON snp."snp_id" = snp_data."id"
                    WHERE {snv_filter}
                ) snp ON seq."id" = snp."sequence_id";
                """
            ).format(
                snv_cols_selection=snv_cols_selection,
                temp_table_name=sql.Identifier(temp_table_name),
                sequence_snv_table=sql.Identifier(sequence_snv_table),
                snv_table=sql.Identifier(snv_table),
                snv_filter=snv_filter,
            ),
            # Variable injections
            {
                "selected_gene": selected_gene,
                "selected_protein": selected_protein,
                **pos_filter_injections,
            },
        )

        res_snv = pd.DataFrame.from_records(
            cur.fetchall(),
            columns=[
                "sequence_id",
                "snp_id",
                "Accession ID",
                "location_id",
                "collection_date",
            ]
            + snv_cols,
        )

        sequence_cols = [
            "id",
            "Accession ID",
            "collection_date",
            "submission_date",
            "location_id",
        ]
        for field in config["metadata_cols"].keys():
            sequence_cols.append(field)
        for grouping in config["group_cols"].keys():
            sequence_cols.append(grouping)

        # Convert all to identifiers
        sequence_cols_expr = sql.SQL(",").join(
            [sql.SQL("q.{}").format(sql.Identifier(col)) for col in sequence_cols]
        )

        join_expr = sql.SQL("")
        if group_key != constants["GROUP_SNV"]:
            sequence_cols.append("color")
            sequence_cols_expr = sql.Composed(
                [sequence_cols_expr, sql.SQL(","), sql.SQL('group_data."color"')]
            )
            join_expr = sql.SQL(
                'JOIN {group_key} group_data ON q.{group_key} = group_data."name"'
            ).format(group_key=sql.Identifier(group_key))

        cur.execute(
            sql.SQL(
                """
                SELECT {sequence_cols} 
                FROM {temp_table_name} q
                {join_expr};
                """
            ).format(
                sequence_cols=sequence_cols_expr,
                temp_table_name=sql.Identifier(temp_table_name),
                join_expr=join_expr,
            )
        )
        res_df = pd.DataFrame.from_records(
            cur.fetchall(), index="id", columns=sequence_cols
        )

        # Clean up
        cur.execute(
            sql.SQL("DROP TABLE IF EXISTS {temp_table_name};").format(
                temp_table_name=sql.Identifier(temp_table_name)
            )
        )

    return res_df, res_snv

