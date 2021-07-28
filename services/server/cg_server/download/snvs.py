# coding: utf-8

"""Download selected SNVs in long form

Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""

import pandas as pd
import psycopg2

from flask import make_response
from psycopg2 import sql

from cg_server.constants import constants
from cg_server.query.selection import (
    build_sequence_query,
    create_location_map_table,
    build_coordinate_filters,
)


def download_snvs(conn, req):

    with conn.cursor() as cur:
        location_ids = req.get("location_ids", None)
        start_date = pd.to_datetime(req.get("start_date", None))
        end_date = pd.to_datetime(req.get("end_date", None))

        subm_start_date = req.get("subm_start_date", "")
        subm_end_date = req.get("subm_end_date", "")
        subm_start_date = (
            None if subm_start_date == "" else pd.to_datetime(subm_start_date)
        )
        subm_end_date = None if subm_end_date == "" else pd.to_datetime(subm_end_date)

        selected_metadata_fields = req.get("selected_metadata_fields", None)
        dna_or_aa = req.get("dna_or_aa", None)
        coordinate_mode = req.get("coordinate_mode", None)
        coordinate_ranges = req.get("coordinate_ranges", None)
        selected_gene = req.get("selected_gene", None)
        selected_protein = req.get("selected_protein", None)

        sequence_query = build_sequence_query(
            location_ids=location_ids,
            start_date=start_date,
            end_date=end_date,
            subm_start_date=subm_start_date,
            subm_end_date=subm_end_date,
            selected_metadata_fields=selected_metadata_fields,
        )

        location_map_table_name = create_location_map_table(cur, location_ids)

        (snv_filter, snv_table) = build_coordinate_filters(
            conn,
            dna_or_aa,
            coordinate_mode,
            coordinate_ranges,
            selected_gene,
            selected_protein,
        )
        sequence_snv_table = "sequence_" + snv_table

        snv_cols = ["snv_name", "pos", "ref", "alt"]
        if dna_or_aa == constants["DNA_OR_AA"]["AA"]:
            snv_cols.append("nt_pos")
            if coordinate_mode == constants["COORDINATE_MODE"]["COORD_GENE"]:
                snv_cols.append("gene")
            else:
                snv_cols.append("protein")
        # Convert to SQL expressions
        snv_cols_expr = sql.SQL(",").join([sql.Identifier(col) for col in snv_cols])

        # welcome to CTE hell
        main_query = sql.SQL(
            """
            WITH "seq" AS (
                {sequence_query}
            ),
            "snp_data" AS (
                SELECT "id", {snv_cols_expr}
                FROM {snv_table}
                {snv_filter}
            )
            SELECT
                "seq"."Accession ID" AS "Accession ID",
                loc_map."name" AS "location",
                "seq"."collection_date",
                COALESCE(snp."snp_id", -1) AS "snp_id",
                {snv_cols_expr}
            FROM {sequence_snv_table} snp
            INNER JOIN "snp_data" ON snp_data.id = snp.snp_id
            RIGHT OUTER JOIN "seq" ON snp.sequence_id = "seq".id
            LEFT OUTER JOIN (
                SELECT "id", "name"
                FROM {location_map_table_name}
            ) loc_map ON "seq"."location_id" = loc_map."id"
            """
        ).format(
            sequence_query=sequence_query,
            location_map_table_name=sql.Identifier(location_map_table_name),
            sequence_snv_table=sql.Identifier(sequence_snv_table),
            snv_cols_expr=snv_cols_expr,
            snv_table=sql.Identifier(snv_table),
            snv_filter=snv_filter,
        )

        cur.execute(main_query)

        res_snv = pd.DataFrame.from_records(
            cur.fetchall(),
            columns=["Accession ID", "location", "collection_date", "snv_id"]
            + snv_cols,
        )

    res = make_response(res_snv.to_csv(index=False), 200, {"Content-Type": "text/csv"},)
    return res
