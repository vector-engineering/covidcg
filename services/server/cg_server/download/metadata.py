# coding: utf-8

"""Download metadata and SNVs from selected sequences

Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""

import pandas as pd
import psycopg2

from flask import make_response
from psycopg2 import sql

from cg_server.config import config
from cg_server.constants import constants
from cg_server.query.selection import select_sequences


def download_metadata(conn, req):

    with conn.cursor() as cur:
        temp_table_name = select_sequences(conn, cur, req)

        # Fields that the user wants
        selected_fields = req.get("selected_fields", [])
        snv_format = req.get("snv_format", constants["SNV_FORMAT"]["POS_REF_ALT"])

        sequence_cols = [
            "Accession ID",
            "collection_date",
            "submission_date",
        ]
        sequence_cols_expr = [
            sql.SQL("q.{}").format(sql.Identifier(col)) for col in sequence_cols
        ]
        metadata_joins = []

        # Location columns
        for col in list(constants["GEO_LEVELS"].values()):
            if col not in selected_fields:
                continue

            sequence_cols.append(col)
            sequence_cols_expr.append(sql.SQL("loc.{}").format(sql.Identifier(col)))

        for grouping in config["group_cols"].keys():
            if grouping not in selected_fields:
                continue

            sequence_cols.append(grouping)
            sequence_cols_expr.append(sql.SQL("q.{}").format(sql.Identifier(grouping)))

        for field in config["metadata_cols"].keys():
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
            metadata_joins.append(
                sql.SQL(
                    """
                    INNER JOIN {metadata_table_name} {metadata_table_name}
                        ON q.{field} = {metadata_table_name}."id"
                    """
                ).format(
                    metadata_table_name=sql.Identifier("metadata_" + field),
                    field=sql.Identifier(field),
                )
            )

        for snp_field in ["dna", "gene_aa", "protein_aa"]:
            if snp_field not in selected_fields:
                continue

            sequence_cols.append(snp_field + "_snp")
            sequence_cols_expr.append(
                sql.SQL("qq.{}").format(sql.Identifier(snp_field + "_snp"))
            )

        # CTEs and evaluating the metadata joins separately
        # from the SNV joins speeds this up by a lot
        snv_agg_field = "snp_str"
        if snv_format == constants["SNV_FORMAT"]["POS_REF_ALT"]:
            snv_agg_field = "snp_str"
        elif snv_format == constants["SNV_FORMAT"]["REF_POS_ALT"]:
            snv_agg_field = "snv_name"

        query = sql.SQL(
            """
            WITH dss AS (
                SELECT
                    q."id" as "sequence_id",
                    array_to_string(array_agg(ds.{snv_agg_field}), ';') as "snp"
                FROM {temp_table_name} q
                INNER JOIN "sequence_dna_snp" sds ON q."id" = sds."sequence_id"
                INNER JOIN "dna_snp" ds ON sds."snp_id" = ds."id"
                GROUP BY q."id"
            ),
            gass AS (
                SELECT
                    q."id" as "sequence_id",
                    array_to_string(array_agg(gas.{snv_agg_field}), ';') as "snp"
                FROM {temp_table_name} q
                INNER JOIN "sequence_gene_aa_snp" sgas ON q."id" = sgas."sequence_id"
                INNER JOIN "gene_aa_snp" gas ON sgas."snp_id" = gas."id"
                GROUP BY q."id"
            ),
            pass AS (
                SELECT
                    q."id" as "sequence_id",
                    array_to_string(array_agg(pas.{snv_agg_field}), ';') as "snp"
                FROM {temp_table_name} q
                INNER JOIN "sequence_protein_aa_snp" spas ON q."id" = spas."sequence_id"
                INNER JOIN "protein_aa_snp" pas ON spas."snp_id" = pas."id"
                GROUP BY q."id"
            ),
            qq AS (
                SELECT
                    q."id",
                    dss."snp" as "dna_snp",
                    gass."snp" as "gene_aa_snp",
                    pass."snp" as "protein_aa_snp"
                FROM {temp_table_name} q
                INNER JOIN dss ON q."id" = dss."sequence_id"
                INNER JOIN gass ON q."id" = gass."sequence_id"
                INNER JOIN pass ON q."id" = pass."sequence_id"
            )
            SELECT
                {sequence_cols_expr}
            FROM {temp_table_name} q
            INNER JOIN "location" loc ON q."location_id" = loc."id"
            {metadata_joins}
            JOIN qq ON qq."id" = q."id"
            """
        ).format(
            snv_agg_field=sql.Identifier(snv_agg_field),
            temp_table_name=sql.Identifier(temp_table_name),
            sequence_cols_expr=sql.SQL(",").join(sequence_cols_expr),
            metadata_joins=sql.SQL("\n").join(metadata_joins),
        )
        # print(query)
        cur.execute(query)

        res_df = pd.DataFrame.from_records(cur.fetchall(), columns=sequence_cols,)

    return make_response(res_df.to_csv(index=False), 200, {"Content-Type": "text/csv"})
