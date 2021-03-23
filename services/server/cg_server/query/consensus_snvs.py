# coding: utf-8

"""Get consensus SNVs from a lineage/clade

Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""

import pandas as pd
import psycopg2

from psycopg2 import sql

from cg_server.query.selection import build_coordinate_filters


def query_consensus_snvs(conn, req, groups):

    group_key = req.get("group_key", None)
    dna_or_aa = req.get("dna_or_aa", None)
    coordinate_mode = req.get("coordinate_mode", None)
    coordinate_ranges = req.get("coordinate_ranges", None)
    selected_gene = req.get("selected_gene", None)
    selected_protein = req.get("selected_protein", None)

    with conn.cursor() as cur:

        (
            snv_cols,
            snv_filter,
            snv_table,
            pos_filter_injections,
        ) = build_coordinate_filters(
            conn, dna_or_aa, coordinate_mode, coordinate_ranges,
        )
        if snv_filter.as_string(conn):
            snv_filter = sql.Composed([sql.SQL(" AND "), snv_filter])

        snv_cols_selection = sql.SQL(",\n").join(
            [sql.SQL("snp_data.{}").format(sql.Identifier(col)) for col in snv_cols]
        )

        consensus_table = "{group_key}_consensus_{snv_table}".format(
            group_key=group_key, snv_table=snv_table
        )

        cur.execute(
            sql.SQL(
                """
                SELECT 
                    consensus."name",
                    {snv_cols_selection}
                FROM {consensus_table} consensus
                JOIN {snv_table} snp_data ON snp_data.id = consensus.snp_id
                WHERE 
                    consensus."name" IN %(groups)s
                    {snv_filter}
                ORDER BY snp_data."pos" ASC
                """
            ).format(
                snv_cols_selection=snv_cols_selection,
                consensus_table=sql.Identifier(consensus_table),
                snv_table=sql.Identifier(snv_table),
                snv_filter=snv_filter,
            ),
            {
                "groups": tuple(groups),
                "selected_gene": selected_gene,
                "selected_protein": selected_protein,
                **pos_filter_injections,
            },
        )

        group_snvs = pd.DataFrame.from_records(
            cur.fetchall(), columns=["group"] + snv_cols,
        )

    return group_snvs

