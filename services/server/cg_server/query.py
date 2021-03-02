# coding: utf-8

import pandas as pd
import uuid

from cg_server.config import config
from cg_server.constants import constants


def build_coordinate_filters(
    dna_or_aa, coordinate_mode, coordinate_ranges, selected_gene, selected_protein
):
    snv_cols = ["snp_str", "snv_name", "color", "pos", "ref", "alt"]
    snv_filter = []
    if dna_or_aa == constants["DNA_OR_AA"]["DNA"]:
        snv_table = "dna_snp"
    elif dna_or_aa == constants["DNA_OR_AA"]["AA"]:
        snv_cols.append("nt_pos")
        if coordinate_mode == constants["COORDINATE_MODES"]["COORD_GENE"]:
            snv_table = "gene_aa_snp"
            snv_cols.append("gene")
            snv_filter.append("snp_data.\"gene\" = '{}'".format(selected_gene))
        elif coordinate_mode == constants["COORDINATE_MODES"]["COORD_PROTEIN"]:
            snv_table = "protein_aa_snp"
            snv_cols.append("protein")
            snv_filter.append("snp_data.\"protein\" = '{}'".format(selected_protein))

    for start, end in coordinate_ranges:
        pos_column = "pos" if dna_or_aa == constants["DNA_OR_AA"]["DNA"] else "nt_pos"
        snv_filter.append('snp_data."{}" >= {}'.format(pos_column, start))
        snv_filter.append('snp_data."{}" <= {}'.format(pos_column, end))

    snv_filter = " AND ".join(snv_filter)

    return snv_cols, snv_filter, snv_table


def select_sequences(cur, req):

    group_key = req.get("group_key", None)
    dna_or_aa = req.get("dna_or_aa", None)
    coordinate_mode = req.get("coordinate_mode", None)

    location_ids = req.get("location_ids", None)
    all_location_ids = sum(location_ids.values(), [])

    start_date = pd.to_datetime(req.get("start_date", None))
    end_date = pd.to_datetime(req.get("end_date", None))

    # Metadata filters will come in the form of a JSON object
    # of { metadata_field: [metadata_values] }
    selected_metadata_fields = req.get("selected_metadata_fields", None)
    metadata_filters = []
    for md_key, md_vals in selected_metadata_fields.items():
        metadata_filters.append(
            '"{field}" IN ({vals})'.format(
                field=md_key, vals=",".join([str(val) for val in md_vals])
            )
        )
    metadata_filters = " AND ".join(metadata_filters)
    if metadata_filters:
        metadata_filters += " AND "

    temp_table_name = "query_" + uuid.uuid4().hex

    cur.execute(
        """
        CREATE TEMP TABLE "{temp_table_name}"
        ON COMMIT DROP
        AS (
            SELECT seq.*
            FROM "sequence" seq
            WHERE
                "collection_date" >= %(start_date)s AND
                "collection_date" <= %(end_date)s AND
                {metadata_filters}
                "location_id" IN %(location_ids)s
        );
        """.format(
            temp_table_name=temp_table_name, metadata_filters=metadata_filters
        ),
        {
            "start_date": start_date,
            "end_date": end_date,
            "location_ids": tuple(all_location_ids),
        },
    )

    return temp_table_name


def query_sequences(conn, req):

    with conn.cursor() as cur:

        temp_table_name = select_sequences(cur, req)

        group_key = req.get("group_key", None)
        dna_or_aa = req.get("dna_or_aa", None)
        coordinate_mode = req.get("coordinate_mode", None)
        coordinate_ranges = req.get("coordinate_ranges", None)
        selected_gene = req.get("selected_gene", None)
        selected_protein = req.get("selected_protein", None)

        snv_cols, snv_filter, snv_table = build_coordinate_filters(
            dna_or_aa,
            coordinate_mode,
            coordinate_ranges,
            selected_gene,
            selected_protein,
        )
        sequence_snv_table = "sequence_" + snv_table

        snv_cols_selection = ",\n".join(
            ['snp_data."{}"'.format(col) for col in snv_cols]
        )

        cur.execute(
            """
            SELECT
                snp."sequence_id", 
                snp."snp_id",
                seq."Accession ID",
                seq."location_id",
                seq."collection_date",
                {snv_cols_selection}
            FROM "{sequence_snv_table}" snp
            RIGHT JOIN "{temp_table_name}" seq ON seq.id = snp.sequence_id 
            INNER JOIN "{snv_table}" snp_data ON snp.snp_id = snp_data.id
            WHERE {snv_filter};
            """.format(
                snv_cols_selection=snv_cols_selection,
                sequence_snv_table=sequence_snv_table,
                temp_table_name=temp_table_name,
                snv_table=snv_table,
                snv_filter=snv_filter,
            )
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
        sequence_cols_expr = ",".join(['q."{}"'.format(col) for col in sequence_cols])

        join_expr = ""
        if group_key != constants["GROUP_SNV"]:
            sequence_cols.append("color")
            sequence_cols_expr += ',group_data."color"'
            join_expr = 'JOIN {group_key} group_data ON q."{group_key}" = group_data."name"'.format(
                group_key=group_key
            )

        cur.execute(
            """
            SELECT {sequence_cols} 
            FROM "{temp_table_name}" q
            {join_expr};
            """.format(
                sequence_cols=sequence_cols_expr,
                temp_table_name=temp_table_name,
                join_expr=join_expr,
            )
        )
        res_df = pd.DataFrame.from_records(
            cur.fetchall(), index="id", columns=sequence_cols
        )

        # Clean up
        cur.execute(
            'DROP TABLE IF EXISTS "{temp_table_name}";'.format(
                temp_table_name=temp_table_name
            )
        )

    return res_df, res_snv


def query_consensus_snvs(conn, req, groups):

    group_key = req.get("group_key", None)
    dna_or_aa = req.get("dna_or_aa", None)
    coordinate_mode = req.get("coordinate_mode", None)
    coordinate_ranges = req.get("coordinate_ranges", None)
    selected_gene = req.get("selected_gene", None)
    selected_protein = req.get("selected_protein", None)

    with conn.cursor() as cur:

        snv_cols, snv_filter, snv_table = build_coordinate_filters(
            dna_or_aa,
            coordinate_mode,
            coordinate_ranges,
            selected_gene,
            selected_protein,
        )
        if snv_filter:
            snv_filter = " AND " + snv_filter
        snv_cols_selection = ",\n".join(
            ['snp_data."{}"'.format(col) for col in snv_cols]
        )

        consensus_table = "{group_key}_consensus_{snv_table}".format(
            group_key=group_key, snv_table=snv_table
        )

        cur.execute(
            """
            SELECT 
                consensus."name",
                {snv_cols_selection}
            FROM "{consensus_table}" consensus
            JOIN "{snv_table}" snp_data ON snp_data.id = consensus.snp_id
            WHERE 
                consensus."name" IN %(groups)s
                {snv_filter}
            ORDER BY snp_data."pos" ASC
            """.format(
                snv_cols_selection=snv_cols_selection,
                consensus_table=consensus_table,
                snv_table=snv_table,
                snv_filter=snv_filter,
            ),
            {"groups": tuple(groups)},
        )

        group_snvs = pd.DataFrame.from_records(
            cur.fetchall(), columns=["group"] + snv_cols,
        )

    return group_snvs

