# coding: utf-8

"""Generate a variant table report

Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""

import io
import pandas as pd

from flask import send_file
from cg_server.config import config
from cg_server.constants import constants
from cg_server.query import (
    build_sequence_location_where_filter,
    build_coordinate_filters,
    get_loc_level_ids,
)
from psycopg2 import sql


def build_variant_table(conn, req):
    """Generate a variant table report

    Args:
        conn (psycopg2.extensions.connection): Database connection
        req (dict): Request object

    Returns:
        pd.DataFrame: Variant table report
    """

    dna_or_aa = req.get("dna_or_aa", None)
    coordinate_mode = req.get("coordinate_mode", None)
    coordinate_ranges = req.get("coordinate_ranges", None)
    selected_gene = req.get("selected_gene", None)
    selected_protein = req.get("selected_protein", None)
    selected_reference = req.get("selected_reference", None)

    if not selected_reference:
        raise Exception("No reference specified")

    with conn.cursor() as cur:
        sequence_where_filter = build_sequence_location_where_filter(
            constants["GROUP_MUTATION"],
            get_loc_level_ids(req),
            req.get("start_date", None),
            req.get("end_date", None),
            req.get("subm_start_date", None),
            req.get("subm_end_date", None),
            req.get("selected_metadata_fields", None),
            req.get("selected_group_fields", None),
            selected_reference,
            req.get("sequence_length", None),
            req.get("percent_ambiguous", None),
        )
        

        (mutation_filter, mutation_table) = build_coordinate_filters(
            conn,
            dna_or_aa,
            coordinate_mode,
            coordinate_ranges,
            selected_gene,
            selected_protein,
        )
        sequence_mutation_table = "sequence_" + mutation_table

        # Fields that the user wants
        selected_fields = req.get("selected_fields", [])
        mutation_format = req.get(
            "mutation_format", constants["MUTATION_FORMAT"]["POS_REF_ALT"]
        )

        # Get grouping columns, metadata columns
        metadata_cols = [
            "Accession ID",
            "collection_date",
            "submission_date",
        ]
        metadata_cols_expr = [
            sql.SQL("m.{}").format(sql.Identifier(col)) for col in metadata_cols
        ]
        joins = []

        for grouping in config["group_cols"].keys():
            if grouping not in selected_fields:
                continue

            metadata_cols.append(grouping)
            metadata_cols_expr.append(sql.SQL("m.{}").format(sql.Identifier(grouping)))

        for field in list(config["metadata_cols"].keys()) + list(
            constants["GEO_LEVELS"].values()
        ):
            if field not in selected_fields:
                continue

            metadata_cols.append(field)
            metadata_cols_expr.append(
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

        query = sql.SQL(
            """
            SELECT
                {metadata_cols_expr},
                sm."reference",
                sm."mutations"
            FROM (
                SELECT
                    sst."sequence_id",
                    sst."reference",
                    (sst.mutations & (
                        SELECT ARRAY_AGG("id")
                        FROM {mutation_table}
                        {mutation_filter}
                    )) as "mutations"
                FROM {sequence_mutation_table} sst
                WHERE {sequence_where_filter}
            ) sm
            INNER JOIN "metadata" m ON sm."sequence_id" = m."sequence_id"
            {joins}
            ORDER BY m."Accession ID" ASC
            """
        ).format(
            metadata_cols_expr=sql.SQL(",").join(metadata_cols_expr),
            sequence_mutation_table=sql.Identifier(sequence_mutation_table),
            sequence_where_filter=sequence_where_filter,
            mutation_table=sql.Identifier(mutation_table),
            joins=sql.SQL("\n").join(joins),
            mutation_filter=mutation_filter,
        )
        # print(query.as_string(conn))
        cur.execute(query)

        res = pd.DataFrame.from_records(
            # cur.fetchall(), columns=metadata_cols + ["mutation_name", "pos"]
            cur.fetchall(),
            columns=metadata_cols + ["reference", "mutations"],
        ).set_index("Accession ID")

        mutation_name_field = "mutation_str"
        if mutation_format == constants["MUTATION_FORMAT"]["POS_REF_ALT"]:
            mutation_name_field = "mutation_str"
        elif mutation_format == constants["MUTATION_FORMAT"]["REF_POS_ALT"]:
            mutation_name_field = "mutation_name"

        # Get mutation ID mappings
        cur.execute(
            sql.SQL(
                """
            SELECT 
                "id",
                SUBSTRING({mutation_name_field} FROM 3) AS "name",
                "pos"
            FROM {mutation_table}
            {mutation_filter}
            """
            ).format(
                mutation_name_field=sql.Identifier(mutation_name_field),
                mutation_table=sql.Identifier(mutation_table),
                mutation_filter=mutation_filter,
            )
        )
        mutation_def = (
            pd.DataFrame.from_records(cur.fetchall(), columns=["id", "name", "pos"])
            .set_index("id")
            .sort_values("pos", ascending=True)
        )

        # print(mutation_def)

        # Create pivot table of mutations
        # Goal: each row = one sequence
        # Also include the metadata (dates, etc)
        # Each column after metadata = one mutation
        # 0 = doesn't have it, 1 = has it

        mutation_id_to_name = mutation_def["name"].to_dict()
        mutation_id_to_pos = mutation_def["pos"].to_dict()

        pivot = (
            pd.pivot_table(
                (
                    res.assign(val=lambda _: 1)
                    .explode("mutations")
                    .assign(
                        mutations=lambda x: x["mutations"]
                        .map(mutation_id_to_name)
                        .fillna("N/A")
                    )
                ),
                index=["Accession ID"],
                columns="mutations",
                values="val",
            )
            .fillna(0)
            .astype(int)
        )

        # print(pivot)

        # Reorder columns by position
        mut_columns_ordered = [n for n in mutation_def["name"] if n in pivot.columns]
        pivot = pivot[mut_columns_ordered]

        # Create column of all mutations concatenated
        def f_mutation_id_to_pos(x):
            return mutation_id_to_pos[x]

        all_mutations = res["mutations"].apply(
            lambda x: ",".join(
                [mutation_id_to_name[_x] for _x in sorted(x, key=f_mutation_id_to_pos)]
            )
        )
        pivot.insert(0, "all_mutations", all_mutations)

        # Join metadata back onto the pivot table
        pivot = res.drop(columns=["mutations"]).join(pivot, how="inner")

        # Create another table with aggregated mutation counts
        # Each row = one mutation
        mut_counts = (
            pivot[mut_columns_ordered]
            .sum(axis=0)
            .rename("counts")
            .to_frame()
            .merge(mutation_def, left_index=True, right_on="name", how="left")
            .reset_index(drop=True)
            .sort_values("pos", ascending=True)
            .drop(columns=["pos"])[["name", "counts"]]  # Reorder columns
            .rename(columns={"name": "mutation_name"})
        )
        mut_counts.insert(0, "reference", selected_reference)

        # print(mut_counts)
        # print(pivot)

    buf = io.BytesIO()
    with pd.ExcelWriter(buf, engine="openpyxl") as writer:
        pivot.to_excel(writer, sheet_name="sequence_mutations")
        mut_counts.to_excel(writer, index=False, sheet_name="mutation_counts")
    buf.seek(0)

    return send_file(
        buf,
        mimetype="application/octet-stream",
        as_attachment=True,
        attachment_filename="variant_table.xlsx",
    )
