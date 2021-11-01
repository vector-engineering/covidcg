# coding: utf-8

"""Generate a mutation and lineage report

Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""

import io
import pandas as pd
import numpy as np
import psycopg2

from psycopg2 import sql
from flask import send_file


def generate_report(conn, req):
    """Generate a Spike mutation and lineage report
    This report will consist of:
    1) Single Spike SNV frequencies, both global and regional
    2) Co-occurring Spike SNV frequencies, both global and regional
    3) PANGO lineage frequencies, both global and regional

    Fetch the more complex data (regional data) from the database,
    and then from this compute the global counts in python.
    Format the regional tables as pivot tables, and do this in python,
    since crosstab() in Postgres makes me sad

    Start and end dates are passed via. URL query params
    Download the full excel file with:
    curl http://[host]/az_report?start_date=YYYY-MM-DD&end_date=YYYY-MM-DD -o out.xlsx
    Dates are in ISO format "YYYY-MM-DD"

    Parameters
    ----------
    conn: psycopg2.connection
    req: flask.request

    Returns
    -------
    Flask send_file response with a xlsx file attachment
    """

    start_date = pd.to_datetime(req.get("start_date", None))
    end_date = pd.to_datetime(req.get("end_date", None))

    with conn.cursor() as cur:

        # First, we'll need the total number of sequences in the database
        # We can infer it from one of the tables, but might as well just
        # fetch it explicitly
        cur.execute(
            """
        SELECT COUNT(*) 
        FROM "metadata"
        WHERE 
            "collection_date" >= %(start_date)s AND
            "collection_date" <= %(end_date)s
        """,
            {"start_date": start_date, "end_date": end_date},
        )
        num_seqs = cur.fetchone()[0]

        # REGIONAL SINGLE SPIKE SNVs
        cur.execute(
            """
            WITH mutation_region_counts AS (
                SELECT
                    "region", "mutation_id", COUNT(*) as "count"
                FROM (
                    SELECT "region", UNNEST("mutations") as "mutation_id"
                    FROM "sequence_gene_aa_snp" seq_snp
                    WHERE
                        "collection_date" >= %(start_date)s AND
                        "collection_date" <= %(end_date)s
                ) muts
                GROUP BY "region", "mutation_id"
            ),
            region_counts AS (
                SELECT "region", COUNT("sequence_id")
                FROM "metadata" m
                WHERE
                    "collection_date" >= %(start_date)s AND
                    "collection_date" <= %(end_date)s
                GROUP BY "region"
            )
            SELECT
                SUBSTRING(snp_def."snv_name" FROM 3) AS "name",
                snp_def."pos",
                snp_def."ref",
                snp_def."alt",
                mr."value",
                mutation_region_counts."count",
                ((mutation_region_counts."count"::REAL / region_counts."count"::REAL) * 100) AS "percent"
            FROM mutation_region_counts
            INNER JOIN "gene_aa_snp" snp_def ON mutation_region_counts."mutation_id" = snp_def."id"
            INNER JOIN region_counts ON region_counts."region" = mutation_region_counts."region"
            INNER JOIN "metadata_region" mr ON region_counts."region" = mr."id"
            WHERE snp_def."gene" = 'S'
            """,
            {"start_date": start_date, "end_date": end_date},
        )
        single_spike_snv_region = pd.DataFrame.from_records(
            cur.fetchall(),
            columns=["snv", "pos", "ref", "alt", "region", "count", "percent"],
        )
        # Pivot
        single_spike_snv_region_pivot = pd.pivot_table(
            single_spike_snv_region,
            index=["snv", "pos", "ref", "alt"],
            values=["count", "percent"],
            columns=["region"],
        ).fillna(0)
        # Collapse column multi-index
        single_spike_snv_region_pivot.columns = [
            "{}_{}".format(a, b)
            for a, b in zip(
                single_spike_snv_region_pivot.columns.get_level_values(0),
                single_spike_snv_region_pivot.columns.get_level_values(1),
            )
        ]
        # Cast count columns to integers
        count_cols = [
            col for col in single_spike_snv_region_pivot.columns if "count" in col
        ]
        for col in count_cols:
            single_spike_snv_region_pivot.loc[:, col] = single_spike_snv_region_pivot[
                col
            ].astype(int)
        single_spike_snv_region_pivot = single_spike_snv_region_pivot.reset_index()

        # Compute sum counts and sort on it
        single_spike_snv_region_pivot.insert(
            4, "sum_counts", single_spike_snv_region_pivot[count_cols].sum(axis=1)
        )
        single_spike_snv_region_pivot = single_spike_snv_region_pivot.sort_values(
            "sum_counts", ascending=False
        )
        # print(single_spike_snv_region_pivot)

        # GLOBAL SPIKE SINGLE SNVs
        single_spike_snv_global = (
            single_spike_snv_region.groupby("snv")
            .agg(
                pos=("pos", "first"),
                ref=("ref", "first"),
                alt=("alt", "first"),
                count=("count", np.sum),
            )
            .sort_values("count", ascending=False)
            .assign(percent=lambda x: x["count"] / num_seqs)
            .reset_index()
        )
        # print(single_spike_snv_global)

        # REGIONAL SPIKE COOC SNVs

        # First get map of IDs -> SNV names
        cur.execute(
            """
            SELECT "id", SUBSTRING("snv_name" FROM 3) AS "name", "pos"
            FROM gene_aa_snp
            WHERE gene = 'S'
            """
        )
        spike_snv_props = pd.DataFrame.from_records(
            cur.fetchall(), columns=["id", "name", "pos"]
        )
        spike_snv_name_map = dict(
            zip(spike_snv_props["id"].values, spike_snv_props["name"].values)
        )
        spike_snv_pos_map = dict(
            zip(spike_snv_props["id"].values, spike_snv_props["pos"].values)
        )

        cur.execute(
            """
            WITH seq_cooc AS (
                SELECT 
                    seq_snp."sequence_id",
                    seq_snp."region",
                    ("mutations" & (
                        SELECT ARRAY_AGG("id")
                        FROM "gene_aa_snp"
                        WHERE "gene" = 'S'
                    )) as "mutations",
                    m."lineage"
                FROM "sequence_gene_aa_snp" seq_snp
                INNER JOIN "metadata" m ON seq_snp."sequence_id" = m."sequence_id"
                WHERE
                    seq_snp."collection_date" >= %(start_date)s AND
                    seq_snp."collection_date" <= %(end_date)s
            ),
            most_common_lineage AS (
                SELECT DISTINCT ON ("mutations")
                    "mutations",
                    "count",
                    "lineage"
                FROM (
                    SELECT 
                        "mutations",
                        "lineage",
                        COUNT("sequence_id") as "count"
                    FROM seq_cooc
                    GROUP BY "mutations", "lineage"
                ) cooc_lineage
                ORDER BY "mutations", "count" DESC
            ),
            region_cooc_counts AS (
                SELECT
                    seq_cooc."region",
                    "mutations",
                    COUNT(seq_cooc."sequence_id") AS "count"
                FROM seq_cooc
                GROUP BY seq_cooc."mutations", seq_cooc."region"
            ),
            region_counts AS (
                SELECT 
                    seq_cooc."region",
                    COUNT(seq_cooc."sequence_id") as "count"
                FROM seq_cooc
                GROUP BY seq_cooc."region"
            )
            SELECT
                mr."value" as "region",
                most_common_lineage."lineage",
                region_cooc_counts."mutations",
                region_cooc_counts."count",
                (region_cooc_counts."count"::REAL / region_counts."count"::REAL) * 100 AS "percent"
            FROM region_cooc_counts
            INNER JOIN region_counts ON region_cooc_counts."region" = region_counts."region"
            INNER JOIN most_common_lineage ON region_cooc_counts."mutations" = most_common_lineage."mutations"
            INNER JOIN metadata_region mr ON region_cooc_counts."region" = mr."id"
            ORDER BY region_cooc_counts."count" DESC
            """,
            {"start_date": start_date, "end_date": end_date},
        )
        cooc_spike_snv_region = pd.DataFrame.from_records(
            cur.fetchall(),
            columns=["region", "lineage", "mutations", "count", "percent"],
        )

        # Sort mutations by position
        cooc_spike_snv_region.loc[:, "mutations"] = cooc_spike_snv_region[
            "mutations"
        ].apply(lambda x: sorted(x, key=lambda _x: spike_snv_pos_map[_x]))

        # Serialize list of mutation IDs
        cooc_spike_snv_region.loc[:, "mutations"] = cooc_spike_snv_region[
            "mutations"
        ].apply(lambda x: ";".join([spike_snv_name_map[_x] for _x in x]))

        # Pivot
        cooc_spike_snv_region_pivot = pd.pivot_table(
            cooc_spike_snv_region,
            index=["mutations", "lineage"],
            values=["count", "percent"],
            columns=["region"],
        ).fillna(0)
        # Collapse column multi-index
        cooc_spike_snv_region_pivot.columns = [
            "{}_{}".format(a, b)
            for a, b in zip(
                cooc_spike_snv_region_pivot.columns.get_level_values(0),
                cooc_spike_snv_region_pivot.columns.get_level_values(1),
            )
        ]
        # Cast count columns to integers
        count_cols = [
            col for col in cooc_spike_snv_region_pivot.columns if "count" in col
        ]
        for col in count_cols:
            cooc_spike_snv_region_pivot.loc[:, col] = cooc_spike_snv_region_pivot[
                col
            ].astype(int)
        cooc_spike_snv_region_pivot = cooc_spike_snv_region_pivot.reset_index()

        # Compute sum counts and sort on it
        cooc_spike_snv_region_pivot.insert(
            2, "sum_counts", cooc_spike_snv_region_pivot[count_cols].sum(axis=1)
        )
        cooc_spike_snv_region_pivot = cooc_spike_snv_region_pivot.sort_values(
            "sum_counts", ascending=False
        )
        # print(cooc_spike_snv_region_pivot)

        # GLOBAL SPIKE COOC SNVs
        cooc_spike_snv_global = (
            cooc_spike_snv_region.groupby("mutations")
            .agg(lineage=("lineage", "first"), count=("count", np.sum),)
            .sort_values("count", ascending=False)
            .assign(percent=lambda x: x["count"] / num_seqs)
            .reset_index()
        )
        # print(cooc_spike_snv_global)

        # REGIONAL LINEAGE COUNTS
        cur.execute(
            """
            WITH region_counts AS (
                SELECT 
                    m."region",
                    COUNT(m."sequence_id") as "count"
                FROM "metadata" m
                WHERE
                    m."collection_date" >= %(start_date)s AND
                    m."collection_date" <= %(end_date)s
                GROUP BY m."region"
            ),
            group_counts AS (
                SELECT
                    m."region",
                    m."lineage",
                    COUNT(m."sequence_id") AS "count"
                FROM "metadata" m
                WHERE
                    m."collection_date" >= %(start_date)s AND
                    m."collection_date" <= %(end_date)s
                GROUP BY m."region", m."lineage"
            )
            SELECT
                mr."value" as "region",
                g."lineage",
                g."count",
                (g."count"::REAL / region_counts."count"::REAL) * 100 AS "percent"
            FROM group_counts g
            INNER JOIN region_counts ON g."region" = region_counts."region"
            INNER JOIN metadata_region mr ON g."region" = mr."id"
            """,
            {"start_date": start_date, "end_date": end_date},
        )
        lineage_region = pd.DataFrame.from_records(
            cur.fetchall(), columns=["region", "lineage", "count", "percent"]
        )
        # Pivot
        lineage_region_pivot = pd.pivot_table(
            lineage_region,
            index=["lineage"],
            values=["count", "percent"],
            columns=["region"],
        ).fillna(0)
        # Collapse column multi-index
        lineage_region_pivot.columns = [
            "{}_{}".format(a, b)
            for a, b in zip(
                lineage_region_pivot.columns.get_level_values(0),
                lineage_region_pivot.columns.get_level_values(1),
            )
        ]
        # Cast count columns to integers
        count_cols = [col for col in lineage_region_pivot.columns if "count" in col]
        for col in count_cols:
            lineage_region_pivot.loc[:, col] = lineage_region_pivot[col].astype(int)
        lineage_region_pivot = lineage_region_pivot.reset_index()

        # Compute sum counts and sort on it
        lineage_region_pivot.insert(
            1, "sum_counts", lineage_region_pivot[count_cols].sum(axis=1)
        )
        lineage_region_pivot = lineage_region_pivot.sort_values(
            "sum_counts", ascending=False
        )

        # print(lineage_region_pivot)

        # GLOBAL LINEAGES
        lineage_global = (
            lineage_region.groupby("lineage")
            .agg(count=("count", np.sum))
            .sort_values("count", ascending=False)
            .assign(percent=lambda x: x["count"] / num_seqs)
            .reset_index()
        )
        # print(lineage_global)

    buf = io.BytesIO()
    with pd.ExcelWriter(buf, engine="openpyxl") as writer:
        single_spike_snv_global.to_excel(writer, index=False, sheet_name="SAV_Global")
        single_spike_snv_region_pivot.to_excel(
            writer, index=False, sheet_name="SAV_Regional"
        )
        cooc_spike_snv_global.to_excel(writer, index=False, sheet_name="Cooc_Global")
        cooc_spike_snv_region_pivot.to_excel(
            writer, index=False, sheet_name="Cooc_Regional"
        )
        lineage_global.to_excel(writer, index=False, sheet_name="Lineage_Global")
        lineage_region_pivot.to_excel(
            writer, index=False, sheet_name="Lineage_Regional"
        )
    buf.seek(0)

    return send_file(
        buf,
        mimetype="application/octet-stream",
        as_attachment=True,
        attachment_filename="test.xlsx",
    )

