# coding: utf-8

"""Get sequences/SNVs from database based on user selections

Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""

import pandas as pd
import psycopg2
import uuid

from psycopg2 import sql
from psycopg2.extras import execute_values

from cg_server.constants import constants


def build_coordinate_filters(
    conn, dna_or_aa, coordinate_mode, coordinate_ranges, selected_gene, selected_protein
):
    """Build SQL 'WHERE' expression to filter for SNVs within the user-defined range
    Only valid in "SNV" mode, not in any other grouping mode

    Parameters
    ----------
    conn: psycopg2.connection
    dna_or_aa: str
        - Enum of constants["DNA_OR_AA"] (DNA mode or AA mode)
    coordinate_mode: str
        - Enum of constants["COORDINATE_MODE"] (NT, GENE AA, or PROTEIN AA)
    coordinate_ranges: list of [int, int]
        - List of ranges, in nucleotide coordinates, for valid SNVs
    selected_gene: str
        - Only used for filtering if in GENE AA mode
    selected_protein: str
        - Only used for filtering if in PROTEIN AA mode

    Returns
    -------
    snv_filter: sql.Composed
        - The "WHERE" clause, in a format-ready sql.Composable form
    snv_table
        - The name of the SNV definition table to use (dna_snp, gene_aa_snp, protein_aa_snp)

    """

    snv_filter = []
    if dna_or_aa == constants["DNA_OR_AA"]["DNA"]:
        snv_table = "dna_snp"
    elif dna_or_aa == constants["DNA_OR_AA"]["AA"]:
        if coordinate_mode == constants["COORDINATE_MODES"]["COORD_GENE"]:
            snv_table = "gene_aa_snp"
            snv_filter.append(
                sql.SQL('"gene" = {gene}').format(gene=sql.Literal(selected_gene))
            )
        elif coordinate_mode == constants["COORDINATE_MODES"]["COORD_PROTEIN"]:
            snv_table = "protein_aa_snp"
            snv_filter.append(
                sql.SQL('"protein" = {protein}').format(
                    protein=sql.Literal(selected_protein)
                )
            )

    snv_filter = sql.SQL(" AND ").join(snv_filter)

    pos_filter = []
    for i in range(len(coordinate_ranges)):
        pos_column = "pos" if dna_or_aa == constants["DNA_OR_AA"]["DNA"] else "nt_pos"
        pos_filter.append(
            sql.SQL(
                """
                ({pos_column} >= {start} AND {pos_column} <= {end})
                """
            ).format(
                pos_column=sql.Identifier(pos_column),
                start=sql.Literal(coordinate_ranges[i][0]),
                end=sql.Literal(coordinate_ranges[i][1]),
            )
        )

    pos_filter = sql.SQL(" OR ").join(pos_filter)

    # Compose final WHERE expression
    snv_filter = [snv_filter]
    # Only combine the snv_filter and pos_filter if the snv_filter exists
    if snv_filter[0].as_string(conn):
        snv_filter.append(sql.SQL(" AND "))

    snv_filter.append(pos_filter)
    snv_filter = sql.Composed(snv_filter)

    # Only add the WHERE clause if the snv_filter exists
    if snv_filter.join("").as_string(conn):
        snv_filter = (sql.SQL("WHERE") + snv_filter).join(" ")

    return snv_filter, snv_table


def create_location_map_table(cur, location_ids):
    """Create a temp table of location name -> location_id
    This is joined onto the results so that we pass to the client
    an aggregated list of group counts per location name

    Parameters
    ----------
    cur: psycopg2.cursor
    location_ids: dict
        - Structured as { "location_name": [location_ids], ... }
          i.e., keys are location names as arbitrary strings, and
          location_ids is a list of integers

    Returns
    -------
    location_map_table_name: str
        - Name of the location map temp table
          Table name has a uuid appended to not collide with other
          ongoing queries
    """

    # Transform into a list of tuples for record insertion
    # [('location', id), ...]
    location_id_to_name_list = []
    for location in location_ids.keys():
        for location_id in location_ids[location]:
            location_id_to_name_list.append((location, location_id))

    location_map_table_name = "location_map_" + uuid.uuid4().hex

    cur.execute(
        sql.SQL(
            """
        CREATE TEMP TABLE {location_map_table_name} (
            "id" INTEGER NOT NULL,
            "name" TEXT  NOT NULL 
        ) ON COMMIT DROP;
        """
        ).format(location_map_table_name=sql.Identifier(location_map_table_name))
    )

    execute_values(
        cur,
        sql.SQL(
            """INSERT INTO {location_map_table_name} ("name", "id") VALUES %s"""
        ).format(location_map_table_name=sql.Identifier(location_map_table_name)),
        location_id_to_name_list,
    )

    # Write raw SQL for debugging purposes
    # print(
    #     """
    #     INSERT INTO "{location_map_table_name}" ("name", "id") VALUES
    #     {location_id_to_name_list}
    #     """.format(
    #         location_map_table_name=location_map_table_name,
    #         location_id_to_name_list=",\n".join(
    #             ["('{}', {})".format(x[0], x[1]) for x in location_id_to_name_list]
    #         ),
    #     )
    # )

    return location_map_table_name


def build_sequence_query(
    location_ids={},
    start_date=None,
    end_date=None,
    subm_start_date=None,
    subm_end_date=None,
    selected_metadata_fields={},
):
    """Build query for filtering sequences based on user's location/date
    selection and selected metadata fields

    Parameters
    ----------
    location_ids: dict
        - Structured as { "location_name": [location_ids], ... }
        - Keys are location names as arbitrary strings
        - Values (location_ids) are a list of integers
    start_date: str
        - Collection sart date, in ISO format (YYYY-MM-DD)
    end_date: str
        - Collection end date, in ISO format (YYYY-MM-DD)
    subm_start_date: str
        - Submission start date, in ISO format (YYYY-MM-DD)
    subm_end_date: str
        - Submission end date, in ISO format (YYYY-MM-DD)
    selected_metadata_fields: dict
        - Structured as { metadata_field: [metadata_values] }
        - Keys are a metadata field, as a string
        - Values are a list of metadata value IDs (integers)

    Returns
    -------
    sequence_query: sql.SQL
        - Ready for format-injection into a main query, as a
          inline table expression or CTE

    """

    # Combine all location IDs into one big array
    all_location_ids = sum(location_ids.values(), [])

    # Construct submission date filters
    if subm_start_date is None and subm_end_date is None:
        submission_date_filter = sql.SQL("")
    else:
        chunks = []
        if subm_start_date is not None:
            chunks.append(
                sql.SQL('"submission_date" >= {}').format(sql.Literal(subm_start_date))
            )
        if subm_end_date is not None:
            chunks.append(
                sql.SQL('"submission_date" <= {}').format(sql.Literal(subm_end_date))
            )

        submission_date_filter = sql.Composed([sql.SQL(" AND ").join(chunks), sql.SQL(" AND ")])

    metadata_filters = []
    # Build dictionary of metadata value tuples to inject
    for md_key, md_vals in selected_metadata_fields.items():
        # Don't process if the list of metadata values is empty
        if not md_vals:
            continue

        metadata_filters.append(
            sql.SQL("{field} IN {vals}").format(
                field=sql.Identifier(md_key),
                vals=sql.Literal(tuple([str(val) for val in md_vals])),
            )
        )

    if metadata_filters:
        metadata_filters = sql.SQL(" AND ").join(metadata_filters)
        metadata_filters = sql.Composed([metadata_filters, sql.SQL(" AND ")])
    else:
        metadata_filters = sql.SQL("")

    sequence_query = sql.SQL(
        """
        SELECT m.*
        FROM "metadata" m
        WHERE
            {metadata_filters}
            "collection_date" >= {start_date} AND
            "collection_date" <= {end_date} AND
            {submission_date_filter}
            "location_id" = ANY({location_ids})
        """
    ).format(
        metadata_filters=metadata_filters,
        start_date=sql.Literal(start_date),
        end_date=sql.Literal(end_date),
        submission_date_filter=submission_date_filter,
        location_ids=sql.Literal(all_location_ids),
    )

    return sequence_query


def create_sequence_temp_table(cur, req):
    """Build the sequence query, run it, and store the results
    in a temporary table

    Parameters
    ----------
    conn: psycopg2.cursor
    req: flask.request

    Returns
    -------
    temp_table_name: str
    """

    location_ids = req.get("location_ids", None)

    start_date = pd.to_datetime(req.get("start_date", None))
    end_date = pd.to_datetime(req.get("end_date", None))

    subm_start_date = req.get("subm_start_date", "")
    subm_end_date = req.get("subm_end_date", "")
    subm_start_date = None if subm_start_date == "" else pd.to_datetime(subm_start_date)
    subm_end_date = None if subm_end_date == "" else pd.to_datetime(subm_end_date)

    selected_metadata_fields = req.get("selected_metadata_fields", None)

    # First store the sequence query into a temp table
    sequence_query = build_sequence_query(
        location_ids=location_ids,
        start_date=start_date,
        end_date=end_date,
        subm_start_date=subm_start_date,
        subm_end_date=subm_end_date,
        selected_metadata_fields=selected_metadata_fields,
    )
    temp_table_name = "sequence_selection_" + uuid.uuid4().hex
    cur.execute(
        sql.SQL(
            """
        CREATE TEMP TABLE {temp_table_name}
        ON COMMIT DROP
        AS ({sequence_query})
        """
        ).format(
            temp_table_name=sql.Identifier(temp_table_name),
            sequence_query=sequence_query,
        )
    )

    return temp_table_name


def query_and_aggregate(conn, req):
    """Select sequences and aggregate results
    If in "SNV" mode, return sequences grouped by co-occurring SNVs
    In any other mode, return sequences grouped by their group column
    (i.e., lineage or clade)

    Sequences are also grouped by location (as defined in the
    location_id map), as well as collection date

    Parameters
    ----------
    conn: psycopg2.connection
    req: flask.request

    Returns
    -------
    res: pandas.DataFrame
        - Dataframe with 4 columns: "location", "collection_date", 
        "group_id", "counts"
        - Each row represents sequences aggregated by location, 
          collection_date, and group_id
        - In "SNV" mode, "group_id" represents co-occurring SNVs,
          and is structured as a list of SNV IDs (integers)
        - In any other mode, "group_id" is a string of the group
          name. i.e., in lineage grouping, group_id is the lineage
          name.
    """

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
        group_key = req.get("group_key", None)
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

        if group_key == constants["GROUP_SNV"]:
            (snv_filter, snv_table) = build_coordinate_filters(
                conn,
                dna_or_aa,
                coordinate_mode,
                coordinate_ranges,
                selected_gene,
                selected_protein,
            )
            sequence_snv_table = "sequence_" + snv_table

            # welcome to CTE hell
            main_query = sql.SQL(
                """
                WITH "seq" AS (
                    {sequence_query}
                ),
                "snp_data" AS (
                    SELECT "id"
                    FROM {snv_table}
                    {snv_filter}
                ),
                "filtered_snp" AS (
                    SELECT
                        "seq"."id" AS "sequence_id",
                        "seq"."collection_date",
                        loc_map."name" AS "location",
                        COALESCE(snp."snp_id", -1) AS "snp_id"
                    FROM {sequence_snv_table} snp
                    INNER JOIN "snp_data" ON snp_data.id = snp.snp_id
                    RIGHT OUTER JOIN "seq" ON snp.sequence_id = "seq".id
                    LEFT OUTER JOIN (
                        SELECT "id", "name"
                        FROM {location_map_table_name}
                    ) loc_map ON "seq"."location_id" = loc_map."id"
                ),
                "snv_list" AS (
                    SELECT 
                        "location",
                        "collection_date",
                        ARRAY_AGG("snp_id") AS "group_id"
                    FROM "filtered_snp"
                    GROUP BY "location", "collection_date", "sequence_id"
                )
                SELECT
                    "location",
                    "collection_date",
                    "group_id",
                    COUNT(*) as "count"
                FROM "snv_list"
                GROUP BY "location", "collection_date", "group_id"
                """
            ).format(
                sequence_query=sequence_query,
                location_map_table_name=sql.Identifier(location_map_table_name),
                sequence_snv_table=sql.Identifier(sequence_snv_table),
                snv_table=sql.Identifier(snv_table),
                snv_filter=snv_filter,
            )

            # print(main_query.as_string(conn))

        else:
            main_query = sql.SQL(
                """
                SELECT 
                    location_map."name",
                    q."collection_date",
                    q."lineage",
                    count(*) as "count"
                FROM (
                    {sequence_query}
                ) q
                LEFT OUTER JOIN (
                    SELECT "id", "name"
                    FROM {location_map_table_name}
                ) location_map ON q."location_id" = location_map."id"
                GROUP BY location_map."name", q."collection_date", q."lineage"
                """
            ).format(
                sequence_query=sequence_query,
                location_map_table_name=sql.Identifier(location_map_table_name),
            )

        # print(main_query.as_string(conn))
        cur.execute(main_query)

        res = pd.DataFrame.from_records(
            cur.fetchall(),
            columns=["location", "collection_date", "group_id", "counts"],
        ).to_json(orient="records")

    return res

