# coding: utf-8

"""Get sequences/mutations from database based on user selections

Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""

import datetime

import pandas as pd
from flask import jsonify
from psycopg2 import sql
from cg_server.config import config
from cg_server.constants import constants


def build_coordinate_filters(
    conn, dna_or_aa, coordinate_mode, coordinate_ranges, selected_gene, selected_protein
):
    """Build SQL 'WHERE' expression to filter for mutations within the user-defined range
    Only valid in "mutation" mode, not in any other grouping mode

    Parameters
    ----------
    conn: psycopg2.connection
    dna_or_aa: str
        - Enum of constants["DNA_OR_AA"] (DNA mode or AA mode)
    coordinate_mode: str
        - Enum of constants["COORDINATE_MODES"] (NT, GENE AA, or PROTEIN AA)
    coordinate_ranges: list of [segment (str), start (int), end (int)]
        - List of ranges, in nucleotide coordinates, for valid mutations
    selected_gene: str
        - Only used for filtering if in GENE AA mode
    selected_protein: str
        - Only used for filtering if in PROTEIN AA mode

    Returns
    -------
    mutation_filter: sql.Composed
        - The "WHERE" clause, in a format-ready sql.Composable form
    mutation_table
        - The name of the mutation definition table to use (dna_mutation, gene_aa_mutation, protein_aa_mutation)

    """

    mutation_filter = []
    if dna_or_aa == constants["DNA_OR_AA"]["DNA"]:
        mutation_table = "dna_mutation"
    elif dna_or_aa == constants["DNA_OR_AA"]["AA"]:
        if coordinate_mode == constants["COORDINATE_MODES"]["COORD_GENE"]:
            mutation_table = "gene_aa_mutation"
            mutation_filter.append(
                sql.SQL('"feature" = {gene}').format(gene=sql.Literal(selected_gene))
            )
        elif coordinate_mode == constants["COORDINATE_MODES"]["COORD_PROTEIN"]:
            mutation_table = "protein_aa_mutation"
            mutation_filter.append(
                sql.SQL('"feature" = {protein}').format(
                    protein=sql.Literal(selected_protein)
                )
            )

    mutation_filter = sql.SQL(" AND ").join(mutation_filter)

    pos_filter = []
    if coordinate_ranges:
        for coord_range in coordinate_ranges:
            range_filter = []
            # For DNA mode, also match the correct segment
            if dna_or_aa == constants["DNA_OR_AA"]["DNA"]:
                range_filter.append(
                    sql.SQL('"segment" = {}').format(sql.Literal(coord_range[0]))
                )
            range_filter.append(
                sql.SQL('"pos" >= {}').format(sql.Literal(coord_range[1]))
            )
            range_filter.append(
                sql.SQL('"pos" <= {}').format(sql.Literal(coord_range[2]))
            )

            pos_filter.append(
                sql.Composed(
                    [sql.SQL("("), sql.SQL(" AND ").join(range_filter), sql.SQL(")")]
                )
            )

    pos_filter = sql.SQL(" OR ").join(pos_filter)

    # Compose final WHERE expression
    mutation_filter = [mutation_filter]
    # Only combine the mutation_filter and pos_filter if the mutation_filter exists
    # AND if the pos_filter exists
    if coordinate_ranges:
        if mutation_filter[0].as_string(conn):
            mutation_filter.append(sql.SQL(" AND "))
        mutation_filter.append(pos_filter)

    mutation_filter = sql.Composed(mutation_filter)

    # Only add the WHERE clause if the mutation_filter exists
    if mutation_filter.join("").as_string(conn):
        mutation_filter = (sql.SQL("WHERE") + mutation_filter).join(" ")

    return mutation_filter, mutation_table


def build_sequence_where_filter(
    group_key,
    start_date=None,
    end_date=None,
    subm_start_date=None,
    subm_end_date=None,
    selected_metadata_fields=None,
    selected_group_fields=None,
    selected_reference=None,
    include_reference=True,
):
    """Build query for filtering sequences based on user's location/date
    selection and selected metadata fields

    Parameters
    ----------
    group_key: str
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
    selected_group_fields: dict
        - Strucutred as { group_key: [group_vals] }
        - Key are group types, i.e., "lineage"
        - Values are a list of group values, i.e., ["B.1.617.2", "BA.1"]
    selected_reference: str
        - Reference name (e.g., "NC_012920.1")
    include_reference: bool
        - Include reference filter (this flag only used to generate
          filter for the coverage query)

    Returns
    -------
    sequence_query: sql.SQL
        - Ready for format-injection into a main query, as a
          inline table expression or CTE

    """

    if not start_date:
        start_date = config["min_date"]
    if not end_date:
        end_date = datetime.date.today().isoformat()

    if not selected_metadata_fields:
        selected_metadata_fields = {}
    if not selected_group_fields:
        selected_group_fields = {}

    # Construct submission date filters
    if not subm_start_date and not subm_end_date:
        submission_date_filter = sql.SQL("")
    else:
        chunks = []
        if subm_start_date:
            chunks.append(
                sql.SQL('"submission_date" >= {}').format(
                    sql.Literal(pd.to_datetime(subm_start_date))
                )
            )
        if subm_end_date:
            chunks.append(
                sql.SQL('"submission_date" <= {}').format(
                    sql.Literal(pd.to_datetime(subm_end_date))
                )
            )

        submission_date_filter = sql.Composed(
            [sql.SQL(" AND "), sql.SQL(" AND ").join(chunks)]
        )

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

    # Only display mutations in the currently selected reference
    if group_key == constants["GROUP_MUTATION"] and include_reference:
        metadata_filters.append(
            sql.SQL('"reference" = {reference}').format(
                reference=sql.Literal(selected_reference)
            )
        )

    if metadata_filters:
        metadata_filters = sql.SQL(" AND ").join(metadata_filters)
        metadata_filters = sql.Composed([metadata_filters, sql.SQL(" AND ")])
    else:
        metadata_filters = sql.SQL("")

    group_filters = []
    for group_key, group_vals in selected_group_fields.items():
        # Skip if no group values provided
        if not group_vals:
            continue

        # Skip if group key is not valid
        if group_key not in config["group_cols"].keys():
            continue

        group_filters.append(
            sql.SQL("{field} IN {vals}").format(
                field=sql.Identifier(group_key),
                vals=sql.Literal(tuple([str(val) for val in group_vals])),
            )
        )

    if group_filters:
        group_filters = sql.SQL(" AND ").join(group_filters)
        group_filters = sql.Composed([group_filters, sql.SQL(" AND ")])
    else:
        group_filters = sql.SQL("")

    sequence_where_filter = sql.SQL(
        """
        {metadata_filters}
        {group_filters}
        "collection_date" >= {start_date} AND "collection_date" <= {end_date}
        {submission_date_filter}
        """
    ).format(
        metadata_filters=metadata_filters,
        group_filters=group_filters,
        start_date=sql.Literal(pd.to_datetime(start_date)),
        end_date=sql.Literal(pd.to_datetime(end_date)),
        submission_date_filter=submission_date_filter,
    )

    return sequence_where_filter


def get_loc_level_ids(req):
    """Get location level IDs from a flask request"""

    res = {}

    for loc_level in constants["GEO_LEVELS"].values():
        loc_ids = req.get(loc_level, None)
        if not loc_ids:
            res[loc_level] = []
        else:
            res[loc_level] = loc_ids

    return res


def build_sequence_location_where_filter(group_key, loc_level_ids, *args, **kwargs):
    """Build query for filtering sequences based on user's location/date
    selection and selected metadata fields - including location

    Parameters
    ----------
    group_key: str
    loc_level_ids: dict
        - key: One of GEO_LEVELS
        - value: list of level IDs
    """

    sequence_where_filter = [build_sequence_where_filter(group_key, *args, **kwargs)]
    loc_where = []
    for loc_level in constants["GEO_LEVELS"].values():
        loc_ids = loc_level_ids[loc_level]
        if not loc_ids:
            continue
        loc_where.append(
            sql.SQL("({loc_level_col} = ANY({loc_ids}))").format(
                loc_level_col=sql.Identifier(loc_level), loc_ids=sql.Literal(loc_ids),
            )
        )

    if loc_where:
        loc_where = sql.Composed(
            [sql.SQL("("), sql.SQL(" OR ").join(loc_where), sql.SQL(")")]
        )
        sequence_where_filter.append(sql.SQL(" AND "))
        sequence_where_filter.append(loc_where)

    sequence_where_filter = sql.Composed(sequence_where_filter)
    return sequence_where_filter


def count_coverage(
    cur,
    sequence_where_filter,
    selected_reference,
    dna_or_aa,
    coordinate_mode,
    coordinate_ranges,
    selected_gene,
    selected_protein,
):
    """Get coverage stats from partial sequences"""

    coverage_filter = [
        sql.SQL('"reference" = {reference}').format(
            reference=sql.Literal(selected_reference)
        )
    ]

    if dna_or_aa == constants["DNA_OR_AA"]["DNA"]:
        coverage_table = "dna_coverage"
        gene_or_protein_col = sql.SQL("")
        gene_or_protein_coalesce = sql.SQL("")
        gene_or_protein_join = sql.SQL("")
        gene_or_protein_partition_by = sql.SQL("")
        gene_or_protein_df_col = []
    else:
        gene_or_protein_col = sql.SQL('"feature",')
        gene_or_protein_coalesce = sql.SQL(
            'COALESCE(start_count."feature", end_count."feature") AS "feature",'
        )
        gene_or_protein_join = sql.SQL(
            'start_count."feature" = end_count."feature" AND '
        )
        gene_or_protein_partition_by = sql.SQL('PARTITION BY "feature"')
        gene_or_protein_df_col = ["feature"]

        if coordinate_mode == constants["COORDINATE_MODES"]["COORD_GENE"]:
            coverage_table = "gene_aa_coverage"
            coverage_filter.append(
                sql.SQL('"feature" = {feature}').format(
                    feature=sql.Literal(selected_gene)
                )
            )

        elif coordinate_mode == constants["COORDINATE_MODES"]["COORD_PROTEIN"]:
            coverage_table = "protein_aa_coverage"
            coverage_filter.append(
                sql.SQL('"feature" = {feature}').format(
                    feature=sql.Literal(selected_protein)
                )
            )

    coverage_filter = sql.SQL(" AND ").join(coverage_filter)

    coverage_query = sql.SQL(
        """
        WITH selected AS (
            SELECT "isolate_id"
            FROM "metadata"
            WHERE {sequence_where_filter}
        ), start_count AS (
            SELECT {gene_or_protein_col} "range_start", COUNT(*) AS "count"
            FROM {coverage_table}
            INNER JOIN SELECTED ON {coverage_table}."isolate_id" = selected."isolate_id"
            WHERE {coverage_filter}
            GROUP BY {gene_or_protein_col} "range_start"
        ), end_count AS (
            SELECT {gene_or_protein_col} "range_end", COUNT(*) AS "count"
            FROM {coverage_table}
            INNER JOIN SELECTED ON {coverage_table}."isolate_id" = selected."isolate_id"
            WHERE {coverage_filter}
            GROUP BY {gene_or_protein_col} "range_end"
        ), start_end_count AS (
            SELECT
                {gene_or_protein_coalesce}
                COALESCE(start_count."range_start", end_count."range_end") AS "ind",
                COALESCE(start_count."count", 0) AS "start",
                COALESCE(end_count."count", 0) AS "end"
            FROM start_count
            FULL JOIN end_count ON 
                {gene_or_protein_join}
                start_count."range_start" = end_count."range_end"
            ORDER BY {gene_or_protein_col} "ind" ASC
        )
        SELECT
            {gene_or_protein_col}
            "ind",
            (SUM("start") OVER ({gene_or_protein_partition_by} ORDER BY "ind" ASC) -
            SUM("end") OVER ({gene_or_protein_partition_by} ORDER BY "ind" ASC))::INTEGER AS "count"
        FROM start_end_count
        ORDER BY {gene_or_protein_col} "ind" ASC
        """
    ).format(
        sequence_where_filter=sequence_where_filter,
        gene_or_protein_col=gene_or_protein_col,
        coverage_table=sql.Identifier(coverage_table),
        coverage_filter=coverage_filter,
        gene_or_protein_coalesce=gene_or_protein_coalesce,
        gene_or_protein_join=gene_or_protein_join,
        gene_or_protein_partition_by=gene_or_protein_partition_by,
    )

    cur.execute(coverage_query)

    coverage_df = pd.DataFrame.from_records(
        cur.fetchall(), columns=gene_or_protein_df_col + ["ind", "count"],
    )

    # print(coverage_df)

    return coverage_df


def query_and_aggregate(conn, req):
    """Select sequences and aggregate results
    If in "mutation" mode, return sequences grouped by co-occurring mutations
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
        - In "mutation" mode, "group_id" represents co-occurring mutations,
          and is structured as a list of mutation IDs (integers)
        - In any other mode, "group_id" is a string of the group
          name. i.e., in lineage grouping, group_id is the lineage
          name.
    """

    group_key = req.get("group_key", None)
    dna_or_aa = req.get("dna_or_aa", None)
    coordinate_mode = req.get("coordinate_mode", None)
    coordinate_ranges = req.get("coordinate_ranges", None)
    selected_gene = req.get("selected_gene", None)
    selected_protein = req.get("selected_protein", None)
    selected_reference = req.get("selected_reference", None)

    with conn.cursor() as cur:

        main_query = []
        for loc_level in constants["GEO_LEVELS"].values():
            loc_ids = req.get(loc_level, None)
            if not loc_ids:
                continue

            sequence_where_filter = build_sequence_where_filter(
                group_key,
                req.get("start_date", None),
                req.get("end_date", None),
                req.get("subm_start_date", None),
                req.get("subm_end_date", None),
                req.get("selected_metadata_fields", None),
                req.get("selected_group_fields", None),
                req.get("selected_reference", None),
            )
            sequence_where_filter = sql.SQL(
                "{prior} AND {loc_level_col} = ANY({loc_ids})"
            ).format(
                prior=sequence_where_filter,
                loc_level_col=sql.Identifier(loc_level),
                loc_ids=sql.Literal(loc_ids),
            )

            if group_key == constants["GROUP_MUTATION"]:
                (mutation_filter, mutation_table) = build_coordinate_filters(
                    conn,
                    dna_or_aa,
                    coordinate_mode,
                    coordinate_ranges,
                    selected_gene,
                    selected_protein,
                )
                sequence_mutation_table = "sequence_" + mutation_table

                main_query.append(
                    sql.SQL(
                        """
                    SELECT
                        sm."location",
                        EXTRACT(EPOCH FROM sm."collection_date"),
                        sm."mutations",
                        COUNT(*) as "count"
                    FROM (
                        SELECT
                            sst."collection_date",
                            locdef."value" AS "location",
                            (sst.mutations & (
                                SELECT ARRAY_AGG("id")
                                FROM {mutation_table}
                                {mutation_filter}
                            )) as "mutations"
                        FROM {sequence_mutation_table} sst
                        INNER JOIN {loc_def_table} locdef ON sst.{loc_level_col} = locdef.id
                        WHERE {sequence_where_filter}
                    ) sm
                    GROUP BY sm.mutations, sm."location", sm.collection_date
                    """
                    ).format(
                        loc_level_col=sql.Identifier(loc_level),
                        loc_def_table=sql.Identifier("metadata_" + loc_level),
                        mutation_table=sql.Identifier(mutation_table),
                        mutation_filter=mutation_filter,
                        sequence_mutation_table=sql.Identifier(sequence_mutation_table),
                        sequence_where_filter=sequence_where_filter,
                    )
                )
            else:
                main_query.append(
                    sql.SQL(
                        """
                    SELECT
                        locdef."value" as "location",
                        EXTRACT(EPOCH FROM m."collection_date"),
                        m.{group_key},
                        COUNT(m."isolate_id") as "count"
                    FROM "metadata" m
                    INNER JOIN {loc_def_table} locdef ON m.{loc_level_col} = locdef."id"
                    WHERE {sequence_where_filter}
                    GROUP BY locdef."value", m."collection_date", m.{group_key}
                    """
                    ).format(
                        group_key=sql.Identifier(group_key),
                        loc_level_col=sql.Identifier(loc_level),
                        loc_def_table=sql.Identifier("metadata_" + loc_level),
                        sequence_where_filter=sequence_where_filter,
                    )
                )

        if not main_query:
            raise Exception("No locations provided")

        main_query = sql.SQL(" UNION ALL ").join(main_query)

        cur.execute(main_query)

        records_df = pd.DataFrame.from_records(
            cur.fetchall(),
            columns=["location", "collection_date", "group_id", "counts"],
        )

        res = {
            "records": records_df.to_dict(orient="records"),
        }

        # Count coverage for partial sequences
        # only do for mutation mode
        if group_key == constants["GROUP_MUTATION"]:
            sequence_location_where_filter = build_sequence_location_where_filter(
                group_key,
                get_loc_level_ids(req),
                req.get("start_date", None),
                req.get("end_date", None),
                req.get("subm_start_date", None),
                req.get("subm_end_date", None),
                req.get("selected_metadata_fields", None),
                req.get("selected_group_fields", None),
                req.get("selected_reference", None),
                include_reference=False,
            )
            coverage_df = count_coverage(
                cur,
                sequence_location_where_filter,
                selected_reference,
                dna_or_aa,
                coordinate_mode,
                coordinate_ranges,
                selected_gene,
                selected_protein,
            )
            res["coverage"] = coverage_df.to_dict(orient="records")

    return jsonify(res)
