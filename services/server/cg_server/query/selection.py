# coding: utf-8

"""Get sequences/mutations from database based on user selections

Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""

import pandas as pd
from psycopg2 import sql
from cg_server.config import config
from cg_server.constants import constants
from cg_server.config import config


def build_coordinate_filters(
    conn, dna_or_aa, coordinate_mode, coordinate_ranges, selected_gene, selected_protein,
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
    coordinate_ranges: list of [int, int]
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
                sql.SQL('"gene" = {gene}').format(gene=sql.Literal(selected_gene))
            )
        elif coordinate_mode == constants["COORDINATE_MODES"]["COORD_PROTEIN"]:
            mutation_table = "protein_aa_mutation"
            mutation_filter.append(
                sql.SQL('"protein" = {protein}').format(
                    protein=sql.Literal(selected_protein)
                )
            )

    mutation_filter = sql.SQL(" AND ").join(mutation_filter)

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
    mutation_filter = [mutation_filter]
    # Only combine the mutation_filter and pos_filter if the mutation_filter exists
    if mutation_filter[0].as_string(conn):
        mutation_filter.append(sql.SQL(" AND "))

    mutation_filter.append(pos_filter)
    mutation_filter = sql.Composed(mutation_filter)

    # Only add the WHERE clause if the mutation_filter exists
    if mutation_filter.join("").as_string(conn):
        mutation_filter = (sql.SQL("WHERE") + mutation_filter).join(" ")

    return mutation_filter, mutation_table


def build_sequence_where_filter(req):
    """Build query for filtering sequences based on user's location/date
    selection and selected metadata fields

    Parameters
    ----------
    req: flask.Request

    Request fields
    --------------
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

    Returns
    -------
    sequence_query: sql.SQL
        - Ready for format-injection into a main query, as a
          inline table expression or CTE

    """

    start_date = pd.to_datetime(req.get("start_date", None))
    end_date = pd.to_datetime(req.get("end_date", None))

    subm_start_date = req.get("subm_start_date", "")
    subm_end_date = req.get("subm_end_date", "")
    subm_start_date = None if subm_start_date == "" else pd.to_datetime(subm_start_date)
    subm_end_date = None if subm_end_date == "" else pd.to_datetime(subm_end_date)

    selected_reference = req.get("selected_reference", None)
    group_key = req.get("group_key", None)

    selected_metadata_fields = req.get("selected_metadata_fields", None)
    selected_group_fields = req.get("selected_group_fields", {})

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

    if selected_reference:
        if config['virus'] == "rsv" and group_key != "genotype":
            metadata_filters.append(
                sql.SQL('"genotype" = {genotype}').format(
                    genotype=sql.Literal(selected_reference)
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
        start_date=sql.Literal(start_date),
        end_date=sql.Literal(end_date),
        submission_date_filter=submission_date_filter,
    )

    return sequence_where_filter


def build_sequence_location_where_filter(req):
    sequence_where_filter = [build_sequence_where_filter(req)]
    loc_where = []
    for loc_level in constants["GEO_LEVELS"].values():
        loc_ids = req.get(loc_level, None)
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

    with conn.cursor() as cur:

        main_query = []
        for loc_level in constants["GEO_LEVELS"].values():
            loc_ids = req.get(loc_level, None)
            if not loc_ids:
                continue

            sequence_where_filter = build_sequence_where_filter(req)
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
                        sm."collection_date",
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
                group_col = ""
                if config["virus"] == "sars2":
                    group_col = "lineage"
                elif config["virus"] == "rsv":
                    group_col = "genotype"

                main_query.append(
                    sql.SQL(
                        """
                    SELECT
                        locdef."value" as "location",
                        m."collection_date",
                        m.{group_key},
                        COUNT(m."sequence_id") as "count"
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
                        group_col=sql.Identifier(group_col)
                    )
                )

        if not main_query:
            raise Exception("No locations provided")

        main_query = sql.SQL(" UNION ALL ").join(main_query)

        cur.execute(main_query)

        res = pd.DataFrame.from_records(
            cur.fetchall(),
            columns=["location", "collection_date", "group_id", "counts"],
        ).to_json(orient="records")

    return res
