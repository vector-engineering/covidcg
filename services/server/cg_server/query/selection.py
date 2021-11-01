# coding: utf-8

"""Get sequences/SNVs from database based on user selections

Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""

import pandas as pd
from psycopg2 import sql
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
        - Enum of constants["COORDINATE_MODES"] (NT, GENE AA, or PROTEIN AA)
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

    selected_metadata_fields = req.get("selected_metadata_fields", None)

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

    if metadata_filters:
        metadata_filters = sql.SQL(" AND ").join(metadata_filters)
        metadata_filters = sql.Composed([metadata_filters, sql.SQL(" AND ")])
    else:
        metadata_filters = sql.SQL("")

    sequence_where_filter = sql.SQL(
        """
        {metadata_filters}
        "collection_date" >= {start_date} AND "collection_date" <= {end_date}
        {submission_date_filter}
        """
    ).format(
        metadata_filters=metadata_filters,
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
                loc_level_col=sql.Identifier(loc_level),
                loc_ids=sql.Literal(loc_ids),
            )
        )
    loc_where = sql.Composed(
        [sql.SQL("("), sql.SQL(" OR ").join(loc_where), sql.SQL(")")]
    )
    sequence_where_filter.append(sql.SQL(" AND "))
    sequence_where_filter.append(loc_where)
    sequence_where_filter = sql.Composed(sequence_where_filter)
    return sequence_where_filter


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
                                FROM {snv_table}
                                {snv_filter}
                            )) as "mutations"
                        FROM {sequence_snv_table} sst 
                        INNER JOIN {loc_def_table} locdef ON sst.{loc_level_col} = locdef.id
                        WHERE {sequence_where_filter}
                    ) sm
                    GROUP BY sm.mutations, sm."location", sm.collection_date
                    """
                    ).format(
                        loc_level_col=sql.Identifier(loc_level),
                        loc_def_table=sql.Identifier("metadata_" + loc_level),
                        snv_table=sql.Identifier(snv_table),
                        snv_filter=snv_filter,
                        sequence_snv_table=sql.Identifier(sequence_snv_table),
                        sequence_where_filter=sequence_where_filter,
                    )
                )
            else:
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
                    GROUP BY locdef."value", m."collection_date", m."lineage"
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

        res = pd.DataFrame.from_records(
            cur.fetchall(),
            columns=["location", "collection_date", "group_id", "counts"],
        ).to_json(orient="records")

    return res
