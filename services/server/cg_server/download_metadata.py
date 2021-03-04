# coding: utf-8

import pandas as pd
import psycopg2

from flask import make_response

from cg_server.config import config
from cg_server.constants import constants
from cg_server.query import select_sequences


def download_metadata(conn, req):

    with conn.cursor() as cur:
        temp_table_name = select_sequences(cur, req)

        sequence_cols = [
            "Accession ID",
            "collection_date",
            "submission_date",
        ]
        sequence_cols_expr = ['q."{}"'.format(col) for col in sequence_cols]
        metadata_joins = []

        # Location columns
        location_cols = ["region", "country", "division", "location"]
        sequence_cols.extend(location_cols)
        sequence_cols_expr.extend(['loc."{}"'.format(col) for col in location_cols])

        for grouping in config["group_cols"].keys():
            sequence_cols.append(grouping)
            sequence_cols_expr.append('q."{}"'.format(grouping))

        for field in config["metadata_cols"].keys():
            sequence_cols.append(field)
            sequence_cols_expr.append(
                """
                metadata_{field}."value" as "{field}"
                """.format(
                    field=field
                )
            )
            metadata_joins.append(
                """
                INNER JOIN "metadata_{field}" metadata_{field} 
                    ON q."database" = metadata_{field}."id"
                """.format(
                    field=field
                )
            )

        # CTEs and evaluating the metadata joins separately
        # from the SNV joins speeds this up by a lot
        query = """
            WITH dss AS (
                SELECT 
                    q."id" as "sequence_id", 
                    array_to_string(array_agg(ds."snp_str"), ';') as "snp"
                FROM "{temp_table_name}" q
                INNER JOIN "sequence_dna_snp" sds ON q."id" = sds."sequence_id"
                INNER JOIN "dna_snp" ds ON sds."snp_id" = ds."id"
                GROUP BY q."id"
            ),
            gass AS (
                SELECT 
                    q."id" as "sequence_id", 
                    array_to_string(array_agg(gas."snp_str"), ';') as "snp"
                FROM "{temp_table_name}" q
                INNER JOIN "sequence_gene_aa_snp" sgas ON q."id" = sgas."sequence_id"
                INNER JOIN "gene_aa_snp" gas ON sgas."snp_id" = gas."id"
                GROUP BY q."id"
            ),
            pass AS (
                SELECT 
                    q."id" as "sequence_id", 
                    array_to_string(array_agg(pas."snp_str"), ';') as "snp"
                FROM "{temp_table_name}" q
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
                FROM "{temp_table_name}" q
                INNER JOIN dss ON q."id" = dss."sequence_id"
                INNER JOIN gass ON q."id" = gass."sequence_id"
                INNER JOIN pass ON q."id" = pass."sequence_id"
            )
            SELECT 
                {sequence_cols_expr}, 
                qq."dna_snp",
                qq."gene_aa_snp",
                qq."protein_aa_snp"
            FROM "{temp_table_name}" q
            INNER JOIN "location" loc ON q."location_id" = loc."id"
            {metadata_joins}
            JOIN qq ON qq."id" = q."id"
            """.format(
            temp_table_name=temp_table_name,
            sequence_cols_expr=",".join(sequence_cols_expr),
            metadata_joins="\n".join(metadata_joins),
        )
        # print(query)
        cur.execute(query)

        res_df = pd.DataFrame.from_records(
            cur.fetchall(),
            columns=sequence_cols + ["dna_snp", "gene_aa_snp", "protein_aa_snp"],
        )

    return make_response(res_df.to_csv(index=False), 200, {"Content-Type": "text/csv"})
