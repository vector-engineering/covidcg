# coding: utf-8

"""Download raw genomes of selected sequences

Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""

import gzip
import psycopg2
import tempfile

from flask import make_response, send_file
from psycopg2 import sql

from cg_server.query import build_sequence_location_where_filter, get_loc_level_ids


def download_genomes(conn, req):
    if req["compress"]:
        fp = tempfile.NamedTemporaryFile(mode="w+b", delete=True, suffix=".fa.gz")
    else:
        fp = tempfile.NamedTemporaryFile(mode="w+b", delete=True, suffix=".fa")

    try:
        with conn.cursor() as cur:
            sequence_where_filter = build_sequence_location_where_filter(
                None,  # req.get("group_key", None),
                get_loc_level_ids(req),
                req.get("start_date", None),
                req.get("end_date", None),
                req.get("subm_start_date", None),
                req.get("subm_end_date", None),
                req.get("selected_metadata_fields", None),
                req.get("selected_group_fields", None),
                req.get("selected_reference", None),
                req.get("sequence_length", None),
                req.get("percent_ambiguous", None),
            )

            cur.execute(
                sql.SQL(
                    """
                    SELECT 
                        m."Accession ID",
                        s."sequence"
                    FROM (
                        SELECT UNNEST("accession_ids") AS "Accession ID"
                        FROM "metadata"
                        WHERE {sequence_where_filter}
                    ) m
                    INNER JOIN "sequence" s on m."Accession ID" = s."Accession ID";
                    """
                ).format(sequence_where_filter=sequence_where_filter)
            )

            if req["compress"]:
                fasta_file = gzip.open(fp, mode="w", compresslevel=6)
            else:
                fasta_file = fp

            counter = 0
            while seqs := cur.fetchmany(1000):
                counter += 1
                for seq in seqs:
                    fasta_file.write(
                        (">" + seq[0] + "\n" + seq[1] + "\n").encode("utf-8")
                    )

            if req["compress"]:
                fasta_file.close()

    except psycopg2.Error:
        # If something went wrong, then first cleanup the temp file
        fp.close()
        return make_response(("Error getting sequences", 500))

    fp.seek(0)
    # send_file() should close the file after it's sent,
    # which should mark it from deletion by the filesystem,
    # so we don't need to do any additional cleanup
    res = send_file(
        fp,
        mimetype=("application/gzip" if req["compress"] else "text/plain"),
        as_attachment=True,
        attachment_filename=("genomes.fa.gz" if req["compress"] else "genomes.fa"),
    )

    return res
