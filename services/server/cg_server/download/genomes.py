# coding: utf-8

"""Download raw genomes of selected sequences

Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""

import gzip
import pandas as pd
import psycopg2
import tempfile

from flask import make_response, send_file
from psycopg2 import sql

# from cg_server.query.selection import select_sequences


def download_genomes(conn, req):

    fp = tempfile.NamedTemporaryFile(mode="w+b", delete=True, suffix=".fa.gz")

    try:
        with conn.cursor() as cur:
            # temp_table_name = select_sequences(conn, cur, req)

            cur.execute(
                sql.SQL(
                    """
                    SELECT q."Accession ID", s."sequence"
                    FROM {temp_table_name} q
                    JOIN "sequence" s on q."id" = s."sequence_id";
                    """
                ).format(temp_table_name=sql.Identifier(temp_table_name))
            )

            fasta_file = gzip.open(fp, mode="wt")
            counter = 0
            while seqs := cur.fetchmany(1000):
                counter += 1
                # print(counter)
                for seq in seqs:
                    fasta_file.write(">" + seq[0] + "\n")
                    fasta_file.write(seq[1])
                    fasta_file.write("\n")
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
        mimetype="application/gzip",
        as_attachment=True,
        attachment_filename="genomes.fa.gz",
    )

    return res
