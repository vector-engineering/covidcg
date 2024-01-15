# coding: utf-8

"""Insert raw genomes into database

Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""

import io
import gzip

from psycopg2 import sql
from pathlib import Path


def flush_chunk(cur, buffer):
    # Drop old temp table
    cur.execute('DROP TABLE IF EXISTS "temp_sequence";')
    # Create new temp table
    cur.execute(
        """
        CREATE TEMP TABLE "temp_sequence" (
            "Accession ID"  TEXT       NOT NULL,
            sequence        TEXT       NOT NULL,
            filename        TEXT       NOT NULL
        );
        """
    )

    # Flush anything queued for writing, and reset pointer
    buffer.flush()
    buffer.seek(0)

    cur.copy_expert(
        """
        COPY "temp_sequence" FROM STDIN WITH (FORMAT CSV);
        """,
        buffer,
    )
    buffer.close()

    cur.execute(
        """
        INSERT INTO "sequence" ("Accession ID", "sequence", "filename")
            SELECT t."Accession ID", t."sequence", t."filename"
            FROM "temp_sequence" t
            INNER JOIN (
                SELECT UNNEST("accession_ids") AS "Accession ID"
                FROM "metadata"
            ) m ON t."Accession ID" = m."Accession ID"
        ON CONFLICT DO NOTHING;
        """
    )

    # Reset buffer
    buffer = io.StringIO()
    return buffer


def insert_sequences(conn, data_path, schema="public"):
    """Insert sequences into database

    Parameters
    ----------
    conn: psycopg2.Connection
    data_path: str
    schema: str

    Returns
    -------
    None
    """

    print("INSERTING SEQUENCES")

    # Get all fasta files
    fasta_path = Path(data_path) / "fasta_processed"

    with conn.cursor() as cur:
        cur.execute(sql.SQL("SET search_path TO {};").format(sql.Identifier(schema)))

        # The SARS-CoV-2 genomes are ~30kb, which seems like it
        # would slow down postgres, but as of PostgreSQL 12, large
        # text fields are stored in background tables so that the
        # rest of the text fields retain fast access. So I think
        # it's ok to just dump the entire genome into one cell
        cur.execute('DROP TABLE IF EXISTS "sequence";')
        cur.execute(
            """
            CREATE TABLE IF NOT EXISTS "sequence" (
                "Accession ID"  TEXT       NOT NULL,
                sequence        TEXT       NOT NULL,
                filename        TEXT       NOT NULL
            );
            """
        )

        # Get all filenames from the database
        # cur.execute(
        #     """
        #     SELECT "filename"
        #     FROM "sequence"
        #     GROUP BY "filename";
        #     """
        # )
        # loaded_files = [record[0] for record in cur.fetchall()]

        # Write everything into a string buffer
        buffer = io.StringIO()
        counter = 0
        chunk_counter = 0
        chunk_size = 10_000

        for i, fasta_file in enumerate(sorted(fasta_path.iterdir())):
            # if fasta_file.name in loaded_files:
            #     print("Skipping: ", fasta_file.name)
            #     continue
            # print("Inserting: ", fasta_file.name)

            with gzip.open(str(fasta_file), "rt") as fp:
                # TODO: rewrite so the entire file isn't in the buffer?
                #       I guess it doesn't really matter since the IO
                #       buffer is gonna all be in memory too, so there's
                #       no way for this to be super memory efficient without
                #       additional chunking
                lines = fp.readlines()
                cur_entry = ""
                cur_seq = ""

                for i, line in enumerate(lines):
                    # Strip whitespace
                    line = line.strip()

                    # Ignore empty lines that aren't the last line
                    if not line and i < (len(lines) - 1):
                        continue

                    # If we're on the last line, but the entry is empty, then skip
                    if not line and i == (len(lines) - 1) and not cur_seq:
                        continue

                    # If not the name of an entry, add this line to the current sequence
                    # (some FASTA files will have multiple lines per sequence)
                    if line[0] != ">":
                        cur_seq = cur_seq + line

                    # Start of another entry = end of the previous entry
                    if line[0] == ">" or i == (len(lines) - 1):
                        # Avoid capturing the first one and pushing an empty sequence
                        if cur_entry and cur_seq:
                            buffer.write(
                                cur_entry + "," + cur_seq + "," + fasta_file.name + "\n"
                            )

                        # Clear the entry and sequence
                        cur_entry = line[1:]
                        # Ignore anything past the first whitespace
                        if cur_entry:
                            cur_entry = cur_entry.split()[0]
                        cur_seq = ""

                    counter += 1

                    if counter >= chunk_size:
                        print(
                            "Writing chunk {} ({} sequences so far)".format(
                                chunk_counter, (chunk_counter + 1) * chunk_size
                            )
                        )
                        buffer = flush_chunk(cur, buffer)
                        counter = 0
                        chunk_counter += 1

        # Create indices
        cur.execute(
            """
            CREATE INDEX "ix_sequence_accession_id" ON "sequence"("Accession ID");
            """
        )
