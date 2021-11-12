#!/usr/bin/env python3
# coding: utf-8

import os
import psycopg2
import sys

from pathlib import Path

# workflow_main/scripts/push_to_database.py
cg_server_path = str(
    Path(__file__).resolve().parent.parent.parent / "services" / "server"
)
# print(cg_server_path)
sys.path.append(cg_server_path)

from cg_server.config import config
from cg_server.db_seed import insert_sequences, seed_database

# workflow_main/scripts/push_to_database.py
project_root = Path(__file__).resolve().parent.parent.parent
data_path = project_root / config["data_folder"]


def main():
    print("Pushing to server")

    connection_options = {
        "dbname": os.environ["POSTGRES_DB"],
        "user": os.environ["POSTGRES_USER"],
        "password": os.environ["POSTGRES_PASSWORD"],
        "host": os.environ["POSTGRES_HOST"],
    }
    if port := os.getenv("POSTGRES_PORT", None):
        connection_options["port"] = port

    conn = psycopg2.connect(**connection_options)

    try:
        with conn.cursor() as cur:
            cur.execute('DROP SCHEMA IF EXISTS "new" CASCADE;')
            cur.execute('CREATE SCHEMA "new";')

        seed_database(conn, schema="new")
        # insert_sequences(conn, data_path, schema="new", filenames_as_dates=True)

        print("Committing changes...", end="", flush=True)
        conn.commit()
        print("done")

        print("Switching schemas...", end="", flush=True)
        with conn.cursor() as cur:
            cur.execute('DROP SCHEMA IF EXISTS "old" CASCADE;')
            cur.execute('ALTER SCHEMA "public" RENAME TO "old";')
            cur.execute('ALTER SCHEMA "new" RENAME TO "public";')
        conn.commit()
        print("done")

    except psycopg2.Error as e:
        raise e

    conn.close()


if __name__ == "__main__":
    main()
