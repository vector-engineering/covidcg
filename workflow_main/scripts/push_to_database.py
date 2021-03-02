#!/usr/bin/env python3
# coding: utf-8

import os
import psycopg2
import sys

from pathlib import Path

cg_server_path = str(
    Path(__file__).absolute().parent.parent.parent / "services" / "server"
)
# print(cg_server_path)
sys.path.append(cg_server_path)

from cg_server.database import seed_database


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

    seed_database(conn)

    conn.close()


if __name__ == "__main__":
    main()
