# coding: utf-8

"""Main flask app

Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""

import os

from cg_server.config import config, project_root
from cg_server.db_seed import seed_database, insert_sequences
from cg_server.db import get_conn_from_pool

from flask import Flask
from flask_gzip import Gzip
from flask_cors import CORS


app = Flask(__name__, static_url_path="", static_folder="dist")
Gzip(app)


cors_domains = ["https://covidcg.org", config["prod_hostname"]]
# Load allowed CORS domains
if os.getenv("FLASK_ENV", "development") == "development":
    # Allow any connections from localhost in development
    cors_domains.append("http://localhost:{}".format(os.getenv("FRONTEND_PORT")))
# Try to import the vm hostname if the file exists.
try:
    from cg_server.cors_conf import vmAddress
    cors_domains.append(vmAddress)
except Exception as e:
    pass
# CORS config
CORS(app, origins=cors_domains)

# Quickly check if our database has been initialized yet
# If not, then let's seed it
# Only allow in development mode
if os.getenv("FLASK_ENV", "development") == "development":
    conn = get_conn_from_pool()
    with conn.cursor() as cur:
        cur.execute(
            """
            SELECT EXISTS (
                SELECT FROM pg_tables
                WHERE  schemaname = 'public'
                AND    tablename  = 'metadata'
            );
            """
        )
        exists = cur.fetchone()[0]

    if not exists:
        print("Seeding DB")
        seed_database(conn)

        insert_sequences(
            conn,
            os.getenv("DATA_PATH", project_root / config["data_folder"]),
            filenames_as_dates=True,
        )
        conn.commit()
