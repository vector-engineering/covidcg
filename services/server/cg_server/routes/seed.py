# coding: utf-8

"""Database seeding (development only)

Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""

import os

from cg_server.app import app
from cg_server.config import config, project_root
from cg_server.db import get_db_connection
from cg_server.db_seed import seed_database, insert_sequences


@app.route("/force_seed")
@get_db_connection()
def force_seed(conn):
    if os.getenv("FLASK_ENV", "development") != "development":
        return "no"

    print("Seeding DB from /force_seed")
    seed_database(conn)
    insert_sequences(
        conn, os.getenv("DATA_PATH", project_root / config["example_data_folder"])
    )

    return "Done!"
