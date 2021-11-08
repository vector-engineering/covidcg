# coding: utf-8

"""Main flask app

Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""

import os
import psycopg2
import re

from flask import Flask, request, make_response
from flask_cors import CORS, cross_origin
from flask_gzip import Gzip
from flask_httpauth import HTTPBasicAuth
from werkzeug.security import generate_password_hash, check_password_hash

from cg_server.config import config
from cg_server.db_seed import seed_database, insert_sequences
from cg_server.download import download_genomes, download_metadata
from cg_server.error import handle_db_errors
from cg_server.dev_only import dev_only
from cg_server.test_connpool import raiseDatabaseError
from cg_server.query import (
    query_initial,
    query_metadata,
    query_country_score,
    query_group_snv_frequencies,
    query_and_aggregate,
    generate_report,
)

from pathlib import Path
from psycopg2 import pool
from cg_server.query.connection_pooling import get_conn_from_pool

# root/services/server/cg_server/main.py
project_root = Path(__file__).parent.parent.parent.parent

app = Flask(__name__, static_url_path="", static_folder="dist")
Gzip(app)

# Load allowed CORS domains
cors_domains = ["https://covidcg.org", config["prod_hostname"]]
if os.getenv("FLASK_ENV", "development") == "development":
    # Allow any connections from localhost in development
    cors_domains.append("http://localhost:{}".format(os.getenv("FRONTEND_PORT")))

# CORS config
CORS(app, origins=cors_domains)

auth = HTTPBasicAuth()

# Load usernames/passwords via. an environment variable,
# as a comma-delimited string, "user1:pass1,user2:pass2"
users = {}
load_users = os.getenv("LOGINS", "")
load_users = [chunk for chunk in load_users.split(",") if chunk]
load_users = [(chunk.split(":")[0], chunk.split(":")[1]) for chunk in load_users]
for username, password in load_users:
    users[username] = generate_password_hash(password)


@auth.verify_password
def verify_password(username, password):
    if username in users and check_password_hash(users.get(username), password):
        return username


connection_options = {
    "dbname": config["postgres_db"],
    "user": os.environ["POSTGRES_USER"],
    "password": os.environ["POSTGRES_PASSWORD"],
    "host": os.environ["POSTGRES_HOST"],
}
if port := os.getenv("POSTGRES_PORT", None):
    connection_options["port"] = port

conn_pool = None
try:
    conn_pool = psycopg2.pool.SimpleConnectionPool(
        1, os.getenv("POSTGRES_MAX_CONN", 20), **connection_options
    )
except Exception:
    # If an error is thrown, check if the database exists
    connection = psycopg2.connect(user=connection_options["user"],
                                  password=connection_options["password"],
                                  host=connection_options["host"],
                                  port=connection_options["port"])
    cur = connection.cursor()
    cur.execute(
        psycopg2.sql.SQL(
            """
            SELECT EXISTS (
                SELECT datname FROM pg_catalog.pg_database
                WHERE datname = {dbname}
            );
            """
        ).format(dbname=psycopg2.sql.Literal(connection_options["dbname"]))
    )

    if not cur.fetchone()[0]:
        connection.commit()
        connection.close()
        # Open a new connection in autocommit to create db
        connection = psycopg2.connect(user=connection_options["user"],
                                      password=connection_options["password"],
                                      host=connection_options["host"],
                                      port=connection_options["port"])
        connection.autocommit = True
        # If the database does not exist, create it
        cur = connection.cursor()
        cur.execute(
            psycopg2.sql.SQL(
                """
                CREATE DATABASE {dbname};
                """
            ).format(dbname=psycopg2.sql.Identifier(connection_options["dbname"]))
            )

    connection.commit()
    connection.close()

    # Create conn_pool
    conn_pool = psycopg2.pool.SimpleConnectionPool(
        1, os.getenv("POSTGRES_MAX_CONN", 20), **connection_options
    )


# Quickly check if our database has been initialized yet
# If not, then let's seed it
# Only allow in development mode
if os.getenv("FLASK_ENV", "development") == "development":
    conn = get_conn_from_pool(connection_options, conn_pool)
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

        if not cur.fetchone()[0]:
            print("Seeding DB")
            seed_database(conn)
            insert_sequences(
                conn,
                os.getenv("DATA_PATH", project_root / config["example_data_folder"]),
                filenames_as_dates=True,
            )

        conn.commit()
        conn_pool.putconn(conn)


@app.route("/seed")
def _seed():
    if os.getenv("FLASK_ENV", "development") != "development":
        return "no"

    conn = get_conn_from_pool(connection_options, conn_pool)
    with conn.cursor() as cur:
        print("Seeding DB")
        seed_database(conn)
        insert_sequences(
            conn,
            os.getenv("DATA_PATH", project_root / config["data_folder"]),
            filenames_as_dates=True,
        )
        conn.commit()
        conn_pool.putconn(conn)

    return "Done!"


@app.route("/")
@auth.login_required(optional=(not config["login_required"]))
def index():
    return app.send_static_file("index.html")


@app.route("/init")
@handle_db_errors(options=connection_options, conn_pool=conn_pool)
def init(conn):
    return query_initial(conn)


@app.route("/country_score", methods=["GET"])
@cross_origin(origins=cors_domains)
@handle_db_errors(options=connection_options, conn_pool=conn_pool)
def _country_score(conn):
    return query_country_score(conn)


@app.route("/data", methods=["GET", "POST"])
@cross_origin(origins=cors_domains)
@handle_db_errors(options=connection_options, conn_pool=conn_pool)
def get_sequences(conn):
    req = request.json
    if not req:
        return make_response(("No filter parameters given", 400))
    # return query_and_aggregate(conn, req)
    return query_and_aggregate(conn, req)


@app.route("/metadata", methods=["GET", "POST"])
@cross_origin(origins=cors_domains)
@handle_db_errors(options=connection_options, conn_pool=conn_pool)
def _get_metadata(conn):
    req = request.json
    return query_metadata(conn, req)


@app.route("/group_snv_frequencies", methods=["POST"])
@cross_origin(origins=cors_domains)
@handle_db_errors(options=connection_options, conn_pool=conn_pool)
def get_group_snv_frequencies(conn):
    req = request.json
    return query_group_snv_frequencies(conn, req)


@app.route("/download_metadata", methods=["POST"])
@cross_origin(origins=cors_domains)
@handle_db_errors(options=connection_options, conn_pool=conn_pool)
def _download_metadata(conn):
    if not config["allow_metadata_download"]:
        return make_response(
            ("Metadata downloads not permitted on this version of COVID CG", 403)
        )
    req = request.json
    return download_metadata(conn, req)


@app.route("/download_genomes", methods=["POST"])
@cross_origin(origins=cors_domains)
@handle_db_errors(options=connection_options, conn_pool=conn_pool)
def _download_genomes(conn):
    if not config["allow_genome_download"]:
        return make_response(
            ("Metadata downloads not permitted on this version of COVID CG", 403)
        )

    req = request.json
    return download_genomes(conn, req)


@app.route("/az_report", methods=["GET", "POST"])
@cross_origin(origins=cors_domains)
@handle_db_errors(options=connection_options, conn_pool=conn_pool)
def _generate_report(conn):
    req = request.args
    return generate_report(conn, req)
