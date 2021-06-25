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
from cg_server.download import download_genomes, download_metadata, download_snvs
from cg_server.error import handle_db_errors
from cg_server.query import (
    query_aggregate_data,
    query_initial,
    query_metadata_fields,
    query_country_score,
    query_group_snv_frequencies,
)

from psycopg2 import pool
from cg_server.query.connection_pooling import get_conn_from_pool

app = Flask(__name__, static_url_path="", static_folder="dist")
Gzip(app)
CORS(app)
auth = HTTPBasicAuth()

# Load usernames/passwords via. an environment variable,
# as a comma-delimited string, "user1:pass1,user2:pass2"
users = {}
load_users = os.getenv("LOGINS", "")
load_users = [chunk for chunk in load_users.split(",") if chunk]
load_users = [(chunk.split(":")[0], chunk.split(":")[1]) for chunk in load_users]
for username, password in load_users:
    users[username] = generate_password_hash(password)


# Load allowed CORS domains
cors_domains = []
if os.getenv("FLASK_ENV", "development") == "development":
    # Allow any connections from localhost in development
    cors_domains.append(re.compile(r"http://localhost"))

# Whitelist the production hostname
cors_domains.append(config["prod_hostname"])


@auth.verify_password
def verify_password(username, password):
    if username in users and check_password_hash(users.get(username), password):
        return username


connection_options = {
    "dbname": os.environ["POSTGRES_DB"],
    "user": os.environ["POSTGRES_USER"],
    "password": os.environ["POSTGRES_PASSWORD"],
    "host": os.environ["POSTGRES_HOST"],
}
if port := os.getenv("POSTGRES_PORT", None):
    connection_options["port"] = port

conn_pool = psycopg2.pool.SimpleConnectionPool(1,
                                               os.environ["POSTGRES_MAX_CONN"],
                                               **connection_options)

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
                os.getenv("DATA_PATH", config["data_folder"]),
                filenames_as_dates=True,
            )

        conn.commit()
        conn_pool.putconn(conn)


@app.route("/")
@auth.login_required(optional=(not config["login_required"]))
def index():
    return app.send_static_file("index.html")


@app.route("/init")
@cross_origin(origins=cors_domains)
@handle_db_errors(conn=get_conn_from_pool(connection_options, conn_pool),
                  options=connection_options,
                  conn_pool=conn_pool)
def init():
    conn = get_conn_from_pool(connection_options, conn_pool)
    return query_initial(conn)


@app.route("/country_score", methods=["GET"])
@cross_origin(origins=cors_domains)
@handle_db_errors(conn=get_conn_from_pool(connection_options, conn_pool),
                  options=connection_options,
                  conn_pool=conn_pool)
def _country_score():
    conn = get_conn_from_pool(connection_options, conn_pool)
    return query_country_score(conn)


@app.route("/data", methods=["GET", "POST"])
@cross_origin(origins=cors_domains)
@handle_db_errors(conn=get_conn_from_pool(connection_options, conn_pool),
                  options=connection_options,
                  conn_pool=conn_pool)
def get_sequences():
    req = request.json
    if not req:
        return make_response(("No filter parameters given", 400))
    conn = get_conn_from_pool(connection_options, conn_pool)
    return query_aggregate_data(conn, req)


@app.route("/metadata_fields", methods=["GET", "POST"])
@cross_origin(origins=cors_domains)
@handle_db_errors(conn=get_conn_from_pool(connection_options, conn_pool),
                  options=connection_options,
                  conn_pool=conn_pool)
def _get_metadata_fields():
    req = request.json
    conn = get_conn_from_pool(connection_options, conn_pool)
    return query_metadata_fields(conn, req)


@app.route("/group_snv_frequencies", methods=["POST"])
@cross_origin(origins=cors_domains)
@handle_db_errors(conn=get_conn_from_pool(connection_options, conn_pool),
                  options=connection_options,
                  conn_pool=conn_pool)
def get_group_snv_frequencies():
    req = request.json
    conn = get_conn_from_pool(connection_options, conn_pool)
    return query_group_snv_frequencies(conn, req)


@app.route("/download_metadata", methods=["POST"])
@cross_origin(origins=cors_domains)
@handle_db_errors(conn=get_conn_from_pool(connection_options, conn_pool),
                  options=connection_options,
                  conn_pool=conn_pool)
def _download_metadata():
    if not config["allow_metadata_download"]:
        return make_response(
            ("Metadata downloads not permitted on this version of COVID CG", 403)
        )
    req = request.json
    conn = get_conn_from_pool(connection_options, conn_pool)
    return download_metadata(conn, req)


@app.route("/download_snvs", methods=["POST"])
@cross_origin(origins=cors_domains)
@handle_db_errors(conn=get_conn_from_pool(connection_options, conn_pool),
                  options=connection_options,
                  conn_pool=conn_pool)
def _download_snvs():
    if not config["allow_metadata_download"]:
        return make_response(
            ("Metadata downloads not permitted on this version of COVID CG", 403)
        )

    req = request.json
    conn = get_conn_from_pool(connection_options, conn_pool)
    return download_snvs(conn, req)


@app.route("/download_genomes", methods=["POST"])
@cross_origin(origins=cors_domains)
@handle_db_errors(conn=get_conn_from_pool(connection_options, conn_pool),
                  options=connection_options,
                  conn_pool=conn_pool)
def _download_genomes():
    if not config["allow_genome_download"]:
        return make_response(
            ("Metadata downloads not permitted on this version of COVID CG", 403)
        )

    req = request.json
    conn = get_conn_from_pool(connection_options, conn_pool)
    return download_genomes(conn, req)
