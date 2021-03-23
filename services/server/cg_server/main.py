# coding: utf-8

"""Main flask app

Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""

import json
import os
import psycopg2
import re

from flask import Flask, request, make_response
from flask_cors import CORS, cross_origin
from flask_gzip import Gzip
from flask_httpauth import HTTPBasicAuth
from werkzeug.security import generate_password_hash, check_password_hash

from cg_server.config import config
from cg_server.database import seed_database
from cg_server.download_metadata import download_metadata
from cg_server.download_genomes import download_genomes
from cg_server.download_snvs import download_snvs
from cg_server.get_metadata_fields import get_metadata_fields
from cg_server.insert_sequences import insert_sequences
from cg_server.query_init import query_init
from cg_server.query_data import query_data

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

conn = psycopg2.connect(**connection_options)

# Quickly check if our database has been initialized yet
# If not, then let's seed it
# Only allow in development mode
if os.getenv("FLASK_ENV", "development") == "development":
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


@app.route("/")
@auth.login_required(optional=(not config["login_required"]))
def index():
    return app.send_static_file("index.html")


@app.route("/init")
@cross_origin(origins=cors_domains)
def init():
    try:
        init = query_init(conn)
    except psycopg2.Error as e:
        conn.rollback()
        return make_response((str(e), 500))
    
    return init


@app.route("/data", methods=["GET", "POST"])
@cross_origin(origins=cors_domains)
def get_sequences():
    req = request.json
    if not req:
        return make_response(("No filter parameters given", 400))

    try:
        res = query_data(conn, req)
    except psycopg2.Error as e:
        conn.rollback()
        return make_response((str(e), 500))
    
    return res


@app.route("/metadata_fields", methods=["GET", "POST"])
@cross_origin(origins=cors_domains)
def _get_metadata_fields():
    req = request.json
    try:
        metadata_fields = get_metadata_fields(conn, req)
    except psycopg2.Error as e:
        conn.rollback()
        return make_response((str(e), 500))

    return metadata_fields


@app.route("/download_metadata", methods=["POST"])
@cross_origin(origins=cors_domains)
def _download_metadata():
    if not config["allow_metadata_download"]:
        return make_response(
            ("Metadata downloads not permitted on this version of COVID CG", 403)
        )

    req = request.json
    return download_metadata(conn, req)


@app.route("/download_snvs", methods=["POST"])
@cross_origin(origins=cors_domains)
def _download_snvs():
    if not config["allow_metadata_download"]:
        return make_response(
            ("Metadata downloads not permitted on this version of COVID CG", 403)
        )

    req = request.json
    return download_snvs(conn, req)


@app.route("/download_genomes", methods=["POST"])
@cross_origin(origins=cors_domains)
def _download_genomes():
    if not config["allow_genome_download"]:
        return make_response(
            ("Metadata downloads not permitted on this version of COVID CG", 403)
        )

    req = request.json
    return download_genomes(conn, req)

