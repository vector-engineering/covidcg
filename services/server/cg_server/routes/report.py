# coding: utf-8

"""AZ downloads

Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""

from flask import request

from cg_server.app import app, cors_domains
from cg_server.db import get_db_connection
from flask_cors import cross_origin

from cg_server.query import generate_report, build_variant_table


@app.route("/az_report", methods=["GET", "POST"])
@cross_origin(origins=cors_domains)
@get_db_connection()
def generate_report_(conn):
    req = request.args
    return generate_report(conn, req)


@app.route("/variant_table", methods=["GET", "POST"])
@cross_origin(origins=cors_domains)
@get_db_connection()
def build_variant_table_(conn):
    req = request.json
    return build_variant_table(conn, req)

