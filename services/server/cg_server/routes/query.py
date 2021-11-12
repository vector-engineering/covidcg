# coding: utf-8

"""Query routes

Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""

from flask import request, make_response
from flask_cors import cross_origin

from cg_server.app import app, cors_domains
from cg_server.db import get_db_connection

from cg_server.query import (
    query_initial,
    query_country_score,
    query_and_aggregate,
    query_metadata,
    query_group_mutation_frequencies,
)


@app.route("/init")
@get_db_connection()
def init(conn):
    return query_initial(conn)


@app.route("/country_score", methods=["GET"])
@cross_origin(origins=cors_domains)
@get_db_connection()
def country_score_(conn):
    return query_country_score(conn)


@app.route("/data", methods=["GET", "POST"])
@cross_origin(origins=cors_domains)
@get_db_connection()
def get_sequences(conn):
    req = request.json
    if not req:
        return make_response(("No filter parameters given", 400))
    # return query_and_aggregate(conn, req)
    return query_and_aggregate(conn, req)


@app.route("/metadata", methods=["GET", "POST"])
@cross_origin(origins=cors_domains)
@get_db_connection()
def get_metadata_(conn):
    req = request.json
    return query_metadata(conn, req)


@app.route("/group_mutation_frequencies", methods=["POST"])
@cross_origin(origins=cors_domains)
@get_db_connection()
def get_group_mutation_frequencies(conn):
    req = request.json
    return query_group_mutation_frequencies(conn, req)
