# coding: utf-8

"""Download routes

Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""

from flask import request, make_response

from cg_server.app import app, cors_domains
from cg_server.config import config
from cg_server.db import get_db_connection
from flask_cors import cross_origin

from cg_server.download import download_metadata, download_genomes


@app.route("/download_metadata", methods=["POST"])
@cross_origin(origins=cors_domains)
@get_db_connection()
def download_metadata_(conn):
    if not config["allow_metadata_download"]:
        return make_response(
            ("Metadata downloads not permitted on this version of COVID CG", 403)
        )
    req = request.json
    return download_metadata(conn, req)


@app.route("/download_genomes", methods=["POST"])
@cross_origin(origins=cors_domains)
@get_db_connection()
def download_genomes_(conn):
    if not config["allow_genome_download"]:
        return make_response(
            ("Metadata downloads not permitted on this version of COVID CG", 403)
        )

    req = request.json
    return download_genomes(conn, req)

