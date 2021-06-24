# coding: utf-8

"""Send country sequencing "leaderboard" data

Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""

import json

from flask import make_response


def query_country_score(conn, conn_pool):
    returnObj = {}
    with conn.cursor() as cur:
        cur.execute(
            """
            SELECT "value"
            FROM "country_score"
            """
        )
        returnObj = {"country_score": cur.fetchone()[0]}
    conn_pool.putconn(conn)
    return returnObj
