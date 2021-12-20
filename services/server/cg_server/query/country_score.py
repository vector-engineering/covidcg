# coding: utf-8

"""Send country sequencing "leaderboard" data

Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""


def query_country_score(conn):
    with conn.cursor() as cur:
        cur.execute(
            """
            SELECT "value"
            FROM "jsons"
            WHERE "key" = 'country_score'
            """
        )
        return {"country_score": cur.fetchone()[0]}
