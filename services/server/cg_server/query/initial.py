# coding: utf-8

"""Send initial data package to front-end
Consists of SNV definitions, etc

Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""


def query_initial(conn):

    with conn.cursor() as cur:
        cur.execute(
            """
            SELECT "key", "value"
            FROM "jsons";
            """
        )
        res = cur.fetchall()

    return {i[0]: i[1] for i in res}
