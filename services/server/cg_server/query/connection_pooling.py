# coding: utf-8

"""Create ConnectionPool and get connection from pool

Author: David Favela - Vector Engineering Team (dfavela@broadinstitute.org)
"""

import psycopg2


def get_conn_from_pool(connection_options, conn_pool):
    if conn_pool:
        return conn_pool.getconn()
    else:
        return psycopg2.connect(**connection_options)
