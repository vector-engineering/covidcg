# coding: utf-8

"""Database connection

Authors:
    - Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
    - David Favela - Vector Engineering Team (dfavela@broadinstitute.org)
"""

import functools
import os
import psycopg2
import sys
import traceback

from flask import make_response
from psycopg2 import pool
from random import random

connection_options = {
    "dbname": os.environ["POSTGRES_DB"],
    "user": os.environ["POSTGRES_USER"],
    "password": os.environ["POSTGRES_PASSWORD"],
    "host": os.environ["POSTGRES_HOST"],
}
if port := os.getenv("POSTGRES_PORT", None):
    connection_options["port"] = port

conn_pool = pool.SimpleConnectionPool(
    1, os.getenv("POSTGRES_MAX_CONN", 20), **connection_options
)


def get_conn_from_pool(connection_options=connection_options, conn_pool=conn_pool):
    if conn_pool:
        return conn_pool.getconn()
    else:
        return psycopg2.connect(**connection_options)


# Print to stderr (for Google Cloud Run error tracking)
# https://stackoverflow.com/a/14981125/4343866
def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


def get_db_connection(options=None, raiseError=False):
    def decorator_db_error(func):
        @functools.wraps(func)
        def wrapper_db_error(*args, **kwargs):
            try:
                get_conn_opts = {}
                if options:
                    get_conn_opts["connection_options"] = options

                conn = get_conn_from_pool(**get_conn_opts)
                if raiseError:
                    if random() < 0.4:
                        print("raisedError")
                        raise psycopg2.Error
                res = func(conn, *args, **kwargs)
            # Catch any database/SQL errors
            except psycopg2.Error as e:
                conn.rollback()
                traceback.print_exc(file=sys.stderr)
                return make_response((str(e), 500))
            # Catch any other error
            # Could just do this instead of the specific DB one,
            # but keep it separate in case we ever need to do
            # something special with DB error handling
            except Exception as e:
                conn.rollback()
                traceback.print_exc(file=sys.stderr)
                return make_response((str(e), 500))

            conn.commit()
            conn_pool.putconn(conn)

            return res

        return wrapper_db_error

    return decorator_db_error
