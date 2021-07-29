# coding: utf-8

"""Request error handling

Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""

import functools
import psycopg2
import sys
import traceback

from flask import make_response
from random import random

from cg_server.query.connection_pooling import get_conn_from_pool


# Print to stderr (for Google Cloud Run error tracking)
# https://stackoverflow.com/a/14981125/4343866
def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


def handle_db_errors(options, conn_pool, raiseError=False):
    def decorator_db_error(func):
        @functools.wraps(func)
        def wrapper_db_error(*args, **kwargs):
            try:
                conn = get_conn_from_pool(options, conn_pool)
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
