# coding: utf-8

"""Request error handling

Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""

import functools
import psycopg2
import sys

from flask import make_response

# Print to stderr (for Google Cloud Run error tracking)
# https://stackoverflow.com/a/14981125/4343866
def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


def handle_db_errors(conn):
    def decorator_db_error(func):
        @functools.wraps(func)
        def wrapper_db_error(*args, **kwargs):
            try:
                res = func(*args, **kwargs)
            # Catch any database/SQL errors
            except psycopg2.Error as e:
                conn.rollback()
                eprint(str(e))
                return make_response((str(e), 500))
            # Catch any other error
            # Could just do this instead of the specific DB one,
            # but keep it separate in case we ever need to do
            # something special with DB error handling
            except Exception as e:
                conn.rollback()
                eprint(str(e))
                return make_response((str(e), 500))

            conn.commit()

            return res

        return wrapper_db_error

    return decorator_db_error
