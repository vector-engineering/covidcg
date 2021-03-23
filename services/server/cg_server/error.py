# coding: utf-8

"""Request error handling

Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""

import functools
import psycopg2

from flask import make_response


def handle_db_errors(conn):
    def decorator_db_error(func):
        @functools.wraps(func)
        def wrapper_db_error(*args, **kwargs):
            try:
                res = func(*args, **kwargs)
            except psycopg2.Error as e:
                conn.rollback()
                return make_response((str(e), 500))
            except Exception as e:
                conn.rollback()
                return make_response((str(e), 500))
            
            conn.commit()
            
            return res

        return wrapper_db_error

    return decorator_db_error
