# coding: utf-8

"""Wrapper for development only endpoints

Author: David Favela - Vector Engineering Team (dfavela@broadinstitute.org)
"""

from functools import wraps
from flask import make_response
import os


def dev_only(func):
    def decorator_dev_only(func):
        @wraps(func)
        def wrapper_dev_only(*args, **kwargs):
            if not os.getenv("FLASK_ENV", "development") == "development":
                return make_response("This is a development only endpoint",
                                     404)

            return func(*args, **kwargs)
        return wrapper_dev_only
    return decorator_dev_only
