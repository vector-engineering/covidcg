# coding: utf-8

"""Drops connection to the database

Author: David Favela - Vector Engineering Team (dfavela@broadinstitute.org)
"""

from functools import wraps


def raiseDatabaseError(kwargs):
    def decorator_test_connpool(func):
        @wraps(func)
        def wrapper_test_connpool(*args, **kwargs):
            return func(*args, **kwargs)
        return wrapper_test_connpool
    return decorator_test_connpool
