# coding: utf-8

"""Configure gunicorn workers

Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""

import multiprocessing

workers = multiprocessing.cpu_count() * 2 + 1
