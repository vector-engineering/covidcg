# coding: utf-8

"""Flask app entrypoint

Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""

from cg_server.app import app

# Trigger authentication
from cg_server.auth import *

from cg_server.routes import *

