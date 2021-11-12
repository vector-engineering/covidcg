# coding: utf-8

"""Main flask app

Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""

import os

from cg_server.config import config

from flask import Flask
from flask_gzip import Gzip
from flask_cors import CORS


app = Flask(__name__, static_url_path="", static_folder="dist")
Gzip(app)

# Load allowed CORS domains
cors_domains = ["https://covidcg.org", config["prod_hostname"]]
if os.getenv("FLASK_ENV", "development") == "development":
    # Allow any connections from localhost in development
    cors_domains.append("http://localhost:{}".format(os.getenv("FRONTEND_PORT")))

# CORS config
CORS(app, origins=cors_domains)

