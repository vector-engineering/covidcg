# coding: utf-8

"""Main flask app

Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""

import os

from flask import Flask, request, make_response
from flask_cors import CORS, cross_origin
from flask_gzip import Gzip
from flask_httpauth import HTTPBasicAuth
from werkzeug.security import generate_password_hash, check_password_hash

from cg_server.config import config
from cg_server.seed import seed

from cassandra.cluster import Cluster
from pathlib import Path

# root/services/server/cg_server/main.py
# root/cassandra_test/cg_server/main.py
project_root = Path(__file__).parent.parent.parent

app = Flask(__name__, static_url_path="", static_folder="dist")
Gzip(app)

# Load allowed CORS domains
cors_domains = ['https://covidcg.org', config['prod_hostname']]
if os.getenv("FLASK_ENV", "development") == "development":
    # Allow any connections from localhost in development
    cors_domains.append("http://localhost:{}".format(os.getenv("FRONTEND_PORT")))

# CORS config
CORS(app, origins=cors_domains)

auth = HTTPBasicAuth()
# Load usernames/passwords via. an environment variable,
# as a comma-delimited string, "user1:pass1,user2:pass2"
users = {}
load_users = os.getenv("LOGINS", "")
load_users = [chunk for chunk in load_users.split(",") if chunk]
load_users = [(chunk.split(":")[0], chunk.split(":")[1]) for chunk in load_users]
for username, password in load_users:
    users[username] = generate_password_hash(password)

@auth.verify_password
def verify_password(username, password):
    if username in users and check_password_hash(users.get(username), password):
        return username

cluster = None
session = None

def connect():
    global cluster, session
    cluster = Cluster()
    session = cluster.connect()

connect()

@app.route("/")
@auth.login_required(optional=(not config["login_required"]))
def index():

    return app.send_static_file("index.html")

seed(session)