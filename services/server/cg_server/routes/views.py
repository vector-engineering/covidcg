# coding: utf-8

"""View routes

Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""

from cg_server.app import app
from cg_server.auth import auth
from cg_server.config import config

@app.route("/")
@auth.login_required(optional=(not config["login_required"]))
def index():
    return app.send_static_file("index.html")
