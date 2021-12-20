# coding: utf-8

"""Request authentication

Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""

import os

from flask_httpauth import HTTPBasicAuth
from werkzeug.security import generate_password_hash, check_password_hash

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
