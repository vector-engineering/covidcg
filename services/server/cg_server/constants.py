# coding: utf-8

import json

constants = {}

# Load constant defs
with open("/opt/defs.json", "r") as fp:
    constants = json.loads(fp.read())

