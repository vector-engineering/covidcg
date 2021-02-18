# coding: utf-8

import json

constants = {}

# Load constant defs
with open("src/constants/defs.json", "r") as fp:
    constants = json.loads(fp.read())

