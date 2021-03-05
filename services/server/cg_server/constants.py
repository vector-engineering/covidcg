# coding: utf-8

"""Load app constants

Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""

import json

constants = {}

# Load constant defs
with open("/opt/constants/defs.json", "r") as fp:
    constants = json.loads(fp.read())

