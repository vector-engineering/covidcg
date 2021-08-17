# coding: utf-8

"""Load app constants

Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""

import json
import os

from pathlib import Path

constants = {}

constants_path = os.getenv(
    "CONSTANTSFILE",
    # root/services/server/cg_server/constants.py --> need to go back 3 levels
    Path(__file__).parent.parent.parent.parent / "src" / "constants" / "defs.json",
)

# Load constant defs
with open(constants_path, "r") as fp:
    constants = json.loads(fp.read())

