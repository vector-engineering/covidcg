# coding: utf-8

"""Create standalone global sequencing vega spec, with inlined data

Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""

import json


def standalone_map_spec(in_spec, data, out_spec):
    # Load vega spec JSON
    with open(in_spec, "r") as fp:
        spec = json.loads(fp.read())

    # Load data JSON
    with open(data, "r") as fp:
        data = json.loads(fp.read())

    # Inject data array into spec object
    spec["data"][0]["values"] = data

    # Save spec with injected data
    with open(out_spec, "w") as fp:
        fp.write(json.dumps(spec))
