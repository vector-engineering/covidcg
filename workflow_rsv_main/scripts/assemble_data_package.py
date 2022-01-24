# coding: utf-8

"""Combine all data into one JSON file

Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""

import gzip
import json
import datetime


def assemble_data_package(
    case_data,
    geo_select_tree,
<<<<<<< HEAD
    group_consensus_mutations,
=======
    group_consensus_snps,
>>>>>>> e6dd8312 (Rsvg workflow main (#420))
    metadata_map,
    data_package_out,
):
    data_package = {"data_date": datetime.date.today().isoformat()}

    with open(case_data, "r") as fp:
        data_package["case_data"] = json.loads(fp.read())
    with open(geo_select_tree, "r") as fp:
        data_package["geo_select_tree"] = json.loads(fp.read())
<<<<<<< HEAD
    with open(group_consensus_mutations, "r") as fp:
        data_package["group_consensus_mutations"] = json.loads(fp.read())
=======
    with open(group_consensus_snps, "r") as fp:
        data_package["group_consensus_snps"] = json.loads(fp.read())
>>>>>>> e6dd8312 (Rsvg workflow main (#420))
    with open(metadata_map, "r") as fp:
        data_package["metadata_map"] = json.loads(fp.read())

    with gzip.open(data_package_out, "wt") as fp:
        fp.write(json.dumps(data_package))
