# coding: utf-8

"""Combine all data into one JSON file

Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""

import gzip
import json
import datetime


def assemble_data_package(
    case_data,
    country_score,
    geo_select_tree,
    global_group_counts,
    group_consensus_snps,
    metadata_map,
    location_map,
    data_package_out,
):
    data_package = {"data_date": datetime.date.today().isoformat()}

    with open(case_data, "r") as fp:
        data_package["case_data"] = json.loads(fp.read())
    with open(country_score, "r") as fp:
        data_package["country_score"] = json.loads(fp.read())
    with open(geo_select_tree, "r") as fp:
        data_package["geo_select_tree"] = json.loads(fp.read())
    with open(global_group_counts, "r") as fp:
        data_package["global_group_counts"] = json.loads(fp.read())
    with open(group_consensus_snps, "r") as fp:
        data_package["group_consensus_snps"] = json.loads(fp.read())
    with open(metadata_map, "r") as fp:
        data_package["metadata_map"] = json.loads(fp.read())
    with open(location_map, "r") as fp:
        data_package["location_map"] = json.loads(fp.read())

    with gzip.open(data_package_out, "wt") as fp:
        fp.write(json.dumps(data_package))
