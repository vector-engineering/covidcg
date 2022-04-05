#!/usr/bin/env python3
# coding: utf-8

"""Combine all data into one JSON file

Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""

import argparse
import gzip
import json
import datetime


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--case-data", type=str, required=True, help="Case data JSON file"
    )
    parser.add_argument(
        "--country-score", type=str, required=True, help="Country score JSON file"
    )
    parser.add_argument(
        "--geo-select-tree", type=str, required=True, help="Geo select tree JSON file"
    )
    parser.add_argument(
        "--global-group-counts",
        type=str,
        required=True,
        help="Global group counts JSON file",
    )
    parser.add_argument(
        "--group-consensus-mutations",
        type=str,
        required=True,
        help="Group consensus mutations JSON file",
    )
    parser.add_argument(
        "--metadata-map", type=str, required=True, help="Metadata map JSON file"
    )
    parser.add_argument(
        "--data-package-out", type=str, required=True, help="Data package output"
    )

    args = parser.parse_args()

    data_package = {"data_date": datetime.date.today().isoformat()}

    with open(args.case_data, "r") as fp:
        data_package["case_data"] = json.loads(fp.read())
    with open(args.country_score, "r") as fp:
        data_package["country_score"] = json.loads(fp.read())
    with open(args.geo_select_tree, "r") as fp:
        data_package["geo_select_tree"] = json.loads(fp.read())
    with open(args.global_group_counts, "r") as fp:
        data_package["global_group_counts"] = json.loads(fp.read())
    with open(args.group_consensus_mutations, "r") as fp:
        data_package["group_consensus_mutations"] = json.loads(fp.read())
    with open(args.metadata_map, "r") as fp:
        data_package["metadata_map"] = json.loads(fp.read())

    with gzip.open(args.data_package_out, "wt") as fp:
        fp.write(json.dumps(data_package))


if __name__ == "__main__":
    main()
