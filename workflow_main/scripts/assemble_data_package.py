# coding: utf-8

"""Combine all data into one JSON file

Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""

import gzip
import json
import datetime


def assemble_data_package(
    case_data,
    clade_snp,
    country_score,
    dna_snp_map,
    gene_aa_snp_map,
    geo_select_tree,
    global_group_counts,
    lineage_snp,
    location_map,
    metadata_map,
    protein_aa_snp_map,
    data_package_out,
):
    data_package = {"data_date": datetime.date.today().isoformat()}

    with open(case_data, "r") as fp:
        data_package["case_data"] = json.loads(fp.read())
    with open(clade_snp, "r") as fp:
        data_package["clade_snp"] = json.loads(fp.read())
    with open(country_score, "r") as fp:
        data_package["country_score"] = json.loads(fp.read())
    with open(dna_snp_map, "r") as fp:
        data_package["dna_snp_map"] = json.loads(fp.read())
    with open(gene_aa_snp_map, "r") as fp:
        data_package["gene_aa_snp_map"] = json.loads(fp.read())
    with open(geo_select_tree, "r") as fp:
        data_package["geo_select_tree"] = json.loads(fp.read())
    with open(global_group_counts, "r") as fp:
        data_package["global_group_counts"] = json.loads(fp.read())
    with open(lineage_snp, "r") as fp:
        data_package["lineage_snp"] = json.loads(fp.read())
    with open(location_map, "r") as fp:
        data_package["location_map"] = json.loads(fp.read())
    with open(metadata_map, "r") as fp:
        data_package["metadata_map"] = json.loads(fp.read())
    with open(protein_aa_snp_map, "r") as fp:
        data_package["protein_aa_snp_map"] = json.loads(fp.read())

    with gzip.open(data_package_out, "wt") as fp:
        fp.write(json.dumps(data_package))
