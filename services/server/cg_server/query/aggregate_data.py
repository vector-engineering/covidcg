# coding: utf-8

"""Get data given filtering conditions

Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""

import json
import pandas as pd

from collections import defaultdict

from cg_server.config import config
from cg_server.constants import constants
from cg_server.query.selection import query_sequences


def query_aggregate_data(conn, req):

    res_df, res_snv = query_sequences(conn, req)

    group_key = req.get("group_key", None)

    group_col = group_key
    if group_key == constants["GROUP_SNV"]:
        group_col = "snp_id"

    # METADATA COUNTS
    # ---------------
    metadata_counts = dict()
    for key in config["metadata_cols"].keys():
        metadata_counts[key] = res_df[key].value_counts().to_dict()

    # Build a map of location_id -> location_group
    location_ids = req.get("location_ids", None)
    location_id_to_labels = defaultdict(list)
    for label, ids in location_ids.items():
        for i in ids:
            location_id_to_labels[i].append(label)

    res_df["location"] = res_df["location_id"].map(location_id_to_labels)
    res_snv["location"] = res_snv["location_id"].map(location_id_to_labels)

    # Versions for when we count each location separately
    res_df_explode_loc = res_df.explode("location")
    res_snv_explode_loc = res_snv.explode("location")

    # COLLAPSED DATA
    # -------------
    if (group_key == constants["GROUP_SNV"] and len(res_snv_explode_loc) == 0) or (
        group_key != constants["GROUP_SNV"] and len(res_df_explode_loc) == 0
    ):
        agg_sequences_location_group_date = pd.DataFrame(
            columns=["collection_date", "location", "group_id", "counts"]
        )
    elif group_key == constants["GROUP_SNV"]:
        agg_sequences_location_group_date = (
            res_snv_explode_loc.groupby(["Accession ID", "location"])
            .agg(
                group_id=("snp_id", tuple),
                collection_date=("collection_date", "first"),
            )
            .groupby(["collection_date", "location", "group_id"])
            .size()
            .rename("counts")
            .reset_index()
        )
    else:
        agg_sequences_location_group_date = (
            res_df_explode_loc.groupby(["collection_date", "location", group_col])
            .size()
            .rename("counts")
            .reset_index()
            .rename(columns={group_col: "group_id"})
        )

    res = """
    {{
        "aggLocationGroupDate": {agg_sequences_location_group_date},
        "metadataCounts": {metadata_counts}
    }}
    """.format(
        agg_sequences_location_group_date=agg_sequences_location_group_date.to_json(
            orient="records"
        ),
        metadata_counts=json.dumps(metadata_counts),
    )

    return res
