# coding: utf-8

"""Get data given filtering conditions

Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""

import itertools
import json
import pandas as pd
import psycopg2
import numpy as np

from collections import OrderedDict, defaultdict

from cg_server.color import (
    snv_colors,
    clade_colors,
)
from cg_server.config import config
from cg_server.constants import constants
from cg_server.query.selection import query_sequences
from cg_server.query.consensus_snvs import query_consensus_snvs


def query_aggregate_data(conn, conn_pool, req):

    res_df, res_snv = query_sequences(conn, req)
    num_sequences = len(res_df)

    group_key = req.get("group_key", None)
    dna_or_aa = req.get("dna_or_aa", None)
    coordinate_mode = req.get("coordinate_mode", None)

    snv_col = ""
    group_col = group_key
    if dna_or_aa == constants["DNA_OR_AA"]["DNA"]:
        snv_col = "dna_snp"
    else:
        if coordinate_mode == constants["COORDINATE_MODES"]["COORD_GENE"]:
            snv_col = "gene_aa_snp"
        elif coordinate_mode == constants["COORDINATE_MODES"]["COORD_PROTEIN"]:
            snv_col = "protein_aa_snp"

    if group_key == constants["GROUP_SNV"]:
        group_col = "snp_id"

    # METADATA COUNTS
    # ---------------
    metadata_counts = dict()
    for key in config["metadata_cols"].keys():
        metadata_counts[key] = res_df[key].value_counts().to_dict()

    # COUNT GROUPS
    # ------------
    group_counts = (
        (res_snv if group_key == constants["GROUP_SNV"] else res_df)
        .groupby(group_col)
        .size()
        .to_dict()
    )
    group_counts = OrderedDict(sorted(group_counts.items(), key=lambda kv: -1 * kv[1]))

    # print('Group counts')
    # print(group_counts)

    # FILTER OUT LOW COUNT GROUPS
    # ---------------------------
    low_count_filter = req.get(
        "low_count_filter", constants["LOW_FREQ_FILTER_TYPES"]["GROUP_COUNTS"]
    )
    max_group_counts = req.get("max_group_counts", 20)
    min_local_counts = req.get("min_local_counts", 10)
    min_global_counts = req.get("min_global_counts", 5)

    # Sort groups and take the top N
    if low_count_filter == constants["LOW_FREQ_FILTER_TYPES"]["GROUP_COUNTS"]:
        # maxGroupCounts
        valid_groups = list(itertools.islice(group_counts, max_group_counts))
    elif low_count_filter == constants["LOW_FREQ_FILTER_TYPES"]["LOCAL_COUNTS"]:
        # minLocalCounts
        valid_groups = [k for k, v in group_counts.items() if v > min_local_counts]
    elif low_count_filter == constants["LOW_FREQ_FILTER_TYPES"]["GLOBAL_COUNTS"]:
        # minGlobalCounts
        _global_group_counts = data["global_group_counts"][group_col]
        valid_groups = [
            k
            for k in group_counts.keys()
            if _global_group_counts[k] > min_global_counts
        ]

    # If we're in SNV mode, always show the reference
    # no matter what
    if group_key == constants["GROUP_SNV"]:
        valid_groups = set(valid_groups)
        valid_groups.add(-1)
        valid_groups = list(valid_groups)

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

    res_df["location"] = res_df["location"].apply(lambda x: x[0])
    res_snv["location"] = res_snv["location"].apply(lambda x: x[0])

    # COLLAPSED DATA
    # -------------
    if (group_key == constants["GROUP_SNV"] and len(res_snv_explode_loc) == 0) or (
        group_key != constants["GROUP_SNV"] and len(res_df_explode_loc) == 0
    ):
        agg_sequences_location_group_date = pd.DataFrame(
            columns=["collection_date", "location", "group_id", "counts"]
        )
        agg_sequences_group_date = pd.DataFrame(
            columns=["collection_date", "group_id", "counts"]
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
        agg_sequences_group_date = (
            res_snv.groupby("Accession ID")
            .agg(
                group_id=("snp_id", tuple),
                collection_date=("collection_date", "first"),
            )
            .groupby(["collection_date", "group_id"])
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
        agg_sequences_group_date = (
            res_df_explode_loc.groupby(["collection_date", group_col])
            .size()
            .rename("counts")
            .reset_index()
            .rename(columns={group_col: "group_id"})
        )

    # COUNTS PER LOCATION
    # -------------------

    counts_per_location = (
        res_df_explode_loc.groupby("location").size().rename("counts").to_dict()
    )
    counts_per_location_date = (
        res_df_explode_loc.groupby(["location", "collection_date"])
        .size()
        .rename("counts")
        .reset_index()
        .rename(columns={"collection_date": "date"})
    )
    counts_per_location_date.insert(
        len(counts_per_location_date.columns),
        "cumulative_count",
        counts_per_location_date.groupby(["location"]).cumsum()["counts"],
    )

    # print("Counts per location")
    # print(counts_per_location)
    # print("Counts per location date")
    # print(counts_per_location_date)

    snv_aggs = {
        "group": ("snp_str", "first"),
        "group_name": ("snv_name", "first"),
        "counts": ("sequence_id", "count"),
        "color": ("color", "first"),
    }
    group_aggs = {
        "group": ("group_id", "first"),
        "group_name": ("group_id", "first"),
        "counts": ("Accession ID", "count"),
        "color": ("color", "first"),
    }
    if group_key == constants["GROUP_SNV"]:
        aggs = snv_aggs
    else:
        aggs = group_aggs

    counts_per_location_date_group = (
        res_snv_explode_loc
        if group_key == constants["GROUP_SNV"]
        else res_df_explode_loc
    ).rename(columns={"collection_date": "date", group_col: "group_id"})

    # Label invalid groups "Other"
    other_inds = (
        counts_per_location_date_group["group_id"]
        .map(dict(zip(valid_groups, np.repeat(False, len(valid_groups)))))
        .fillna(True)
    )
    counts_per_location_date_group.loc[other_inds, "group_id"] = constants["GROUPS"][
        "OTHER_GROUP"
    ]
    counts_per_location_date_group.loc[other_inds, ["color"]] = "#AAA"
    if group_key == constants["GROUP_SNV"]:
        counts_per_location_date_group.loc[
            other_inds, ["snp_str", "snv_name"]
        ] = constants["GROUPS"]["OTHER_GROUP"]

    counts_per_location_date_group = (
        counts_per_location_date_group.groupby(["location", "date", "group_id"])
        .agg(**aggs)
        .reset_index()
    )

    if group_key == constants["GROUP_SNV"]:
        # Label reference group
        ref_inds = counts_per_location_date_group["group_id"] == -1
        counts_per_location_date_group.loc[
            ref_inds, ["group", "group_name"]
        ] = constants["GROUPS"]["REFERENCE_GROUP"]
        counts_per_location_date_group.loc[ref_inds, "color"] = snv_colors[0]

    # Add location counts?
    counts_per_location_date_group["location_counts"] = counts_per_location_date_group[
        "location"
    ].map(counts_per_location)

    # print("Counts per location date group")
    # print(counts_per_location_date_group)

    # AGGREGATE BY GROUP AND DATE (COLLAPSE LOCATION)
    # -----------------------------------------------

    # Run the same operations for location_date_group, except this time
    # with the collapsed locations (one per sequence) instead of the
    # exploded ones (multiple locations per sequence possible)

    counts_per_date_group = (
        res_snv if group_key == constants["GROUP_SNV"] else res_df
    ).rename(columns={"collection_date": "date", group_col: "group_id"})

    # Label invalid groups "Other"
    other_inds = (
        counts_per_date_group["group_id"]
        .map(dict(zip(valid_groups, np.repeat(False, len(valid_groups)))))
        .fillna(True)
    )
    counts_per_date_group.loc[other_inds, "group_id"] = constants["GROUPS"][
        "OTHER_GROUP"
    ]
    counts_per_date_group.loc[other_inds, ["color"]] = "#AAA"
    if group_key == constants["GROUP_SNV"]:
        counts_per_date_group.loc[other_inds, ["snp_str", "snv_name"]] = constants[
            "GROUPS"
        ]["OTHER_GROUP"]

    counts_per_date_group = (
        counts_per_date_group.groupby(["date", "group_id"]).agg(**aggs).reset_index()
    )

    if group_key == constants["GROUP_SNV"]:
        # Label reference group
        ref_inds = counts_per_date_group["group_id"] == -1
        counts_per_date_group.loc[ref_inds, ["group", "group_name"]] = constants[
            "GROUPS"
        ]["REFERENCE_GROUP"]
        counts_per_date_group.loc[ref_inds, "color"] = snv_colors[0]

    # print("Counts per date & group")
    # print(counts_per_date_group)

    # COUNT GROUPS
    # ------------
    group_counts_after_low_freq = (
        counts_per_date_group.groupby("group")
        .agg(
            counts=("counts", "sum"),
            color=("color", "first"),
            group_name=("group_name", "first"),
        )
        .reset_index()
    )
    # print("Group counts after low freq collapse")
    # print(group_counts_after_low_freq)

    # AGGREGATE BY GROUP (COLLAPSE LOCATION + DATE)
    # ---------------------------------------------

    counts_per_group = (
        counts_per_date_group.groupby("group")
        .agg(
            group_id=("group_id", "first"),
            group_name=("group_name", "first"),
            counts=("counts", "sum"),
            color=("color", "first"),
        )
        .reset_index()
    )
    counts_per_group["percent"] = counts_per_group["counts"] / len(res_df)

    # print("Counts per group")
    # print(counts_per_group)

    # CHANGING POSITIONS
    # ------------------

    if group_key == constants["GROUP_SNV"]:

        def parse_snv_str(snv_str):
            chunks = snv_str.split("|")
            if dna_or_aa == constants["DNA_OR_AA"]["DNA"]:
                if len(chunks) < 3:
                    return (-1, None, None)
                return (int(chunks[0]), chunks[1], chunks[2])
            elif dna_or_aa == constants["DNA_OR_AA"]["AA"]:
                if len(chunks) < 4:
                    return (None, -1, None, None)
                return (chunks[0], int(chunks[1]), chunks[2], chunks[3])

        # Join SNV details (gene/protein, pos, ref, alt)
        if dna_or_aa == constants["DNA_OR_AA"]["DNA"]:
            snv_details = pd.DataFrame(
                [
                    [pos, ref, alt]
                    for pos, ref, alt in counts_per_group["group"]
                    .apply(parse_snv_str)
                    .values
                ],
                columns=["pos", "ref", "alt"],
                index=counts_per_group.index,
            )
        elif dna_or_aa == constants["DNA_OR_AA"]["AA"]:
            snv_details = pd.DataFrame(
                [
                    [gene_protein, pos, ref, alt]
                    for gene_protein, pos, ref, alt in counts_per_group["group"]
                    .apply(parse_snv_str)
                    .values
                ],
                columns=[
                    "gene"
                    if coordinate_mode == constants["COORDINATE_MODES"]["COORD_GENE"]
                    else "protein"
                ]
                + ["pos", "ref", "alt"],
                index=counts_per_group.index,
            )

        counts_per_group = pd.concat([counts_per_group, snv_details], axis=1)

    else:

        # Add the reference row
        counts_per_group = pd.concat(
            [
                counts_per_group,
                pd.DataFrame(
                    # group, group_id, group_name, counts, color, percent
                    [
                        (
                            constants["GROUPS"]["REFERENCE_GROUP"],
                            constants["GROUPS"]["REFERENCE_GROUP"],
                            constants["GROUPS"]["REFERENCE_GROUP"],
                            np.nan,
                            clade_colors[0],
                            np.nan,
                        )
                    ],
                    columns=counts_per_group.columns,
                ),
            ],
            axis=0,
            ignore_index=True,
        )

        group_snvs = query_consensus_snvs(conn, req, counts_per_group["group"].values)

        # Fill in reference bases/residues
        for pos, ref in zip(group_snvs["pos"].values, group_snvs["ref"].values):
            pos = int(pos)
            counts_per_group["pos_" + str(pos)] = ref

        # Fill in alternative bases/residues
        counts_per_group = counts_per_group.set_index("group")
        for i, row in group_snvs.iterrows():
            counts_per_group.loc[row["group"], "pos_" + str(row["pos"])] = row["alt"]
        counts_per_group = counts_per_group.reset_index()

    # print("Counts per group (changing positions)")
    # print(counts_per_group)

    res = """
    {{
        "aggSequencesLocationGroupDate": {agg_sequences_location_group_date},
        "aggSequencesGroupDate": {agg_sequences_group_date},
        "numSequences": {num_sequences},
        "dataAggLocationGroupDate": {counts_per_location_date_group},
        "dataAggGroupDate": {counts_per_date_group},
        "metadataCounts": {metadata_counts},
        "countsPerLocation": {counts_per_location},
        "countsPerLocationDate": {counts_per_location_date},
        "validGroups": {valid_groups},
        "dataAggGroup": {counts_per_group},
        "groupCounts": {group_counts_after_low_freq}
    }}
    """.format(
        agg_sequences_location_group_date=agg_sequences_location_group_date.to_json(
            orient="records"
        ),
        agg_sequences_group_date=agg_sequences_group_date.to_json(orient="records"),
        num_sequences=num_sequences,
        counts_per_location_date_group=counts_per_location_date_group.to_json(
            orient="records"
        ),
        counts_per_date_group=counts_per_date_group.to_json(orient="records"),
        metadata_counts=json.dumps(metadata_counts),
        counts_per_location=json.dumps(counts_per_location),
        counts_per_location_date=counts_per_location_date.to_json(orient="records"),
        valid_groups=json.dumps(valid_groups),
        counts_per_group=counts_per_group.to_json(orient="records"),
        group_counts_after_low_freq=group_counts_after_low_freq.to_json(
            orient="values"
        ),
    )

    conn_pool.putconn(conn)
    return res
