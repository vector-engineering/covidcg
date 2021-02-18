# coding: utf-8

import itertools
import json
import os
import pandas as pd
import numpy as np

from collections import defaultdict, Counter, OrderedDict
from functools import partial
from flask import Flask, request, make_response
from flask_cors import CORS
from pathlib import Path
from yaml import load, dump

from time import sleep

try:
    from yaml import CLoader as Loader, CDumper as Dumper
except ImportError:
    from yaml import Loader, Dumper

from flask_server.color import (
    get_categorical_colormap,
    get_snv_colors,
    snv_colors,
    clade_colors,
)
from flask_server.load_data import load_data
from flask_server.RepeatedTimer import RepeatedTimer


app = Flask(__name__, static_url_path="", static_folder="dist")
CORS(app)

# Load app configuration
with open("config/config_genbank.yaml", "r") as fp:
    config = load(fp.read(), Loader=Loader)

print(config)

# Load constant defs
with open("src/constants/defs.json", "r") as fp:
    constants = json.loads(fp.read())

data = load_data(config["data_package_url"])


def filter_coordinate_ranges(
    _df, dna_or_aa, coordinate_mode, coordinate_ranges, selected_gene, selected_protein
):
    # Filter for positions within the coordinate ranges
    include_mask = pd.Series(True, index=_df.index)
    # The range is inclusive on both sides, i.e., [start, end]
    for start, end in coordinate_ranges:
        if dna_or_aa == constants["DNA_OR_AA"]["DNA"]:
            include_mask = include_mask & (_df["pos"] >= start)
            include_mask = include_mask & (_df["pos"] <= end)
        elif dna_or_aa == constants["DNA_OR_AA"]["AA"]:
            include_mask = include_mask & (_df["nt_pos"] >= start)
            include_mask = include_mask & (_df["nt_pos"] <= end)
            if coordinate_mode == constants["COORDINATE_MODES"]["COORD_GENE"]:
                include_mask = include_mask & (_df["gene"] == selected_gene["name"])
            elif coordinate_mode == constants["COORDINATE_MODES"]["COORD_PROTEIN"]:
                include_mask = include_mask & (
                    _df["protein"] == selected_protein["name"]
                )

    _df = _df.loc[include_mask, :]
    return _df


@app.route("/")
def index():
    return app.send_static_file("index.html")


@app.route("/init")
def init():
    return {
        "data_date": data["data_date"],
        "num_sequences": data["num_sequences"],
        "country_score": data["country_score"],
        "geo_select_tree": data["geo_select_tree"],
        "metadata_map": data["metadata_map"],
    }


@app.route("/data", methods=["GET", "POST"])
def get_sequences():
    print("args", request.args)
    print("form", request.form)
    print("json", request.json)

    req = request.json
    if not req:
        return make_response(("No filter parameters given", 400))

    # REQUIRED?
    # lineage, clade, SNV
    group_key = req.get("group_key", constants["GROUP_SNV"])
    dna_or_aa = req.get("dna_or_aa", constants["DNA_OR_AA"]["AA"])
    coordinate_mode = req.get(
        "coordinate_mode", constants["COORDINATE_MODES"]["COORD_GENE"]
    )

    snv_col = ""
    group_col = group_key
    if dna_or_aa == constants["DNA_OR_AA"]["DNA"]:
        snp_df = data["dna_snp"]
        snv_col = "dna_snp"
    else:
        if coordinate_mode == constants["COORDINATE_MODES"]["COORD_GENE"]:
            snp_df = data["gene_aa_snp"]
            snv_col = "gene_aa_snp"
        elif coordinate_mode == constants["COORDINATE_MODES"]["COORD_PROTEIN"]:
            snp_df = data["protein_aa_snp"]
            snv_col = "protein_aa_snp"

    if group_key == constants["GROUP_SNV"]:
        group_col = snv_col

    res_df = data["case_data"]

    # FILTER BY LOCATION
    # ------------------
    # location_ids should be a dictionary, of:
    # { "location_group": [location_ids...], ... }
    # Each group of location IDs will be grouped together when
    # counting sequences per location
    location_ids = req.get("location_ids", None)
    if location_ids is not None:

        if type(location_ids) is not dict:
            return make_response(("Invalid format for location_ids", 400))

        all_location_ids = sum(location_ids.values(), [])
        res_df = res_df.loc[data["case_data"]["location_id"].isin(all_location_ids), :]

    # FILTER BY COORDINATE RANGE
    # --------------------------
    coordinate_ranges = req.get("coordinate_ranges", None)
    selected_gene = req.get("selected_gene", None)
    selected_protein = req.get("selected_protein", None)

    if coordinate_ranges is not None:

        snv_field = "dna_snp"
        snv_cols = ["pos"]
        if dna_or_aa == constants["DNA_OR_AA"]["AA"]:
            snv_cols.append("nt_pos")
            if coordinate_mode == constants["COORDINATE_MODES"]["COORD_GENE"]:
                snv_field = "gene_aa_snp"
                snv_cols.append("gene")
            elif coordinate_mode == constants["COORDINATE_MODES"]["COORD_PROTEIN"]:
                snv_field = "protein_aa_snp"
                snv_cols.append("protein")

        # Only keep SNVs within the given NT/AA coordinate ranges

        # Get all SNVs from the selected sequences, then
        # join onto reference dataframe to get SNV positions
        res_snv = (
            res_df[[snv_field]]
            .explode(snv_field)
            .join(snp_df[snv_cols], on=snv_field, how="inner")
        )
        res_snv = filter_coordinate_ranges(
            res_snv,
            dna_or_aa,
            coordinate_mode,
            coordinate_ranges,
            selected_gene,
            selected_protein,
        )

        # Collapse back into a dataframe by Accession IDs, with
        # SNVs as lists per Accession ID
        res_snv_list = res_snv.groupby("Accession ID")[[snv_field]].agg(list)
        res_df = res_df.drop(columns=[snv_field]).join(res_snv_list, how="left")
        # res_df[snv_field].fillna({i: [] for i in res_df.index}, inplace=True)
        res_df.loc[:, snv_field] = res_df[snv_field].fillna("").apply(list)

    # INITIAL METADATA COUNTS
    # -----------------------
    num_sequences_before_metadata_filtering = len(res_df)
    pre_metadata_counts = dict()
    for key in config["metadata_cols"].keys():
        pre_metadata_counts[key] = res_df[key].value_counts().to_dict()

    # FILTER BY METADATA FIELDS
    # -------------------------
    # Metadata filters will come in the form of a JSON object
    # of { metadata_field: metadata_value }
    selected_metadata_fields = req.get("selected_metadata_fields", None)
    if selected_metadata_fields is not None:
        include_mask = pd.Series(True, index=res_df.index)
        for md_key, md_vals in selected_metadata_fields.items():
            for md_val in md_vals:
                include_mask = include_mask & (res_df[md_key] == md_val)
        res_df = res_df.loc[include_mask, :]

    # FILTER BY PATIENT AGE RANGE (special case)
    # age_range = req.get("age_range", None)
    # if age_range is not None:

    # METADATA COUNTS AFTER FILTERING
    # -------------------------------
    post_metadata_counts = dict()
    for key in config["metadata_cols"].keys():
        post_metadata_counts[key] = res_df[key].value_counts().to_dict()

    res_snv = res_snv.join(res_df[["collection_date", "location_id"]], how="right")
    res_snv[snv_field].fillna(-1, inplace=True)

    print("After metadata filtering")
    print(res_df)
    print(res_snv)

    # COUNT GROUPS
    # ------------
    group_counts = (
        (res_snv if group_key == constants["GROUP_SNV"] else res_df)
        .groupby(group_col)
        .size()
        .to_dict()
    )
    group_counts = OrderedDict(sorted(group_counts.items(), key=lambda kv: -1 * kv[1]))

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
    else:
        return make_response(("Invalid low_count_filter", 400))

    # If we're in SNV mode, always show the reference
    # no matter what
    if group_key == constants["GROUP_SNV"]:
        valid_groups = set(valid_groups)
        valid_groups.add(-1)
        valid_groups = list(valid_groups)

    print("Valid groups")
    print(valid_groups)
    # print(res_df.head())

    # COUNTS PER LOCATION
    # -------------------

    # Build a map of location_id -> location_group
    location_id_to_label = {}
    for label, ids in location_ids.items():
        for i in ids:
            location_id_to_label[i] = label

    res_df["location"] = res_df["location_id"].map(location_id_to_label)
    res_snv["location"] = res_snv["location_id"].map(location_id_to_label)

    counts_per_location = res_df.groupby("location").size().rename("counts").to_dict()
    counts_per_location_date = (
        res_df.groupby(["location", "collection_date"])
        .size()
        .rename("counts")
        .reset_index()
        .rename(columns={"collection_date": "date"})
    )
    counts_per_location_date["cumulative_count"] = counts_per_location_date.groupby(
        ["location"]
    ).cumsum()

    # print(counts_per_location)
    print("Counts per location date")
    print(counts_per_location_date)

    counts_per_location_date_group = (
        (res_snv if group_key == constants["GROUP_SNV"] else res_df)
        .groupby(["location", "collection_date", group_col])
        .size()
        .rename("counts")
        .reset_index()
        .rename(columns={"collection_date": "date", group_col: "group_id"})
    )

    if group_key == constants["GROUP_SNV"]:
        counts_per_location_date_group = counts_per_location_date_group.join(
            snp_df[["snp_str", "snv_name", "color"]], on="group_id", how="left"
        ).rename(columns={"snp_str": "group", "snv_name": "group_name"})

        # Label reference group
        ref_inds = counts_per_location_date_group["group_id"] == -1
        counts_per_location_date_group.loc[
            ref_inds, ["group", "group_name"]
        ] = constants["GROUPS"]["REFERENCE_GROUP"]
        counts_per_location_date_group.loc[ref_inds, "color"] = snv_colors[0]

    else:
        counts_per_location_date_group["group"] = counts_per_location_date_group[
            "group_id"
        ]
        counts_per_location_date_group["group_name"] = counts_per_location_date_group[
            "group_id"
        ]
        counts_per_location_date_group["color"] = counts_per_location_date_group[
            "group"
        ].map(data["group_colormaps"][group_key])

    # Label invalid groups "Other"
    other_inds = (
        counts_per_location_date_group["group_id"]
        .map(dict(zip(valid_groups, np.repeat(False, len(valid_groups)))))
        .fillna(True)
    )
    counts_per_location_date_group.loc[other_inds, ["group", "group_name"]] = constants[
        "GROUPS"
    ]["OTHER_GROUP"]
    counts_per_location_date_group.loc[other_inds, ["color"]] = "#AAA"

    # Add location counts?
    counts_per_location_date_group["location_counts"] = counts_per_location_date_group[
        "location"
    ].map(counts_per_location)

    print(counts_per_location_date_group)

    # AGGREGATE BY GROUP AND DATE (COLLAPSE LOCATION)
    # -----------------------------------------------

    counts_per_date_group = (
        counts_per_location_date_group.groupby(["group", "date"])
        .agg(
            group_id=("group_id", "first"),
            counts=("counts", "sum"),
            group_name=("group_name", "first"),
            color=("color", "first"),
        )
        .reset_index()
    )
    print("Counts per date & group")
    print(counts_per_date_group)

    # FILTER BY DATE RANGE
    # --------------------
    start_date = req.get("start_date", None)
    end_date = req.get("end_date", None)

    if start_date is not None and end_date is not None:
        # In the form of [start, end]
        # where start/end should be in milliseconds from Unix epoch
        # i.e., from Date().getTime() in JavaScript
        start_date = pd.to_datetime(start_date)
        end_date = pd.to_datetime(end_date)

        res_df = res_df.loc[
            (res_df["collection_date"] >= start_date)
            & (res_df["collection_date"] <= end_date),
            :,
        ]
        counts_per_date_group = counts_per_date_group.loc[
            (counts_per_date_group["date"] >= start_date)
            & (counts_per_date_group["date"] <= end_date),
            :,
        ]

        print("Counts per date & group (after date filtering)")
        print(counts_per_date_group)

    # COUNT GROUPS
    # ------------
    group_counts_after_date_filter = (
        counts_per_date_group.groupby("group")
        .agg(
            counts=("counts", "sum"),
            color=("color", "first"),
            group_name=("group_name", "first"),
        )
        .reset_index()
    )
    print(group_counts_after_date_filter)

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

    print("Counts per group")
    print(counts_per_group)

    # CHANGING POSITIONS
    # ------------------

    if group_key == constants["GROUP_SNV"]:
        # Use SNV ids to get SNV data
        counts_per_group = counts_per_group.join(
            snp_df[snv_cols + ["ref", "alt"]], on="group_id", how="left"
        )
        counts_per_group.reset_index(drop=True, inplace=True)
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

        # Get the consensus SNVs for each group
        counts_per_group[snv_col] = (
            counts_per_group["group"]
            .map(data["group_consensus_snvs"][group_key])
            .apply(lambda x: x[snv_col + "_ids"] if type(x) is dict else [])
        )

        # Get SNV data for all of the SNV ids
        group_snvs = (
            pd.DataFrame(
                {snv_col: sorted(set(sum(counts_per_group[snv_col].values, [])))}
            )
            .set_index(snv_col)
            .join(snp_df, how="left")
        )

        # Filter for only the SNVs in the specified coordinate ranges
        group_snvs = filter_coordinate_ranges(
            group_snvs,
            dna_or_aa,
            coordinate_mode,
            coordinate_ranges,
            selected_gene,
            selected_protein,
        )
        group_snvs = group_snvs.sort_values("pos")

        # Fill in reference bases/residues
        for i, pos, ref in zip(
            group_snvs.index.values, group_snvs["pos"].values, group_snvs["ref"].values
        ):
            pos = int(pos)
            counts_per_group["pos_" + str(pos)] = ref

        # Fill in alternative bases/residues
        for i, snvs in zip(
            counts_per_group.index.values, counts_per_group[snv_col].values
        ):
            snvs = [snv for snv in snvs if snv in group_snvs.index]
            for snv in snvs:
                pos = group_snvs.at[snv, "pos"]
                alt = group_snvs.at[snv, "alt"]
                counts_per_group.loc[i, "pos_" + str(pos)] = alt

        counts_per_group.drop(columns=[snv_col], inplace=True)

    print("Counts per group (changing positions)")
    print(counts_per_group)

    # SELECTED SNVS
    # -------------

    # if group_key == 'snv':

    res = """
    {{
        "filteredCaseData": {filtered_case_data},
        "dataAggLocationGroupDate": {counts_per_location_date_group},
        "dataAggGroupDate": {counts_per_date_group},
        "metadataCounts": {pre_metadata_counts},
        "metadataCountsAfterFiltering": {post_metadata_counts},
        "numSequencesBeforeMetadataFiltering": {num_sequences_before_metadata_filtering},
        "countsPerLocation": {counts_per_location},
        "countsPerLocationDate": {counts_per_location_date},
        "validGroups": {valid_groups},
        "countsPerGroup": {group_counts},
        "dataAggGroup": {counts_per_group},
        "countsPerGroupDateFiltered": {group_counts_after_date_filter}
    }}
    """.format(
        filtered_case_data=res_df.drop(
            columns=list(config["metadata_cols"].keys())
        ).to_json(orient="records"),
        counts_per_location_date_group=counts_per_location_date_group.to_json(
            orient="records"
        ),
        counts_per_date_group=counts_per_date_group.to_json(orient="records"),
        pre_metadata_counts=json.dumps(pre_metadata_counts),
        post_metadata_counts=json.dumps(post_metadata_counts),
        num_sequences_before_metadata_filtering=num_sequences_before_metadata_filtering,
        counts_per_location=json.dumps(counts_per_location),
        counts_per_location_date=counts_per_location_date.to_json(orient="records"),
        valid_groups=json.dumps({k: 1 for k in valid_groups}),
        group_counts=json.dumps(group_counts),
        counts_per_group=counts_per_group.to_json(orient="records"),
        group_counts_after_date_filter=group_counts_after_date_filter.to_json(
            orient="values"
        ),
    )

    return res


def reload_data():
    print("Reloading data")
    data = load_data(config["data_package_url"])


if os.environ.get("FLASK_ENV", None) == "production":
    rt = RepeatedTimer(config["data_package_refresh_freq"], reload_data)
# try:
#    sleep(1)  # your long-running job goes here...
# finally:
# rt.stop()  # better in a try/finally block to make sure the program ends!

