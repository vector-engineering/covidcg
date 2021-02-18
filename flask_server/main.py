# coding: utf-8

import json
import os
import pandas as pd
import numpy as np

from functools import partial
from flask import Flask, request, make_response
from flask_cors import CORS
from flask_gzip import Gzip
from pathlib import Path

from time import sleep

from flask_server.color import (
    get_categorical_colormap,
    get_snv_colors,
    snv_colors,
    clade_colors,
)
from flask_server.config import config
from flask_server.constants import constants
from flask_server.load_data import load_data
from flask_server.filter_data import (
    filter_data,
    filter_coordinate_ranges,
    get_grouping_options,
)
from flask_server.RepeatedTimer import RepeatedTimer


app = Flask(__name__, static_url_path="", static_folder="dist")
Gzip(app)
CORS(app)

data = load_data(config["data_package_url"])


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


@app.route("/download_sequences", methods=["POST"])
def download_sequences():
    req = request.json
    res_df, _, _, _ = filter_data(data, req)

    # Join metadata definitions
    for key in config["metadata_cols"].keys():
        res_df[key] = res_df[key].map(data["metadata_map"][key])

    # Join locations
    res_df = res_df.join(data["location_map"], on="location_id")

    # Join SNV names
    snv_cols = ["dna_snp", "gene_aa_snp", "protein_aa_snp"]
    for snv_col in snv_cols:
        snp_df = data[snv_col]
        res_snv = (
            res_df[[snv_col]]
            .explode(snv_col)  # Expand to one SNV per row
            .join(snp_df[["snv_name"]], on=snv_col, how="inner")
            .reset_index()
            .groupby("Accession ID")["snv_name"]  # Collapse back by Accession ID
            .agg(lambda x: ";".join(list(x)))
            .rename(snv_col)
        )
        res_df = res_df.drop(columns=[snv_col]).join(res_snv)

    return make_response(res_df.to_csv(), 200, {"Content-Type": "text/csv"})


@app.route("/data", methods=["GET", "POST"])
def get_sequences():
    # print("args", request.args)
    # print("form", request.form)
    # print("json", request.json)

    req = request.json
    if not req:
        return make_response(("No filter parameters given", 400))

    res_df, res_snv, metadata_counts, valid_groups = filter_data(data, req)

    group_key = req.get("group_key", None)
    dna_or_aa = req.get("dna_or_aa", None)
    coordinate_mode = req.get("coordinate_mode", None)

    group_col, snv_col, snp_df = get_grouping_options(
        data, group_key, dna_or_aa, coordinate_mode
    )

    # COUNTS PER LOCATION
    # -------------------

    # Build a map of location_id -> location_group
    location_ids = req.get("location_ids", None)
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

    # print("Counts per location")
    # print(counts_per_location)
    # print("Counts per location date")
    # print(counts_per_location_date)

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

    # print("Counts per location date group")
    # print(counts_per_location_date_group)

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
    # print("Counts per date & group")
    # print(counts_per_date_group)

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
    # print('Group counts after date filter')
    # print(group_counts_after_date_filter)

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
        # Use SNV ids to get SNV data
        # counts_per_group = counts_per_group.join(
        #     snp_df[snv_cols + ["ref", "alt"]], on="group_id", how="left"
        # )
        # counts_per_group.reset_index(drop=True, inplace=True)
        pass
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

        coordinate_ranges = req.get("coordinate_ranges", None)
        selected_gene = req.get("selected_gene", None)
        selected_protein = req.get("selected_protein", None)

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

    # print("Counts per group (changing positions)")
    # print(counts_per_group)

    # COLLAPSE DATA
    # -------------

    if group_key == constants["GROUP_SNV"]:
        res_df[group_col] = res_df[group_col].apply(tuple)
    agg_sequences = (
        res_df.groupby(["collection_date", "location", group_col])
        .size()
        .rename("counts")
        .reset_index()
        .rename(columns={group_col: "group_id"})
    )

    res = """
    {{
        "aggSequences": {agg_sequences},
        "numSequences": {num_sequences},
        "dataAggLocationGroupDate": {counts_per_location_date_group},
        "dataAggGroupDate": {counts_per_date_group},
        "metadataCounts": {metadata_counts},
        "countsPerLocation": {counts_per_location},
        "countsPerLocationDate": {counts_per_location_date},
        "validGroups": {valid_groups},
        "dataAggGroup": {counts_per_group},
        "groupCounts": {group_counts_after_date_filter}
    }}
    """.format(
        agg_sequences=agg_sequences.to_json(orient="records"),
        num_sequences=agg_sequences["counts"].sum(),
        counts_per_location_date_group=counts_per_location_date_group.to_json(
            orient="records"
        ),
        counts_per_date_group=counts_per_date_group.to_json(orient="records"),
        metadata_counts=json.dumps(metadata_counts),
        counts_per_location=json.dumps(counts_per_location),
        counts_per_location_date=counts_per_location_date.to_json(orient="records"),
        valid_groups=json.dumps({k: 1 for k in valid_groups}),
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

