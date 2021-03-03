# coding: utf-8

import itertools
import json
import os
import pandas as pd
import psycopg2
import numpy as np

from collections import OrderedDict
from functools import partial
from flask import Flask, jsonify, request, make_response
from flask_cors import CORS
from flask_gzip import Gzip
from flask_httpauth import HTTPBasicAuth
from werkzeug.security import generate_password_hash, check_password_hash

from cg_server.color import (
    get_categorical_colormap,
    get_snv_colors,
    snv_colors,
    clade_colors,
)
from cg_server.config import config
from cg_server.constants import constants
from cg_server.database import seed_database
from cg_server.query import query_sequences, query_consensus_snvs, select_sequences
from cg_server.query_init import query_init, query_metadata_map

app = Flask(__name__, static_url_path="", static_folder="dist")
Gzip(app)
CORS(app)
auth = HTTPBasicAuth()

# Load usernames/passwords via. an environment variable,
# as a comma-delimited string, "user1:pass1,user2:pass2"
users = {}
load_users = os.getenv("LOGINS", "")
load_users = [
    (chunk.split(":")[0], chunk.split(":")[1]) for chunk in load_users.split(",")
]
for username, password in load_users:
    users[username] = generate_password_hash(password)


@auth.verify_password
def verify_password(username, password):
    if username in users and check_password_hash(users.get(username), password):
        return username


connection_options = {
    "dbname": os.environ["POSTGRES_DB"],
    "user": os.environ["POSTGRES_USER"],
    "password": os.environ["POSTGRES_PASSWORD"],
    "host": os.environ["POSTGRES_HOST"],
}
if port := os.getenv("POSTGRES_PORT", None):
    connection_options["port"] = port

conn = psycopg2.connect(**connection_options)


@app.route("/")
@auth.login_required(optional=(not config["login_required"]))
def index():
    return app.send_static_file("index.html")


@app.route("/seed")
def seed_db():
    # Only allow in development mode
    if os.getenv("FLASK_ENV", "development") == "development":
        print("Seeding DB")
        seed_database(conn)
        return "success"
    else:
        return make_response(("Forbidden in production", 400))


@app.route("/init")
def init():
    return query_init(conn)


@app.route("/data", methods=["GET", "POST"])
def get_sequences():
    req = request.json
    if not req:
        return make_response(("No filter parameters given", 400))

    res_df, res_snv = query_sequences(conn, req)

    # print(res_df)
    # print(res_snv)

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

    snv_aggs = {}
    if group_key == constants["GROUP_SNV"]:
        snv_aggs["counts"] = ("sequence_id", "count")
        snv_aggs["group"] = ("snp_str", "first")
        snv_aggs["group_name"] = ("snv_name", "first")
        snv_aggs["color"] = ("color", "first")
    else:
        snv_aggs["counts"] = ("Accession ID", "count")
        snv_aggs["color"] = ("color", "first")

    counts_per_location_date_group = (
        (res_snv if group_key == constants["GROUP_SNV"] else res_df)
        .groupby(["location", "collection_date", group_col])
        .agg(**snv_aggs)
        .reset_index()
        .rename(columns={"collection_date": "date", group_col: "group_id"})
    )

    if group_key == constants["GROUP_SNV"]:
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

    if group_key != constants["GROUP_SNV"]:

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

    # COLLAPSE DATA
    # -------------

    if group_key == constants["GROUP_SNV"]:
        agg_sequences = (
            res_snv.groupby("Accession ID")
            .agg(
                group_id=("snp_id", tuple),
                collection_date=("collection_date", "first"),
                location=("location", "first"),
            )
            .groupby(["collection_date", "location", "group_id"])
            .size()
            .rename("counts")
            .reset_index()
        )
    else:
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
        valid_groups=json.dumps(valid_groups),
        counts_per_group=counts_per_group.to_json(orient="records"),
        group_counts_after_date_filter=group_counts_after_date_filter.to_json(
            orient="values"
        ),
    )

    return res


@app.route("/download_metadata", methods=["POST"])
def download_metadata():
    req = request.json

    with conn.cursor() as cur:
        temp_table_name = select_sequences(cur, req)

        sequence_cols = [
            "Accession ID",
            "collection_date",
            "submission_date",
        ]
        sequence_cols_expr = ['q."{}"'.format(col) for col in sequence_cols]
        metadata_joins = []

        # Location columns
        location_cols = ["region", "country", "division", "location"]
        sequence_cols.extend(location_cols)
        sequence_cols_expr.extend(['loc."{}"'.format(col) for col in location_cols])

        for grouping in config["group_cols"].keys():
            sequence_cols.append(grouping)
            sequence_cols_expr.append('q."{}"'.format(grouping))

        for field in config["metadata_cols"].keys():
            sequence_cols.append(field)
            sequence_cols_expr.append(
                """
                metadata_{field}."value" as "{field}"
                """.format(
                    field=field
                )
            )
            metadata_joins.append(
                """
                JOIN "metadata_{field}" metadata_{field} 
                    ON q."database" = metadata_{field}."id"
                """.format(
                    field=field
                )
            )

        cur.execute(
            """
            SELECT {sequence_cols_expr}
            FROM "{temp_table_name}" q
            JOIN "location" loc ON q."location_id" = loc."id"
            {metadata_joins}
            """.format(
                temp_table_name=temp_table_name,
                sequence_cols_expr=",".join(sequence_cols_expr),
                metadata_joins="\n".join(metadata_joins),
            )
        )

        res_df = pd.DataFrame.from_records(cur.fetchall(), columns=sequence_cols)

    return make_response(res_df.to_csv(index=False), 200, {"Content-Type": "text/csv"})


@app.route("/download_snvs", methods=["POST"])
def download_snvs():
    req = request.json
    res_df, res_snv = query_sequences(conn, req)
    return make_response(
        res_snv.drop(
            columns=[
                "sequence_id",
                "snp_id",
                "location_id",
                "collection_date",
                "snp_str",
                "color",
            ]
        ).to_csv(index=False),
        200,
        {"Content-Type": "text/csv"},
    )

