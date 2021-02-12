# coding: utf-8

import itertools
import json
import pandas as pd
import numpy as np

import urllib.request
import gzip

from collections import defaultdict, Counter, OrderedDict
from functools import partial
from flask import Flask, request, make_response
from pathlib import Path
from yaml import load, dump

from threading import Timer
from time import sleep

try:
    from yaml import CLoader as Loader, CDumper as Dumper
except ImportError:
    from yaml import Loader, Dumper

from .color import get_categorical_colormap, get_snv_colors, snv_colors

app = Flask(__name__)

# Load app configuration
with open("config/config_genbank.yaml", "r") as fp:
    config = load(fp.read(), Loader=Loader)

print(config)

with urllib.request.urlopen(config["data_package_url"]) as f:
    print("Download complete")
    f = gzip.decompress(f.read())
    print("Decompression complete")
    f = json.loads(f)
    print("Loaded into memory")

case_data = (
    pd.DataFrame.from_dict(f["case_data"], orient="columns")
    .set_index("Accession ID")
    .rename(
        columns={
            "dna_snp_str": "dna_snp",
            "gene_aa_snp_str": "gene_aa_snp",
            "protein_aa_snp_str": "protein_aa_snp",
        }
    )
)
case_data["collection_date"] = pd.to_datetime(case_data["collection_date"])
# print(case_data)

# Load metadata map
metadata_map = f["metadata_map"]
# print(list(metadata_map.keys()))

# Load global group counts
global_group_counts = f["global_group_counts"]
# Convert SNV keys from strings to integers
global_group_counts["dna_snp"] = {
    int(k): v for k, v in global_group_counts["dna_snp"].items()
}
global_group_counts["gene_aa_snp"] = {
    int(k): v for k, v in global_group_counts["gene_aa_snp"].items()
}
global_group_counts["protein_aa_snp"] = {
    int(k): v for k, v in global_group_counts["protein_aa_snp"].items()
}
# print(global_group_counts)

# Group consensus SNVs
group_consensus_snvs = f["group_consensus_snps"]
# print(group_consensus_snps)

# Build colormaps
sequence_groups = list(group_consensus_snvs.keys())
group_colormaps = dict()
for group in sequence_groups:
    group_colormaps[group] = get_categorical_colormap(
        list(group_consensus_snvs[group].keys())
    )


dna_snp = (
    pd.DataFrame.from_dict(metadata_map["dna_snp"], orient="index", columns=["id"])
    .reset_index()
    .rename(columns={"index": "snp_str"})
)
dna_snp = dna_snp.join(
    dna_snp["snp_str"]
    .str.split("|", expand=True)
    .rename(columns={0: "pos", 1: "ref", 2: "alt"})
).set_index("id")
dna_snp.loc[:, "pos"] = dna_snp["pos"].astype(int)
dna_snp["color"] = get_snv_colors(dna_snp.index.values)


def format_nt_snv(snv_str):
    # Don't do this if it's a special group
    #   if (Object.values(GROUPS).includes(snvStr)) {
    #     return snvStr;
    #   }

    # Print as REF POS ALT
    # i.e., 23403|A|G -> A23403G
    chunks = snv_str.split("|")
    return "{}{}{}".format(chunks[1], chunks[0], chunks[2])


dna_snp["snv_name"] = dna_snp["snp_str"].apply(format_nt_snv)

# print(dna_snp)

gene_aa_snp = (
    pd.DataFrame.from_dict(metadata_map["gene_aa_snp"], orient="index", columns=["id"])
    .reset_index()
    .rename(columns={"index": "snp_str"})
)
gene_aa_snp = gene_aa_snp.join(
    gene_aa_snp["snp_str"]
    .str.split("|", expand=True)
    .rename(columns={0: "gene", 1: "pos", 2: "ref", 3: "alt"})
).set_index("id")
gene_aa_snp.loc[:, "pos"] = gene_aa_snp["pos"].astype(int)
gene_aa_snp["color"] = get_snv_colors(gene_aa_snp.index.values)

protein_aa_snp = (
    pd.DataFrame.from_dict(
        metadata_map["protein_aa_snp"], orient="index", columns=["id"]
    )
    .reset_index()
    .rename(columns={"index": "snp_str"})
)
protein_aa_snp = protein_aa_snp.join(
    protein_aa_snp["snp_str"]
    .str.split("|", expand=True)
    .rename(columns={0: "protein", 1: "pos", 2: "ref", 3: "alt"})
).set_index("id")
protein_aa_snp.loc[:, "pos"] = protein_aa_snp["pos"].astype(int)
protein_aa_snp["color"] = get_snv_colors(protein_aa_snp.index.values)


def format_aa_snv(snv_str):
    # Don't do this if it's a special group
    #   if (Object.values(GROUPS).includes(snvStr)) {
    #     return snvStr;
    #   }

    # Print as GENE/PROTEIN · REF POS ALT
    # i.e., S|614|D|G -> S · D614G
    chunks = snv_str.split("|")
    return "{} · {}{}{}".format(chunks[0], chunks[2], chunks[1], chunks[3])


gene_aa_snp["snv_name"] = gene_aa_snp["snp_str"].apply(format_aa_snv)
protein_aa_snp["snv_name"] = protein_aa_snp["snp_str"].apply(format_aa_snv)

# print(gene_aa_snp)
# print(protein_aa_snp)

"""
class RepeatedTimer(object):
    def __init__(self, interval, function, *args, **kwargs):
        self._timer = None
        self.interval = interval
        self.function = function
        self.args = args
        self.kwargs = kwargs
        self.is_running = False
        self.start()

    def _run(self):
        self.is_running = False
        self.start()
        self.function(*self.args, **self.kwargs)

    def start(self):
        if not self.is_running:
            self._timer = Timer(self.interval, self._run)
            self._timer.start()
            self.is_running = True

    def stop(self):
        self._timer.cancel()
        self.is_running = False

msg = {"a": 0}

def hello(name):
    # print("Hello {}!".format(name))
    msg["a"] += 1
"""


def filter_coordinate_ranges(
    _df, dna_or_aa, coordinate_mode, coordinate_ranges, selected_gene, selected_protein
):
    # Filter for positions within the coordinate ranges
    include_mask = pd.Series(True, index=_df.index)
    # The range is inclusive on both sides, i.e., [start, end]
    for start, end in coordinate_ranges:
        include_mask = include_mask & (_df["pos"] >= start)
        include_mask = include_mask & (_df["pos"] <= end)
        if dna_or_aa == "AA" and coordinate_mode == "GENE":
            include_mask = include_mask & (_df["gene"] == selected_gene["gene"])
        elif dna_or_aa == "AA" and coordinate_mode == "PROTEIN":
            include_mask = include_mask & (
                _df["protein"] == selected_protein["protein"]
            )

    _df = _df.loc[include_mask, :]
    return _df


@app.route("/")
def hello_world():
    print("args", request.args)
    print("form", request.form)
    print("json", request.json)

    req = request.json
    if not req:
        return make_response(("No filter parameters given", 400))

    # REQUIRED?
    # lineage, clade, SNV
    group_key = req.get("group_key", "snv")
    dna_or_aa = req.get("dna_or_aa", "DNA")
    coordinate_mode = req.get("coordinate_mode", "GENE")

    snv_col = ""
    group_col = group_key
    if dna_or_aa == "DNA":
        snp_df = dna_snp
        snv_col = "dna_snp"
    else:
        if coordinate_mode == "GENE":
            snp_df = gene_aa_snp
            snv_col = "gene_aa_snp"
        else:
            snp_df = protein_aa_snp
            snv_col = "protein_aa_snp"

    if group_key == "snv":
        group_col = snv_col

    res_df = case_data

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
        res_df = res_df.loc[case_data["location_id"].isin(all_location_ids), :]

    # FILTER BY COORDINATE RANGE
    # --------------------------
    coordinate_ranges = req.get("coordinate_ranges", None)
    selected_gene = req.get("selected_gene", None)
    selected_protein = req.get("selected_protein", None)

    if coordinate_ranges is not None:

        snv_field = "dna_snp"
        snv_cols = ["pos"]
        if dna_or_aa == "AA":
            if coordinate_mode == "GENE":
                snv_field = "gene_aa_snp"
                snv_cols.append("gene")
            else:
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
    pre_metadata_counts = dict()
    for key in config["metadata_cols"].keys():
        pre_metadata_counts[key] = res_df[key].value_counts().to_dict()

    # FILTER BY METADATA FIELDS
    # -------------------------
    # Metadata filters will come in the form of a JSON object
    # of { metadata_field: metadata_value }
    metadata_filters = req.get("metadata_filters", None)
    if metadata_filters is not None:
        include_mask = pd.Series(True, index=res_df.index)
        for md_key, md_val in metadata_filters.items():
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
    # print(res_snv.head())

    # COUNT GROUPS
    # ------------
    group_counts = (
        (res_snv if group_key == "snv" else res_df).groupby(group_col).size().to_dict()
    )
    group_counts = OrderedDict(sorted(group_counts.items(), key=lambda kv: -1 * kv[1]))

    # FILTER OUT LOW COUNT GROUPS
    # ---------------------------
    low_count_filter = req.get("low_count_filter", "GROUP_COUNTS")
    max_group_counts = req.get("max_group_counts", 20)
    min_local_counts = req.get("min_local_counts", 10)
    min_global_counts = req.get("min_global_counts", 5)

    # Sort groups and take the top N
    if low_count_filter == "GROUP_COUNTS":
        # maxGroupCounts
        valid_groups = list(itertools.islice(group_counts, max_group_counts))
    elif low_count_filter == "LOCAL_COUNTS":
        # minLocalCountsToShow
        valid_groups = [k for k, v in group_counts.items() if v > min_local_counts]
    elif low_count_filter == "GLOBAL_COUNTS":
        # minGlobalCountsToShow
        _global_group_counts = global_group_counts[group_col]
        valid_groups = [
            k
            for k in group_counts.keys()
            if _global_group_counts[k] > min_global_counts
        ]
    else:
        return make_response(("Invalid low_count_filter", 400))

    # If we're in SNV mode, always show the reference
    # no matter what
    if group_key == "snv":
        valid_groups = set(valid_groups)
        valid_groups.add(-1)
        valid_groups = list(valid_groups)

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
    # print(counts_per_location)
    # print(counts_per_location_date)

    counts_per_location_date_group = (
        (res_snv if group_key == "snv" else res_df)
        .groupby(["location", "collection_date", group_col])
        .size()
        .rename("counts")
        .reset_index()
        .rename(columns={"collection_date": "date", group_col: "group_id"})
    )

    if group_key == "snv":
        counts_per_location_date_group = counts_per_location_date_group.join(
            snp_df[["snp_str", "snv_name", "color"]], on="group_id", how="left"
        ).rename(columns={"snp_str": "group", "snv_name": "group_name"})

        # Label reference group
        ref_inds = counts_per_location_date_group["group_id"] == -1
        counts_per_location_date_group.loc[
            ref_inds, ["group", "group_name"]
        ] = "Reference"
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
        ].map(group_colormaps[group_key])

    # Label invalid groups "Other"
    other_inds = (
        counts_per_location_date_group["group_id"]
        .map(dict(zip(valid_groups, np.repeat(False, len(valid_groups)))))
        .fillna(True)
    )
    counts_per_location_date_group.loc[other_inds, ["group", "group_name"]] = "Other"
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
    print(counts_per_date_group)

    # FILTER BY DATE RANGE
    # --------------------
    date_range = req.get("date_range", None)

    if date_range is not None:
        # In the form of [start, end]
        # where start/end should be in milliseconds from Unix epoch
        # i.e., from Date().getTime() in JavaScript
        start_date = pd.to_datetime(date_range[0], unit="ms")
        end_date = pd.to_datetime(date_range[1], unit="ms")

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

        print(counts_per_date_group)
        print(res_df.head())

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
    print(counts_per_group)

    # CHANGING POSITIONS
    # ------------------

    if group_key == 'snv':

        # Use SNV ids to get SNV data
        counts_per_group = (
            counts_per_group
            .join(snp_df[snv_cols + ['ref', 'alt']], on='group_id', how='left')
        )

        # Process "Reference" and "Other" group separately
        group_pos_inds = (counts_per_group['group'] != 'Other') & (counts_per_group['group'] != 'Reference')
        # Wipe out all SNV-specific data for these groups
        # Set position to -1 so they show up first
        counts_per_group.loc[~group_pos_inds, snv_cols + ['ref', 'alt']] = np.nan
        counts_per_group.loc[~group_pos_inds, 'pos'] = -1

        counts_per_group.sort_values('pos', inplace=True)
        counts_per_group.loc[:, 'pos'] = counts_per_group['pos'].astype(int)

        # Create position columns, fill with reference values
        for pos, ref in zip(
                counts_per_group.loc[group_pos_inds, 'pos'].values, 
                counts_per_group.loc[group_pos_inds, 'ref'].values
            ):
            counts_per_group[pos] = ref

        # Fill in alternative values
        for i, pos, alt in zip(
                counts_per_group.index[group_pos_inds].values, 
                counts_per_group.loc[group_pos_inds, 'pos'].values, 
                counts_per_group.loc[group_pos_inds, 'alt'].values
            ):
            counts_per_group.loc[i, pos] = alt
            
        counts_per_group.reset_index(drop=True, inplace=True)
        
    else:
        
        # Add the reference row
        counts_per_group = pd.concat([
            counts_per_group,
            pd.DataFrame(
                # group, group_id, group_name, counts, color, percent
                [('Reference', 'Reference', 'Reference', np.nan, clade_colors[0], np.nan)],
                columns=counts_per_group.columns
            )
        ], axis=0, ignore_index=True)

        # Get the consensus SNVs for each group
        counts_per_group[snv_col] = (
            counts_per_group['group']
            .map(group_consensus_snvs[group_key])
            .apply(lambda x: x[snv_col + '_ids'] if type(x) is dict else [])
        )
        
        # Get SNV data for all of the SNV ids
        group_snvs = (
            pd.DataFrame({
                snv_col: sorted(set(sum(counts_per_group[snv_col].values, [])))
            })
            .set_index(snv_col)
            .join(snp_df, how='left')
        )
        
        # Filter for only the SNVs in the specified coordinate ranges
        group_snvs = filter_coordinate_ranges(
            group_snvs,
            dna_or_aa, 
            coordinate_mode, 
            coordinate_ranges, 
            selected_gene, 
            selected_protein
        )
        group_snvs = group_snvs.sort_values('pos')
        
        # Fill in reference bases/residues
        for pos, ref in zip(group_snvs['pos'].values, group_snvs['ref'].values):
            counts_per_group[pos] = ref
            
        # Fill in alternative bases/residues
        for i, snvs in zip(counts_per_group.index.values, counts_per_group[snv_col].values):
            snvs = [snv for snv in snvs if snv in group_snvs.index]
            for snv in snvs:
                pos = group_snvs.at[snv, 'pos']
                alt = group_snvs.at[snv, 'alt']
                counts_per_group.loc[i, pos] = alt
        
        counts_per_group.drop(columns=[snv_col], inplace=True)

    print(counts_per_group)


    return res_df.head(10).to_json(orient="records")


# print("starting...")
# rt = RepeatedTimer(1, hello, "World")  # it auto-starts, no need of rt.start()
# try:
#     sleep(100)  # your long-running job goes here...
# finally:
#     rt.stop()  # better in a try/finally block to make sure the program ends!
