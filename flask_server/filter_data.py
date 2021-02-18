# coding: utf-8

import itertools
import pandas as pd
import numpy as np

from collections import OrderedDict

from flask_server.config import config
from flask_server.constants import constants


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
                include_mask = include_mask & (_df["gene"] == selected_gene)
            elif coordinate_mode == constants["COORDINATE_MODES"]["COORD_PROTEIN"]:
                include_mask = include_mask & (_df["protein"] == selected_protein)

    _df = _df.loc[include_mask, :]
    return _df


def get_grouping_options(data, group_key, dna_or_aa, coordinate_mode):
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

    return group_col, snv_col, snp_df


def filter_data(data, req):

    # COLLECT OPTIONS
    # ---------------

    group_key = req.get("group_key", None)
    dna_or_aa = req.get("dna_or_aa", None)
    coordinate_mode = req.get("coordinate_mode", None)

    group_col, snv_col, snp_df = get_grouping_options(
        data, group_key, dna_or_aa, coordinate_mode
    )

    res_df = data["case_data"]

    # FILTER BY LOCATION
    # ------------------
    # location_ids should be a dictionary, of:
    # { "location_group": [location_ids...], ... }
    # Each group of location IDs will be grouped together when
    # counting sequences per location
    location_ids = req.get("location_ids", None)
    if location_ids is not None:
        all_location_ids = sum(location_ids.values(), [])
        res_df = res_df.loc[data["case_data"]["location_id"].isin(all_location_ids), :]

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

    # FILTER BY COORDINATE RANGE
    # --------------------------
    coordinate_ranges = req.get("coordinate_ranges", None)
    selected_gene = req.get("selected_gene", None)
    selected_protein = req.get("selected_protein", None)

    if coordinate_ranges is not None:

        snv_cols = ["pos"]
        if dna_or_aa == constants["DNA_OR_AA"]["AA"]:
            snv_cols.append("nt_pos")
            if coordinate_mode == constants["COORDINATE_MODES"]["COORD_GENE"]:
                snv_cols.append("gene")
            elif coordinate_mode == constants["COORDINATE_MODES"]["COORD_PROTEIN"]:
                snv_cols.append("protein")

        # Only keep SNVs within the given NT/AA coordinate ranges

        # Get all SNVs from the selected sequences, then
        # join onto reference dataframe to get SNV positions
        res_snv = (
            res_df[[snv_col]]
            .explode(snv_col)
            .join(snp_df[snv_cols], on=snv_col, how="inner")
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
        res_snv_list = res_snv.groupby("Accession ID")[[snv_col]].agg(list)
        res_df = res_df.drop(columns=[snv_col]).join(res_snv_list, how="left")
        # res_df[snv_col].fillna({i: [] for i in res_df.index}, inplace=True)
        res_df.loc[:, snv_col] = res_df[snv_col].fillna("").apply(list)

    # INITIAL METADATA COUNTS
    # -----------------------
    metadata_counts = dict()
    for key in config["metadata_cols"].keys():
        metadata_counts[key] = res_df[key].value_counts().to_dict()

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

    res_snv = res_snv.join(res_df[["collection_date", "location_id"]], how="right")
    res_snv[snv_col].fillna(-1, inplace=True)
    # print("After metadata filtering")
    # print(res_df)
    # print(res_snv)

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

    # print("Valid groups")
    # print(valid_groups)
    # print(res_df.head())

    return res_df, res_snv, metadata_counts, valid_groups
