#!/usr/bin/env python3
# coding: utf-8

"""Global surveillance data for the home page

Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""

import argparse
import datetime
import json
import pandas as pd
import numpy as np

from scipy.stats import linregress
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--case-data", type=str, required=True, help="Path to case data JSON file",
    )

    parser.add_argument(
        "--location-map",
        type=str,
        required=True,
        help="Path to location map JSON file",
    )

    parser.add_argument(
        "--metadata-map", type=str, required=True, help="Path to metadata map JSON file"
    )

    parser.add_argument(
        "--start-date-days-ago",
        type=int,
        default=90,
        help="Number of days before today to filter data. Default: 90",
    )
    parser.add_argument(
        "--start-date",
        type=str,
        default=None,
        help="Start date for filtering data in ISO format (YYYY-MM-DD). Overrides --start-date-days-ago if defined. Default: None",
    )
    parser.add_argument(
        "--end-date-days-ago",
        type=int,
        default=30,
        help="Number of days before today to cut off data prior to regressions. Default: 0",
    )

    parser.add_argument(
        "--min-combo-count",
        type=int,
        default=50,
        help="Minimum counts for a spike SNV combo to be included in the dataset",
    )
    parser.add_argument(
        "--min-single-count",
        type=int,
        default=50,
        help="Minimum counts for a single spike SNV to be included in the dataset",
    )

    parser.add_argument(
        "-o", "--output", type=str, required=True, help="Path to output directory",
    )

    args = parser.parse_args()

    case_data = pd.read_json(args.case_data)
    location_map = pd.read_json(args.location_map)
    with open(args.metadata_map, "r") as fp:
        metadata_map = json.loads(fp.read())

    out_path = Path(args.output)

    # Get Spike SNVs
    gene_aa_snp = pd.DataFrame.from_dict(
        metadata_map["gene_aa_snp"], orient="index"
    ).reset_index()
    gene_aa_snp.columns = ["snv", "snv_id"]
    spike_snvs = gene_aa_snp.loc[gene_aa_snp["snv"].str.match(r"^S")]
    spike_snv_map = spike_snvs["snv"].to_dict()

    df = case_data[
        ["Accession ID", "collection_date", "lineage", "gene_aa_snp_str", "location_id"]
    ].join(location_map, on="location_id")

    # Filter for only SNVs in spike
    valid_snv_ids = spike_snvs["snv_id"].values
    df["spike_aa_snv"] = df["gene_aa_snp_str"].apply(
        lambda x: tuple([snv for snv in x if snv in valid_snv_ids])
    )

    # Filter for valid regions
    valid_regions = [
        "Africa",
        "Asia",
        "Europe",
        "North America",
        "Oceania",
        "South America",
    ]
    df = df.loc[df["region"].isin(valid_regions)]

    df["collection_date"] = pd.to_datetime(df["collection_date"])
    df["collection_week"] = df["collection_date"].dt.to_period("W")

    location_counts = (
        df.groupby(["region", "collection_week"])
        .size()
        .rename("location_counts")
        .reset_index()
    )

    if args.start_date:
        start_date_iso = args.start_date
    else:
        start_date_iso = (
            datetime.date.today() - datetime.timedelta(days=args.start_date_days_ago)
        ).isoformat()

    # End date only used for the regression
    end_date_iso = (
        datetime.date.today() - datetime.timedelta(days=args.end_date_days_ago)
    ).isoformat()

    # LINEAGE DATA
    lineage_counts = (
        df.loc[df["collection_date"] >= pd.to_datetime(start_date_iso)]
        .groupby(["region", "lineage", "collection_week"])
        .size()
        .reset_index()
        .rename(columns={"lineage": "group", 0: "counts"})
    )

    # SPIKE SNV COMBO DATA

    spike_combo_snv_counts = (
        df.loc[df["collection_date"] >= pd.to_datetime(start_date_iso)]
        .groupby(["region", "spike_aa_snv", "collection_week"])
        .size()
        .reset_index()
        .rename(columns={0: "counts", "spike_aa_snv": "group"})
    )
    spike_combo_snv_freq = (
        df.groupby("spike_aa_snv").size().sort_values(ascending=False)
    )
    spike_combo_snv_counts = spike_combo_snv_counts.loc[
        spike_combo_snv_counts["group"].isin(
            spike_combo_snv_freq.index[
                spike_combo_snv_freq >= args.min_combo_count
            ].values
        )
    ]

    # SPIKE SINGLE SNV DATA
    spike_single_snv_counts = (
        df.loc[df["collection_date"] >= pd.to_datetime(start_date_iso)]
        .groupby(["region", "spike_aa_snv", "collection_week"])
        .size()
        .reset_index()
        .rename(columns={0: "counts", "spike_aa_snv": "group"})
    )
    spike_single_snv_counts.loc[:, "group"] = spike_single_snv_counts["group"].apply(
        list
    )
    spike_single_snv_counts = (
        spike_single_snv_counts.explode("group")
        .assign(group=lambda x: x["group"].fillna(-1))
        .groupby(["region", "collection_week", "group"])
        .agg(counts=("counts", np.sum))
        .reset_index()
    )
    spike_single_snv_freq = (
        spike_single_snv_counts.groupby("group")["counts"]
        .sum()
        .sort_values(ascending=False)
    )
    spike_single_snv_counts = spike_single_snv_counts.loc[
        spike_single_snv_counts["group"].isin(
            spike_single_snv_freq.index[spike_single_snv_freq >= args.min_single_count]
        )
    ]

    # CALCULATE PERCENTAGES

    def calculate_percentages(_df):
        _df = (
            _df.set_index(["region", "collection_week"])
            .join(location_counts.set_index(["region", "collection_week"]))
            .reset_index()
        )

        _df["percent"] = _df["counts"] / _df["location_counts"]
        _df["collection_week"] = _df["collection_week"].dt.start_time

        return _df

    lineage_counts = calculate_percentages(lineage_counts)
    spike_combo_snv_counts = calculate_percentages(spike_combo_snv_counts)
    spike_single_snv_counts = calculate_percentages(spike_single_snv_counts)

    # SNV IDS TO SNV NAMES
    def snv_ids_to_name(ids):
        if len(ids) == 0:
            return "Reference"

        snvs = []
        for snv_id in ids:
            split = spike_snv_map[snv_id].split("|")
            # Tuple of position (for sorting), and the pretty SNV name
            snvs.append((int(split[1]), split[2] + split[1] + split[3]))

        # Sort by position
        snvs = sorted(snvs, key=lambda x: x[0])

        return ",".join([snv[1] for snv in snvs])

    def snv_id_to_name(snv_id):
        if snv_id == -1:
            return "Reference"

        split = spike_snv_map[snv_id].split("|")
        # Tuple of position (for sorting), and the pretty SNV name
        return split[2] + split[1] + split[3]

    spike_combo_snv_counts.loc[:, "group"] = spike_combo_snv_counts["group"].apply(
        snv_ids_to_name
    )
    spike_single_snv_counts.loc[:, "group"] = spike_single_snv_counts["group"].apply(
        snv_id_to_name
    )

    lineage_counts.insert(0, "type", "lineage")
    spike_combo_snv_counts.insert(0, "type", "spike_combo")
    spike_single_snv_counts.insert(0, "type", "spike_single")
    all_counts = pd.concat(
        [lineage_counts, spike_combo_snv_counts, spike_single_snv_counts],
        axis=0,
        ignore_index=True,
    )
    all_counts.to_csv(out_path / "group_counts2.csv", index=False)

    # DO REGRESSIONS

    def group_regression(_df):
        min_date = _df["collection_week"].dt.date.min().isoformat()
        _df = _df.set_index("collection_week").reindex(
            pd.date_range(min_date, end_date_iso, freq="7D", closed="left"),
            fill_value=0,
        )

        if len(_df) == 0:
            return (0, 0, 1, 0, 0)

        slope, intercept, r, pval, err = linregress(
            np.arange(0, len(_df)), _df["percent"]
        )

        # if r > 0.5:
        #    print(lineage, slope, r, pval, _df['counts'].sum())

        return (slope, r, pval, _df["counts"].sum(), _df["percent"].max())

    def calculate_trends(_df):
        regression_df = (
            _df.groupby(["region", "group"])[["collection_week", "counts", "percent"]]
            .apply(group_regression)
            .rename("res")
            .reset_index()
        )

        (
            regression_df["slope"],
            regression_df["r"],
            regression_df["pval"],
            regression_df["counts"],
            regression_df["max_percent"],
        ) = zip(*regression_df["res"])
        regression_df.drop(columns=["res"], inplace=True)

        return regression_df

    lineage_regression = calculate_trends(lineage_counts)
    spike_combo_snv_regression = calculate_trends(spike_combo_snv_counts)
    spike_single_snv_regression = calculate_trends(spike_single_snv_counts)

    lineage_regression.insert(0, "type", "lineage")
    spike_combo_snv_regression.insert(0, "type", "spike_combo")
    spike_single_snv_regression.insert(0, "type", "spike_single")
    all_regression = pd.concat(
        [lineage_regression, spike_combo_snv_regression, spike_single_snv_regression],
        axis=0,
        ignore_index=True,
    )
    all_regression.to_csv(out_path / "group_regression2.csv", index=False)


if __name__ == "__main__":
    main()
