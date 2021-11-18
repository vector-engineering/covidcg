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
        help="Number of days before today to cut off data prior to regressions. Default: 30",
    )
    parser.add_argument(
        "--end-date",
        type=str,
        default=None,
        help="End date for filtering data in ISO format (YYYY-MM-DD). Overrides --end-date-days-ago if defined. Default: None",
    )

    parser.add_argument(
        "--min-combo-count",
        type=int,
        default=50,
        help="Minimum counts for a spike mutation combo to be included in the dataset",
    )
    parser.add_argument(
        "--min-single-count",
        type=int,
        default=50,
        help="Minimum counts for a single spike mutation to be included in the dataset",
    )

    parser.add_argument(
        "-o", "--output", type=str, required=True, help="Path to output directory",
    )

    args = parser.parse_args()

    case_data = pd.read_json(args.case_data)
    with open(args.metadata_map, "r") as fp:
        metadata_map = json.loads(fp.read())

    # Join region onto case_data
    case_data.loc[:, "region"] = case_data["region"].map(
        {int(k): v for k, v in metadata_map["region"].items()}
    )

    out_path = Path(args.output)

    # Get Spike mutations
    gene_aa_mutation = pd.DataFrame.from_dict(
        metadata_map["gene_aa_mutation"], orient="index"
    ).reset_index()
    gene_aa_mutation.columns = ["mutation", "mutation_id"]
    spike_mutations = gene_aa_mutation.loc[
        gene_aa_mutation["mutation"].str.match(r"^S")
    ]
    spike_mutation_map = spike_mutations["mutation"].to_dict()

    df = case_data[
        ["Accession ID", "collection_date", "genotype", "gene_aa_mutation_str", "region"]
    ]

    # Filter for only mutations in spike
    valid_mutation_ids = spike_mutations["mutation_id"].values
    df["spike_aa_mutation"] = df["gene_aa_mutation_str"].apply(
        lambda x: tuple([mut for mut in x if mut in valid_mutation_ids])
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
    if args.end_date:
        end_date_iso = args.end_date
    else:
        end_date_iso = (
            datetime.date.today() - datetime.timedelta(days=args.end_date_days_ago)
        ).isoformat()

    # GENOTYPE DATA
    genotype_counts = (
        df.loc[df["collection_date"] >= pd.to_datetime(start_date_iso)]
        .groupby(["region", "genotype", "collection_week"])
        .size()
        .reset_index()
        .rename(columns={"genotype": "group", 0: "counts"})
    )

    # SPIKE MUTATION COMBO DATA

    spike_combo_mutation_counts = (
        df.loc[df["collection_date"] >= pd.to_datetime(start_date_iso)]
        .groupby(["region", "spike_aa_mutation", "collection_week"])
        .size()
        .reset_index()
        .rename(columns={0: "counts", "spike_aa_mutation": "group"})
    )
    spike_combo_mutation_freq = (
        df.groupby("spike_aa_mutation").size().sort_values(ascending=False)
    )
    spike_combo_mutation_counts = spike_combo_mutation_counts.loc[
        spike_combo_mutation_counts["group"].isin(
            spike_combo_mutation_freq.index[
                spike_combo_mutation_freq >= args.min_combo_count
            ].values
        )
    ]

    # SPIKE SINGLE MUTATION DATA
    spike_single_mutation_counts = (
        df.loc[df["collection_date"] >= pd.to_datetime(start_date_iso)]
        .groupby(["region", "spike_aa_mutation", "collection_week"])
        .size()
        .reset_index()
        .rename(columns={0: "counts", "spike_aa_mutation": "group"})
    )
    spike_single_mutation_counts.loc[:, "group"] = spike_single_mutation_counts[
        "group"
    ].apply(list)
    spike_single_mutation_counts = (
        spike_single_mutation_counts.explode("group")
        .assign(group=lambda x: x["group"].fillna(-1))
        .groupby(["region", "collection_week", "group"])
        .agg(counts=("counts", np.sum))
        .reset_index()
    )
    spike_single_mutation_freq = (
        spike_single_mutation_counts.groupby("group")["counts"]
        .sum()
        .sort_values(ascending=False)
    )
    spike_single_mutation_counts = spike_single_mutation_counts.loc[
        spike_single_mutation_counts["group"].isin(
            spike_single_mutation_freq.index[
                spike_single_mutation_freq >= args.min_single_count
            ]
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

    genotype_counts = calculate_percentages(genotype_counts)
    spike_combo_mutation_counts = calculate_percentages(spike_combo_mutation_counts)
    spike_single_mutation_counts = calculate_percentages(spike_single_mutation_counts)

    # MUTATION IDS TO MUTATION NAMES
    def mutation_ids_to_name(ids):
        if len(ids) == 0:
            return "Reference"

        mutations = []
        for mutation_id in ids:
            split = spike_mutation_map[mutation_id].split("|")
            # Tuple of position (for sorting), and the pretty mutation name
            mutations.append((int(split[1]), split[2] + split[1] + split[3]))

        # Sort by position
        mutations = sorted(mutations, key=lambda x: x[0])

        return ",".join([mut[1] for mut in mutations])

    def mutation_id_to_name(mutation_id):
        if mutation_id == -1:
            return "Reference"

        split = spike_mutation_map[mutation_id].split("|")
        # Tuple of position (for sorting), and the pretty mutation name
        return split[2] + split[1] + split[3]

    spike_combo_mutation_counts.loc[:, "group"] = spike_combo_mutation_counts[
        "group"
    ].apply(mutation_ids_to_name)
    spike_single_mutation_counts.loc[:, "group"] = spike_single_mutation_counts[
        "group"
    ].apply(mutation_id_to_name)

    genotype_counts.insert(0, "type", "genotype")
    spike_combo_mutation_counts.insert(0, "type", "spike_combo")
    spike_single_mutation_counts.insert(0, "type", "spike_single")
    all_counts = pd.concat(
        [genotype_counts, spike_combo_mutation_counts, spike_single_mutation_counts],
        axis=0,
        ignore_index=True,
    )
    all_counts.to_json(out_path / "group_counts2.json", orient='records')

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
        #    print(genotype, slope, r, pval, _df['counts'].sum())

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

    genotype_regression = calculate_trends(genotype_counts)
    spike_combo_mutation_regression = calculate_trends(spike_combo_mutation_counts)
    spike_single_mutation_regression = calculate_trends(spike_single_mutation_counts)

    genotype_regression.insert(0, "type", "genotype")
    spike_combo_mutation_regression.insert(0, "type", "spike_combo")
    spike_single_mutation_regression.insert(0, "type", "spike_single")
    all_regression = pd.concat(
        [
            genotype_regression,
            spike_combo_mutation_regression,
            spike_single_mutation_regression,
        ],
        axis=0,
        ignore_index=True,
    )
    all_regression.to_json(out_path / "group_regression2.json", orient='records')


if __name__ == "__main__":
    main()
