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
        "--isolate-data", type=str, required=True, help="Path to isolate data CSV file",
    )
    parser.add_argument(
        "--metadata-map", type=str, required=True, help="Metadata map JSON file"
    )

    parser.add_argument(
        "--group-col", type=str, required=True, help="Column name for grouping"
    )
    parser.add_argument(
        "--group-references",
        type=str,
        nargs="*",
        default=None,
        help="References to use for each group",
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
        "--period",
        type=str,
        default="W",
        help="Aggregation period. W = week, M = month, Y = year. Default: W",
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

    isolate_df = pd.read_csv(
        args.isolate_data,
        usecols=[
            "isolate_id",
            "reference",
            "collection_date",
            args.group_col,
            "region",
        ],
    )
    out_path = Path(args.output)

    with open(args.metadata_map, "r") as fp:
        metadata_map = json.loads(fp.read())

    # Join locations
    isolate_df.loc[:, "region"] = isolate_df["region"].map(
        {int(k): v for k, v in metadata_map["region"].items()}
    )
    isolate_df.loc[isolate_df["region"].isna(), "region"] = None

    # Filter for group references, if defined
    if args.group_references:
        valid_group_reference_pair = pd.Series(
            data=False, index=isolate_df.index.values
        )
        for pair in args.group_references:
            group = pair.split("=")[0]
            reference = pair.split("=")[1]
            valid_group_reference_pair = valid_group_reference_pair | (
                (isolate_df[args.group_col] == group)
                & (isolate_df["reference"] == reference)
            )

        isolate_df = isolate_df.loc[valid_group_reference_pair, :]

    # Filter for valid regions
    valid_regions = [
        "Africa",
        "Asia",
        "Europe",
        "North America",
        "Oceania",
        "South America",
    ]
    isolate_df = isolate_df.loc[isolate_df["region"].isin(valid_regions)]

    isolate_df["collection_date"] = pd.to_datetime(isolate_df["collection_date"])
    isolate_df["collection_period"] = isolate_df["collection_date"].dt.to_period(
        args.period
    )

    location_counts = (
        isolate_df.groupby(["region", "collection_period"])
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

    # LINEAGE DATA
    lineage_counts = (
        isolate_df.loc[isolate_df["collection_date"] >= pd.to_datetime(start_date_iso)]
        .groupby(["region", args.group_col, "collection_period"])
        .size()
        .reset_index()
        .rename(columns={args.group_col: "group", 0: "counts"})
    )

    # CALCULATE PERCENTAGES

    def calculate_percentages(_df):
        _df = (
            _df.set_index(["region", "collection_period"])
            .join(location_counts.set_index(["region", "collection_period"]))
            .reset_index()
        )

        _df["percent"] = _df["counts"] / _df["location_counts"]
        _df["collection_period"] = _df["collection_period"].dt.start_time

        return _df

    lineage_counts = calculate_percentages(lineage_counts)

    # lineage_counts.insert(0, "type", "lineage")
    lineage_counts.to_json(out_path / "group_counts2.json", orient="records", indent=2)

    # DO REGRESSIONS

    def group_regression(_df):
        min_date = _df["collection_period"].dt.date.min().isoformat()
        _df = _df.set_index("collection_period").reindex(
            pd.date_range(min_date, end_date_iso, freq="7D", closed="left"),
            fill_value=0,
        )

        if len(_df) == 0:
            return (0, 0, 1, 0, 0)

        slope, _, r, pval, _ = linregress(np.arange(0, len(_df)), _df["percent"])

        # if r > 0.5:
        #    print(lineage, slope, r, pval, _df['counts'].sum())

        return (slope, r, pval, _df["counts"].sum(), _df["percent"].max())

    def calculate_trends(_df):
        regression_df = (
            _df.groupby(["region", "group"])[["collection_period", "counts", "percent"]]
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

    # lineage_regression.insert(0, "type", "lineage")
    lineage_regression.to_json(
        out_path / "group_regression2.json", orient="records", indent=2
    )


if __name__ == "__main__":
    main()
