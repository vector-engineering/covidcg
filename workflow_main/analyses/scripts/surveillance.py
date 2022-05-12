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
        "--metadata", type=str, required=True, help="Path to metadata CSV file",
    )

    parser.add_argument(
        "--group-col", type=str, required=True, help="Column name for grouping"
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

    case_data = pd.read_csv(args.metadata)

    out_path = Path(args.output)

    df = case_data[["Accession ID", "collection_date", args.group_col, "region"]]

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

    # LINEAGE DATA
    lineage_counts = (
        df.loc[df["collection_date"] >= pd.to_datetime(start_date_iso)]
        .groupby(["region", args.group_col, "collection_week"])
        .size()
        .reset_index()
        .rename(columns={args.group_col: "group", 0: "counts"})
    )

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

    # lineage_counts.insert(0, "type", "lineage")
    lineage_counts.to_json(out_path / "group_counts2.json", orient="records")

    # DO REGRESSIONS

    def group_regression(_df):
        min_date = _df["collection_week"].dt.date.min().isoformat()
        _df = _df.set_index("collection_week").reindex(
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

    # lineage_regression.insert(0, "type", "lineage")
    lineage_regression.to_json(out_path / "group_regression2.json", orient="records")


if __name__ == "__main__":
    main()
