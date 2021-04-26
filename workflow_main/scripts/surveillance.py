#!/usr/bin/env python3
# coding: utf-8

"""Global surveillance data for the home page

Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""

import argparse
import datetime
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
        "-o", "--output", type=str, required=True, help="Path to output directory",
    )

    args = parser.parse_args()

    case_data = pd.read_json(args.case_data)
    location_map = pd.read_json(args.location_map)
    out_path = Path(args.output)

    df = case_data[["Accession ID", "collection_date", "lineage", "location_id"]].join(
        location_map, on="location_id"
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

    start_date_iso = (datetime.date.today() - datetime.timedelta(days=90)).isoformat()
    # End date only used for the regression
    end_date_iso = (datetime.date.today() - datetime.timedelta(days=30)).isoformat()
    lineage_counts = (
        df.loc[df["collection_date"] >= pd.to_datetime(start_date_iso)]
        .groupby(["region", "lineage", "collection_week"])
        .size()
        .rename("counts")
        .reset_index()
    )

    lineage_counts = (
        lineage_counts.set_index(["region", "collection_week"])
        .join(location_counts.set_index(["region", "collection_week"]))
        .reset_index()
    )

    lineage_counts["percent"] = (
        lineage_counts["counts"] / lineage_counts["location_counts"]
    )
    lineage_counts["collection_week"] = lineage_counts["collection_week"].dt.start_time
    lineage_counts = lineage_counts.rename(columns={"lineage": "group"})

    lineage_counts.to_json(str(out_path / "group_counts.json"), orient="records")

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

    lineage_regression = (
        lineage_counts.groupby(["region", "group"])[
            ["collection_week", "counts", "percent"]
        ]
        .apply(group_regression)
        .rename("res")
        .reset_index()
    )
    (
        lineage_regression["slope"],
        lineage_regression["r"],
        lineage_regression["pval"],
        lineage_regression["counts"],
        lineage_regression["max_percent"],
    ) = zip(*lineage_regression["res"])
    lineage_regression.drop(columns=["res"], inplace=True)

    lineage_regression.to_json(
        str(out_path / "group_regression.json"), orient="records"
    )


if __name__ == "__main__":
    main()
