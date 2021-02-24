#!/usr/bin/env python3
# coding: utf-8

"""Global sequencing data for the home page

Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""

import argparse
import pandas as pd
import numpy as np
import json

from pathlib import Path


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--case-data", type=str, required=True, help="Path to case data CSV file",
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

    out_path = Path(args.output)

    # Load case counts by country
    case_count_df = pd.read_csv(
        "https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv"
    )
    case_count_df.rename(columns={"Country/Region": "country"}, inplace=True)

    # Upgrade some province/states to country/regions
    upgrade_provinces = [
        "Hong Kong",
        "Macau",
        "Faroe Islands",
        "Greenland",
        "French Guiana",
        "French Polynesia",
        "Guadeloupe",
        "Martinique",
        "Mayotte",
        "New Caledonia",
        "Reunion",
        "Saint Barthelemy",
        "Saint Pierre and Miquelon",
        "St Martin",
        "Aruba",
        "Bonaire, Sint Eustatius and Saba",
        "Curacao",
        "Sint Maarten",
        "Anguilla",
        "Bermuda",
        "British Virgin Islands",
        "Cayman Islands",
        "Falkland Islands (Malvinas)",
        "Gibraltar",
        "Isle of Man",
        "Channel Islands",
        "Montserrat",
        "Turks and Caicos Islands",
        "American Samoa",
        "Guam",
        "Northern Mariana Islands",
        "Virgin Islands",
        "Puerto Rico",
    ]
    upgrade_province_inds = case_count_df["Province/State"].isin(upgrade_provinces)
    case_count_df.loc[upgrade_province_inds, "country"] = case_count_df.loc[
        upgrade_province_inds, "Province/State"
    ]

    # Group by country/region
    case_count_df = (
        case_count_df.drop(columns=["Lat", "Long"])
        .groupby("country")
        .agg(np.sum)
        .reset_index()
    )
    # Unpivot table
    case_count_df = pd.melt(
        case_count_df,
        id_vars=["country"],
        var_name="date",
        value_name="cumulative_cases",
    )
    # Convert date strings to datetime objects
    case_count_df["date"] = pd.to_datetime(case_count_df["date"])
    case_count_df["month"] = case_count_df["date"].dt.to_period("M")

    JHU_rename_map = {
        "US": "USA",
        "Congo (Kinshasa)": "DRC",
        "Congo (Brazzaville)": "Republic of the Congo",
        "Korea, South": "South Korea",
        "Taiwan*": "Taiwan",
        "Burma": "Myanmar",
        #    "Aruba": "Netherlands",
        #    "Bonaire, Sint Eustatius and Saba": "Netherlands",
        #    "Curacao": "Netherlands",
        #    "Sint Maarten": "Netherlands",
        #    "British Virgin Islands": "United Kingdom",
        #    "Channel Islands": "United Kingdom",
        #    "Cayman Islands": "United Kingdom",
        #    "Gibraltar": "United Kingdom",
        #    "Isle of Man": "United Kingdom",
        #    "Montserrat": "United Kingdom",
        #    "Turks and Caicos Islands": "United Kingdom",
        #    "Falkland Islands (Malvinas)": "United Kingdom",
        #    "Diamond Princess": "Japan",
        #    "Faroe Islands": "Denmark",
        #    "French Polynesia": "France",
        #    "Guadeloupe": "France",
        #    "Martinique": "France",
        #    "Mayotte": "France",
        #    "Reunion": "France",
        #    "New Caledonia": "France",
        #    "Saint Barthelemy": "France",
        #    "Saint Pierre and Miquelon": "France",
        #    "St Martin": "France",
        #    "St Martin": "Saint Martin",
        #    "MS Zaandam": "USA",
        #    "Marshall Islands": "USA",
        #    "Macau": "China",
    }

    def rename_countries(country):
        if country in JHU_rename_map.keys():
            return JHU_rename_map[country]
        else:
            return country

    case_count_df["country"] = case_count_df["country"].apply(rename_countries)

    case_count_df = (
        case_count_df.groupby(["country", "month"])["cumulative_cases"]
        .agg(np.max)
        .reset_index()
    )

    case_count_df["month"] = case_count_df["month"].dt.start_time
    case_count_df.to_json(str(out_path / "case_count.json"), orient="records")

    case_df = pd.read_csv(args.case_data, index_col="Accession ID")
    case_df = case_df[["collection_date", "submission_date", "location_id"]]
    location_map = pd.read_json(args.location_map)
    case_df = case_df.join(location_map, on="location_id", how="left")

    case_df["collection_date"] = pd.to_datetime(
        case_df["collection_date"], errors="coerce"
    )
    case_df["submission_date"] = pd.to_datetime(
        case_df["submission_date"], errors="coerce"
    )

    # Remove failed date parsing
    case_df = case_df.loc[
        (~pd.isnull(case_df["collection_date"]))
        & (~pd.isnull(case_df["submission_date"]))
    ]

    # Only take dates from 2019-12-15
    case_df = case_df.loc[case_df["collection_date"] > pd.to_datetime("2019-12-15")]

    # Calculate time deltas
    case_df["turnaround_days"] = (
        case_df["submission_date"] - case_df["collection_date"]
    ).dt.days
    # Extract month
    case_df["month"] = case_df["collection_date"].dt.to_period("M")
    case_df["submission_month"] = case_df["submission_date"].dt.to_period("M")

    # Remove invalid submission dates (negative turnaround times)
    case_df = case_df.loc[case_df["turnaround_days"] >= 0]

    # Upgrade provinces to countries
    upgrade_inds = case_df["division"].isin(upgrade_provinces)
    case_df.loc[upgrade_inds, "country"] = case_df.loc[upgrade_inds, "division"]

    sequences_per_month = (
        case_df.reset_index()
        .groupby(["country", "month"])["Accession ID"]
        .size()
        .rename({"Palestine": "West Bank and Gaza"})
        .rename("new_sequences")
        .reset_index()
    )
    sequences_per_month["month"] = sequences_per_month["month"].dt.start_time
    sequences_per_month.to_json(
        str(out_path / "sequences_per_month.json"), orient="records"
    )

    turnaround_per_month = (
        case_df.reset_index()
        .groupby(["country", "submission_month"])["turnaround_days"]
        .agg(
            q5=lambda x: np.quantile(x, 0.05),
            q25=lambda x: np.quantile(x, 0.25),
            q50=lambda x: np.quantile(x, 0.50),
            q75=lambda x: np.quantile(x, 0.75),
            q95=lambda x: np.quantile(x, 0.95),
        )
        .reset_index()
    )
    turnaround_per_month["submission_month"] = turnaround_per_month[
        "submission_month"
    ].dt.start_time

    turnaround_per_month.to_json(
        str(out_path / "turnaround_per_month.json"), orient="records"
    )

    # Load UID ISO FIPS lookup table
    iso_lookup_df = pd.read_csv(
        "https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/UID_ISO_FIPS_LookUp_Table.csv"
    )
    # Upgrade provinces to country/regions
    upgrade_inds = iso_lookup_df["Province_State"].isin(upgrade_provinces)
    iso_lookup_df.rename(columns={"Country_Region": "country"}, inplace=True)
    iso_lookup_df.loc[upgrade_inds, "country"] = iso_lookup_df.loc[
        upgrade_inds, "Province_State"
    ]

    # Only take countries, then set as the index
    iso_lookup_df = (
        iso_lookup_df.loc[
            (upgrade_inds & pd.isnull(iso_lookup_df["Admin2"]))
            | (pd.isnull(iso_lookup_df["Province_State"]))
        ]
        .set_index("country")
        .rename(JHU_rename_map)
        .drop(columns=["Combined_Key", "Province_State", "Admin2", "FIPS"])
    )

    # New Caledonia -> France
    # Falklands -> UK

    add_countries = [
        {"UID": 260, "country": "Fr. S. Antarctic Lands"},
        {"UID": 795, "country": "Turkmenistan"},
        {"UID": 10, "country": "Antarctica"},
        {"UID": 408, "country": "North Korea"},
        {"UID": 90, "country": "Solomon Islands"},
        {"UID": 548, "country": "Vanuatu"},
        # GISAID really wants French Guiana separate from France,
        # so in my custom geojson I made French Guiana ID: -98
        {"UID": -98, "country": "French Guiana"},
        # Northern Cyprus
        {"UID": -99, "country": "Northern Cyprus"},
    ]

    iso_lookup_df = pd.concat(
        [
            iso_lookup_df,
            pd.DataFrame.from_records(
                add_countries,
                columns=["country"] + iso_lookup_df.columns.tolist(),
                index="country",
            ),
        ],
        axis=0,
    )

    iso_lookup_df.reset_index().to_json(
        str(out_path / "iso_lookup.json"), orient="records"
    )


if __name__ == "__main__":
    main()

