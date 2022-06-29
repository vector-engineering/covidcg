#!/usr/bin/env python3
# coding: utf-8

"""Global sequencing data for the home page

Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""

import argparse
import json
import pandas as pd
import numpy as np

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

    isolate_df = pd.read_csv(
        args.isolate_data,
        usecols=[
            "isolate_id",
            "collection_date",
            "submission_date",
            "country",
            "division",
        ],
    )

    isolate_df.drop_duplicates("isolate_id", inplace=True)
    isolate_df.set_index("isolate_id", inplace=True)

    # Load location metadata mappings
    with open(args.metadata_map, "r") as fp:
        metadata_map = json.loads(fp.read())

    # Join locations onto isolate_df
    loc_levels = ["country", "division"]
    for loc_level in loc_levels:
        isolate_df.loc[:, loc_level] = isolate_df[loc_level].map(
            {int(k): v for k, v in metadata_map[loc_level].items()}
        )
        isolate_df.loc[isolate_df[loc_level].isna(), loc_level] = None

    isolate_df["collection_date"] = pd.to_datetime(
        isolate_df["collection_date"], errors="coerce"
    )
    isolate_df["submission_date"] = pd.to_datetime(
        isolate_df["submission_date"], errors="coerce"
    )

    # Remove failed date parsing
    isolate_df = isolate_df.loc[
        (~pd.isnull(isolate_df["collection_date"]))
        & (~pd.isnull(isolate_df["submission_date"]))
    ]

    # Only take dates from 2019-12-15
    isolate_df = isolate_df.loc[
        isolate_df["collection_date"] > pd.to_datetime("2019-12-15")
    ]

    # Calculate time deltas
    isolate_df["turnaround_days"] = (
        isolate_df["submission_date"] - isolate_df["collection_date"]
    ).dt.days
    # Extract month
    isolate_df["month"] = isolate_df["collection_date"].dt.to_period("M")
    isolate_df["submission_month"] = isolate_df["submission_date"].dt.to_period("M")

    # Remove invalid submission dates (negative turnaround times)
    isolate_df = isolate_df.loc[isolate_df["turnaround_days"] >= 0]

    # Upgrade provinces to countries
    upgrade_inds = isolate_df["division"].isin(upgrade_provinces)
    isolate_df.loc[upgrade_inds, "country"] = isolate_df.loc[upgrade_inds, "division"]

    sequences_per_month = (
        isolate_df.reset_index()
        .groupby(["country", "month"])["isolate_id"]
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
        isolate_df.reset_index()
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

