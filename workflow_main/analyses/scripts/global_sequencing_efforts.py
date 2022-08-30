# coding: utf-8

"""Merge sequence data with case data, calculate sequencing efforts per country
and merge in some geological data

Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""

import argparse
import json
import numpy as np
import pandas as pd


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--isolate-data", type=str, required=True, help="Path to isolate data CSV file",
    )
    parser.add_argument(
        "--metadata-map", type=str, required=True, help="Metadata map JSON file"
    )
    parser.add_argument(
        "--output", type=str, required=True, help="Path to output JSON file",
    )

    args = parser.parse_args()

    # Load case counts by country
    case_count_df = pd.read_csv(
        "https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv"
    )

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
    case_count_df.loc[upgrade_province_inds, "Country/Region"] = case_count_df.loc[
        upgrade_province_inds, "Province/State"
    ]

    # Group by country/region
    case_count_df = (
        case_count_df.drop(columns=["Lat", "Long"])
        .groupby("Country/Region")
        .agg(np.sum)
        .reset_index()
    )
    # Unpivot table
    case_count_df = pd.melt(
        case_count_df,
        id_vars=["Country/Region"],
        var_name="date",
        value_name="cumulative_cases",
    )
    # Convert date strings to datetime objects
    case_count_df["date"] = pd.to_datetime(case_count_df["date"])
    case_count_df["month"] = case_count_df["date"].dt.to_period("M")

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
    isolate_df["year_month"] = isolate_df["collection_date"].dt.to_period("M")

    # Remove invalid submission dates (negative turnaround times)
    isolate_df = isolate_df.loc[isolate_df["turnaround_days"] >= 0]

    # Upgrade provinces to countries
    upgrade_inds = isolate_df["division"].isin(upgrade_provinces)
    isolate_df.loc[upgrade_inds, "country"] = isolate_df.loc[upgrade_inds, "division"]

    # Load UID ISO FIPS lookup table
    iso_lookup_df = pd.read_csv(
        "https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/UID_ISO_FIPS_LookUp_Table.csv"
    )
    # Upgrade provinces to country/regions
    upgrade_inds = iso_lookup_df["Province_State"].isin(upgrade_provinces)
    iso_lookup_df.loc[upgrade_inds, "Country_Region"] = iso_lookup_df.loc[
        upgrade_inds, "Province_State"
    ]

    # Only take countries, then set as the index
    iso_lookup_df = (
        iso_lookup_df.loc[
            (upgrade_inds & pd.isnull(iso_lookup_df["Admin2"]))
            | (pd.isnull(iso_lookup_df["Province_State"]))
        ]
        .set_index("Country_Region")
        .rename(
            {
                "US": "USA",
                "Congo (Kinshasa)": "Democratic Republic of the Congo",
                "Congo (Brazzaville)": "Republic of the Congo",
                "Korea, South": "South Korea",
                "Taiwan*": "Taiwan",
                "Burma": "Myanmar",
            }
        )
    )

    # Combine everything together
    country_df = (
        isolate_df
        # .loc[
        #     (nextmeta_df["date"] > pd.to_datetime("2020-01-01")) &
        #     (nextmeta_df["date"] < pd.to_datetime("2020-07-01"))
        # ]
        .reset_index()
        .groupby("country")
        .agg(
            median_turnaround_days=pd.NamedAgg(
                column="turnaround_days", aggfunc=np.median
            ),
            min_turnaround_days=pd.NamedAgg(column="turnaround_days", aggfunc=np.min),
            max_turnaround_days=pd.NamedAgg(column="turnaround_days", aggfunc=np.max),
            num_sequences=pd.NamedAgg(column="isolate_id", aggfunc="count"),
        )
        .rename({"Palestine": "West Bank and Gaza"})
        .join(
            case_count_df.groupby("Country/Region")["cumulative_cases"]
            .agg(np.max)
            .rename(
                {
                    "US": "USA",
                    "Congo (Kinshasa)": "Democratic Republic of the Congo",
                    "Congo (Brazzaville)": "Republic of the Congo",
                    "Korea, South": "South Korea",
                    "Taiwan*": "Taiwan",
                    "Burma": "Myanmar",
                }
            )
        )
        .join(iso_lookup_df, how="right")
        .reset_index()
        .rename(columns={"index": "country", "cumulative_cases": "cases"})
    )

    # Fill some column"s missing values with 0
    country_df["num_sequences"] = country_df["num_sequences"].fillna(0)
    country_df["sequences_per_case"] = (
        country_df["num_sequences"] / country_df["cases"]
    ).fillna(0)

    # Only take some columns
    country_df = country_df.loc[
        :,
        [
            "UID",
            "Country_Region",
            "median_turnaround_days",
            "min_turnaround_days",
            "max_turnaround_days",
            "num_sequences",
            "cases",
            "sequences_per_case",
        ],
    ]

    # Write to disk
    # First write JSON to string
    country_df_str = country_df.to_json(orient="records")
    # Manually add some missing records
    country_df_str = country_df_str[:-1] + (
        ',{"UID":260,"Country_Region":"Fr. S. Antarctic Lands","median_turnaround_days":null,"min_turnaround_days":null,"max_turnaround_days":null,"num_sequences":null,"cases":null,"sequences_per_case":null}'
        + ',{"UID":795,"Country_Region":"Turkmenistan","median_turnaround_days":null,"min_turnaround_days":null,"max_turnaround_days":null,"num_sequences":null,"cases":null,"sequences_per_case":null}'
        + ',{"UID":10,"Country_Region":"Antarctica","median_turnaround_days":null,"min_turnaround_days":null,"max_turnaround_days":null,"num_sequences":null,"cases":null,"sequences_per_case":null}'
        + ',{"UID":408,"Country_Region":"North Korea","median_turnaround_days":null,"min_turnaround_days":null,"max_turnaround_days":null,"num_sequences":null,"cases":null,"sequences_per_case":null}'
        + ',{"UID":90,"Country_Region":"Solomon Islands","median_turnaround_days":null,"min_turnaround_days":null,"max_turnaround_days":null,"num_sequences":null,"cases":null,"sequences_per_case":null}'
        + ',{"UID":548,"Country_Region":"Vanuatu","median_turnaround_days":null,"min_turnaround_days":null,"max_turnaround_days":null,"num_sequences":null,"cases":null,"sequences_per_case":null}'
        +
        # GISAID really wants French Guiana separate from france, so in my custom geojson I made French Guiana ID: -98
        ',{"UID":-98,"Country_Region":"French Guiana","median_turnaround_days":null,"min_turnaround_days":null,"max_turnaround_days":null,"num_sequences":null,"cases":null,"sequences_per_case":null}'
        +
        # Northern Cyprus
        ',{"UID":-99,"Country_Region":"Northern Cyprus","median_turnaround_days":null,"min_turnaround_days":null,"max_turnaround_days":null,"num_sequences":null,"cases":null,"sequences_per_case":null}'
        + "]"
    )
    with open(args.output, "w") as fp:
        fp.write(country_df_str)


if __name__ == "__main__":
    main()
