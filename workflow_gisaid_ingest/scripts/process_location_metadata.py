#!/usr/bin/env python3
# coding: utf-8

"""Clean location metadata
Build hierarchical location tree for the UI

Author: Albert Chen (Deverman Lab - Broad Institute)
"""

import json
import numpy as np
import pandas as pd

from functools import reduce


def process_location_metadata(df, location_corretions):

    print("Processing location data...", end="", flush=True)
    # Location data is stored in one column, "region / country / division / location"
    location_df = (
        # Unset index for faster filtering
        df.reset_index()["covv_location"]
        .str.split("/", expand=True)
        .iloc[:, :4]  # Only take 4 columns
        # Rename columns
        .rename(columns={0: "region", 1: "country", 2: "division", 3: "location"})
        .applymap(lambda x: x.strip() if x else x)
        # Placeholder for missing values, so that it will still
        # be caught by groupby() later on
        .fillna(-1)
    )

    # With sparse data, sometimes none of the sequences fed in will have a
    # "division" or "location" entry.
    # Make them manually now, if they don't already exist
    if "division" not in location_df.columns:
        location_df["division"] = -1
    if "location" not in location_df.columns:
        location_df["location"] = -1

    # Clean location data
    location_df = clean_location_data(location_df, location_corretions)

    # Re-create index for the join
    location_df["Accession ID"] = df.index.values
    location_df = location_df.set_index("Accession ID")

    # Join values to master dataframe
    # Do an inner join... but none of the rows should've been filtered out
    # so this doesn't really matter
    df = df.join(location_df, how="inner")

    print("done")

    return df


def clean_location_data(location_df, location_corretions):
    """Fix typos, unify nomenclature in location data
    """

    # Load rules
    location_correction_df = pd.read_csv(location_corretions, comment="#")
    # region_pattern,country_pattern,division_pattern,location_pattern,out_region,out_country,out_division,out_location,comment

    for i, rule in location_correction_df.iterrows():
        # print(rule)
        input_rule = {
            "region": rule["region_pattern"],
            "country": rule["country_pattern"],
            "division": rule["division_pattern"],
            "location": rule["location_pattern"],
        }
        output_rule = {
            "region": rule["out_region"],
            "country": rule["out_country"],
            "division": rule["out_division"],
            "location": rule["out_location"],
        }

        # Get matching entries for the input rule
        # by creating a logical mask
        # Start out with matching everything
        loc_mask = pd.Series(np.repeat(True, len(location_df)))
        for key in input_rule.keys():
            if type(input_rule[key]) is not str or not input_rule[key]:
                continue

            vals = input_rule[key].split("|")
            # Make it a list if it's just a single value
            if type(vals) is not list:
                vals = [vals]
            vals = [str(val) for val in vals]

            # Turn each value into a logical mask
            vals = [location_df[key] == v for v in vals]
            # Combine logical masks with logical ORs, and merge into the master mask with AND
            loc_mask = loc_mask & reduce(lambda x, y: (x | y), vals)

        # Set the output rules on the matching entries from loc_mask
        for out_key in output_rule.keys():
            if (
                type(output_rule[out_key]) is not str
                and type(output_rule[out_key]) is not int
            ):
                continue
            location_df.loc[loc_mask, out_key] = output_rule[out_key]

    # Done
    return location_df
