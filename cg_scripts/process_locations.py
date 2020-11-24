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
from cg_scripts.util import human_format


def process_location_metadata(case_df, location_corretions):

    # Unset index for faster filtering
    case_df = case_df.reset_index()

    print("Processing location data...", end="", flush=True)
    # Location data is stored in one column, "region / country / division / location"
    location_df = (
        case_df["covv_location"]
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

    # Create complete location column from the separate parts
    # This time with no padding
    location_df["loc_str"] = location_df["region"].str.cat(
        [
            location_df["country"].astype(str),
            location_df["division"].astype(str),
            location_df["location"].astype(str),
        ],
        sep="/",
    )

    location_map_df = (
        location_df["loc_str"]
        .drop_duplicates()
        .sort_values(ignore_index=True)
        .reset_index()  # Produce an index column
    )
    # Location data is stored in one column, "region/country/division/location"
    location_map_df = pd.concat(
        [
            location_map_df,
            (
                location_map_df["loc_str"]
                .str.split("/", expand=True)
                .iloc[:, :4]  # Only take 4 columns
                # Rename columns
                .rename(
                    columns={0: "region", 1: "country", 2: "division", 3: "location"}
                )
            ),
        ],
        axis=1,
    )

    # Map location IDs back to taxon_locations dataframe
    location_df["location_id"] = location_df["loc_str"].map(
        pd.Series(
            location_map_df["index"].values, index=location_map_df["loc_str"].values,
        )
    )

    # Take subset of columns, re-index
    location_df = location_df[
        ["location_id", "region", "country", "division", "location"]
    ].reset_index(drop=True)
    print("done")

    # Re-apply Accession ID index
    location_df["Accession ID"] = case_df["Accession ID"]
    location_df = location_df.set_index("Accession ID")

    return location_df, location_map_df


def clean_location_data(location_df, location_corretions):
    """Fix typos, unify nomenclature in location data
    """

    # Load rules
    location_correction_df = pd.read_csv(location_corretions)
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


def build_select_tree(
    location_df, location_map_df, emoji_map_file="static_data/country_to_emoji.xls"
):
    """Build tree for ReactDropdownTreeSelect

    data
    Type: Object or Array

    Data for rendering the tree select items. The object requires the following structure:

    {
    label,          // required: Checkbox label
    value,          // required: Checkbox value
    children,       // optional: Array of child objects
    checked,        // optional: Initial state of checkbox. if true, checkbox is selected and corresponding pill is rendered.
    disabled,       // optional: Selectable state of checkbox. if true, the checkbox is disabled and the node is not selectable.
    expanded,       // optional: If true, the node is expanded (children of children nodes are not expanded by default unless children nodes also have expanded: true).
    className,      // optional: Additional css class for the node. This is helpful to style the nodes your way
    tagClassName,   // optional: Css class for the corresponding tag. Use this to add custom style the pill corresponding to the node.
    actions,        // optional: An array of extra action on the node (such as displaying an info icon or any custom icons/elements)
    dataset,        // optional: Allows data-* attributes to be set on the node and tag elements
    isDefaultValue, // optional: Indicate if a node is a default value. When true, the dropdown will automatically select the node(s) when there is no other selected node. Can be used on more than one node.
    ...             // optional: Any extra properties that you'd like to receive during `onChange` event
    }
    The action object requires the following structure:

    {
    className, // required: CSS class for the node. e.g. `fa fa-info`
    title,     // optional: HTML tooltip text
    text,      // optional: Any text to be displayed. This is helpful to pass ligatures if you're using ligature fonts
    ...        // optional: Any extra properties that you'd like to receive during `onChange` event
    }
    An array renders a tree with multiple root level items whereas an object renders a tree with a single root element (e.g. a Select All root node).


    Example:
    const data = {
    label: 'search me',
    value: 'searchme',
    children: [
        {
        label: 'search me too',
        value: 'searchmetoo',
        children: [
            {
            label: 'No one can get me',
            value: 'anonymous',
            },
        ],
        },
    ],
    }
    """

    # Set unspecified locations to None so that they don't get
    # caught up in the groupby
    location_df.loc[location_df["region"] == "-1", "region"] = None
    location_df.loc[location_df["country"] == "-1", "country"] = None
    location_df.loc[location_df["division"] == "-1", "division"] = None
    location_df.loc[location_df["location"] == "-1", "location"] = None

    # Count sequences per grouping level
    region_counts = dict(location_df.groupby("region")["location_id"].count())
    country_counts = dict(
        location_df.groupby(["region", "country"])["location_id"].count()
    )
    division_counts = dict(
        location_df.groupby(["region", "country", "division"])["location_id"].count()
    )
    location_counts = dict(
        location_df.groupby(["region", "country", "division", "location"])[
            "location_id"
        ].count()
    )

    # Load country -> emoji map
    emoji_map = pd.read_excel(emoji_map_file, skiprows=1)
    # Expand country aliases, remove whitespace from each alias
    emoji_map["aliases"] = (
        emoji_map["aliases"].str.split(",").apply(lambda x: [y.strip() for y in x])
    )

    # Root node
    select_tree = {"label": "All", "value": "All", "children": []}

    for i, loc in location_map_df.iterrows():
        # Add region node
        if loc["region"] == "-1":
            continue

        region_node = [
            c for c in select_tree["children"] if c["value"] == loc["region"]
        ]
        if region_node:
            region_node = region_node[0]
        else:
            region_node = {
                "label": loc["region"],
                "value": loc["region"],
                "level": "region",
                "location_id": loc["index"],
                "actions": [
                    {
                        "className": "fa fa-info",
                        "title": str(region_counts[loc["region"]]) + " sequences",
                        "text": human_format(region_counts[loc["region"]]),
                    }
                ],
                "children": [],
            }
            select_tree["children"].append(region_node)

        # Add country --> region
        if loc["country"] == "-1":
            continue

        country_node = [
            c for c in region_node["children"] if c["value"] == loc["country"]
        ]
        if country_node:
            country_node = country_node[0]
        else:

            # Look for an emoji for this country
            country_emoji = ""
            emoji_entry = emoji_map.loc[
                emoji_map["aliases"].apply(lambda x: loc["country"] in x), :
            ]
            # Fill the country emoji, if it exists
            if len(emoji_entry) == 1:
                country_emoji = emoji_entry.iat[0, 1] + " "

            country_node = {
                "label": country_emoji + loc["country"],
                "value": loc["country"],
                "region": loc["region"],
                "level": "country",
                "location_id": loc["index"],
                "actions": [
                    {
                        "className": "fa fa-info",
                        "title": str(country_counts[(loc["region"], loc["country"])])
                        + " sequences",
                        "text": human_format(
                            country_counts[(loc["region"], loc["country"])]
                        ),
                    }
                ],
                "children": [],
            }
            region_node["children"].append(country_node)

        # Add division --> country
        if loc["division"] == "-1":
            continue

        division_node = [
            c for c in country_node["children"] if c["value"] == loc["division"]
        ]
        if division_node:
            division_node = division_node[0]
        else:
            division_node = {
                "label": loc["division"],
                "value": loc["division"],
                "region": loc["region"],
                "country": loc["country"],
                "level": "division",
                "location_id": loc["index"],
                "actions": [
                    {
                        "className": "fa fa-info",
                        "title": str(
                            division_counts[
                                (loc["region"], loc["country"], loc["division"])
                            ]
                        )
                        + " sequences",
                        "text": human_format(
                            division_counts[
                                (loc["region"], loc["country"], loc["division"])
                            ]
                        ),
                    }
                ],
                "children": [],
            }
            country_node["children"].append(division_node)

        # Add location --> division
        if loc["location"] == "-1":
            continue

        location_node = [
            c for c in division_node["children"] if c["value"] == loc["location"]
        ]
        if location_node:
            location_node = location_node[0]
        else:
            location_node = {
                "label": loc["location"],
                "value": loc["location"],
                "region": loc["region"],
                "country": loc["country"],
                "division": loc["division"],
                "level": "location",
                "location_id": loc["index"],
                "actions": [
                    {
                        "className": "fa fa-info",
                        "title": str(
                            location_counts[
                                (
                                    loc["region"],
                                    loc["country"],
                                    loc["division"],
                                    loc["location"],
                                )
                            ]
                        )
                        + " sequences",
                        "text": human_format(
                            location_counts[
                                (
                                    loc["region"],
                                    loc["country"],
                                    loc["division"],
                                    loc["location"],
                                )
                            ]
                        ),
                    }
                ],
                "children": [],
            }
            division_node["children"].append(location_node)

    # print(loc_tree.nodes)
    # print(loc_tree.in_edges('New York'))
    return select_tree


# if __name__ == "__main__":
#     patient_meta_df = load_patient_metadata()
#     location_df, location_map_df = process_location_metadata(patient_meta_df)
