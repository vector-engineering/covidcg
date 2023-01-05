#!/usr/bin/env python3
# coding: utf-8

"""Build tree for ReactDropdownTreeSelect

Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""

import argparse
import json
import pandas as pd

from util import human_format


def main():
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

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--isolate-data", type=str, required=True, help="Isolate data CSV file"
    )
    parser.add_argument(
        "--metadata-map", type=str, required=True, help="Metadata map file"
    )
    parser.add_argument("--emoji-map", type=str, required=True, help="Emoji map file")
    parser.add_argument("--out", type=str, required=True, help="Output file")

    args = parser.parse_args()

    loc_levels = ["region", "country", "division", "location"]
    df = pd.read_csv(
        args.isolate_data, usecols=["isolate_id"] + loc_levels
    ).drop_duplicates("isolate_id")

    with open(args.metadata_map, "r") as fp:
        metadata_map = json.loads(fp.read())

    loc_level_id_cols = [col + "_id" for col in loc_levels]
    for i, loc_level in enumerate(loc_levels):
        # Make separate column for integer ID
        df[loc_level_id_cols[i]] = df[loc_level]
        # Map integer IDs to strings
        mmap = {int(k): v for k, v in metadata_map[loc_level].items()}
        df.loc[:, loc_level] = df[loc_level].apply(
            lambda x: mmap[x] if x >= 0 else None
        )

    # Get unique locations
    unique_locations_df = df.drop_duplicates(loc_level_id_cols)

    # Count sequences per grouping level
    region_counts = dict(df.groupby("region_id")["isolate_id"].count())
    country_counts = dict(df.groupby(["region_id", "country_id"])["isolate_id"].count())
    division_counts = dict(
        df.groupby(["region_id", "country_id", "division_id"])["isolate_id"].count()
    )
    location_counts = dict(
        df.groupby(["region_id", "country_id", "division_id", "location_id"])[
            "isolate_id"
        ].count()
    )

    # Load country -> emoji map
    emoji_map = pd.read_excel(args.emoji_map, skiprows=1)
    # Expand country aliases, remove whitespace from each alias
    emoji_map["aliases"] = (
        emoji_map["aliases"].str.split(",").apply(lambda x: [y.strip() for y in x])
    )

    # Root node
    select_tree = {"label": "All", "value": "All", "value_txt": "All", "children": []}

    for i, loc in unique_locations_df.iterrows():
        # Add region node
        if loc["region_id"] == -1:
            continue

        region_node = [
            c for c in select_tree["children"] if c["value"] == loc["region_id"]
        ]
        if region_node:
            region_node = region_node[0]
        else:
            region_node = {
                "label": loc["region"],
                "value": loc["region_id"],
                "value_txt": loc["region"],
                "level": "region",
                "actions": [
                    {
                        "className": "fa fa-info",
                        "title": str(region_counts[loc["region_id"]]) + " sequences",
                        "text": human_format(region_counts[loc["region_id"]]),
                    }
                ],
                "children": [],
            }
            select_tree["children"].append(region_node)

        # Add country --> region
        if loc["country_id"] == -1:
            continue

        country_node = [
            c for c in region_node["children"] if c["value"] == loc["country_id"]
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
                "value": loc["country_id"],
                "value_txt": loc["country"],
                "region": loc["region_id"],
                "level": "country",
                "actions": [
                    {
                        "className": "fa fa-info",
                        "title": str(
                            country_counts[(loc["region_id"], loc["country_id"])]
                        )
                        + " sequences",
                        "text": human_format(
                            country_counts[(loc["region_id"], loc["country_id"])]
                        ),
                    }
                ],
                "children": [],
            }
            region_node["children"].append(country_node)

        # Add division --> country
        if loc["division_id"] == -1:
            continue

        division_node = [
            c for c in country_node["children"] if c["value"] == loc["division_id"]
        ]
        if division_node:
            division_node = division_node[0]
        else:
            division_node = {
                "label": loc["division"],
                "value": loc["division_id"],
                "value_txt": loc["division"],
                "region": loc["region_id"],
                "country": loc["country_id"],
                "level": "division",
                "actions": [
                    {
                        "className": "fa fa-info",
                        "title": str(
                            division_counts[
                                (
                                    loc["region_id"],
                                    loc["country_id"],
                                    loc["division_id"],
                                )
                            ]
                        )
                        + " sequences",
                        "text": human_format(
                            division_counts[
                                (
                                    loc["region_id"],
                                    loc["country_id"],
                                    loc["division_id"],
                                )
                            ]
                        ),
                    }
                ],
                "children": [],
            }
            country_node["children"].append(division_node)

        # Add location --> division
        if loc["location_id"] == -1:
            continue

        location_node = [
            c for c in division_node["children"] if c["value"] == loc["location_id"]
        ]
        if location_node:
            location_node = location_node[0]
        else:
            location_node = {
                "label": loc["location"],
                "value": loc["location_id"],
                "value_txt": loc["location"],
                "region": loc["region_id"],
                "country": loc["country_id"],
                "division": loc["division_id"],
                "level": "location",
                "actions": [
                    {
                        "className": "fa fa-info",
                        "title": str(
                            location_counts[
                                (
                                    loc["region_id"],
                                    loc["country_id"],
                                    loc["division_id"],
                                    loc["location_id"],
                                )
                            ]
                        )
                        + " sequences",
                        "text": human_format(
                            location_counts[
                                (
                                    loc["region_id"],
                                    loc["country_id"],
                                    loc["division_id"],
                                    loc["location_id"],
                                )
                            ]
                        ),
                    }
                ],
                "children": [],
            }
            division_node["children"].append(location_node)

    with open(args.out, "w") as fp:
        fp.write(json.dumps(select_tree))


if __name__ == "__main__":
    main()
