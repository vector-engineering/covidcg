# coding: utf-8

"""Build tree for ReactDropdownTreeSelect

Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""

import json
import pandas as pd

from scripts.util import human_format


def build_location_tree(case_data, location_map, emoji_map_file, geo_select_tree_out):
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

    df = pd.read_csv(case_data).set_index("Accession ID")
    with open(location_map, "r") as fp:
        location_map_df = pd.DataFrame(json.loads(fp.read()))

    # Join location data back to main dataframe
    df = df.join(location_map_df, on="location_id")

    # Set unspecified locations to None so that they don't get
    # caught up in the groupby
    df.loc[df["region"] == "-1", "region"] = None
    df.loc[df["country"] == "-1", "country"] = None
    df.loc[df["division"] == "-1", "division"] = None
    df.loc[df["location"] == "-1", "location"] = None

    # Count sequences per grouping level
    region_counts = dict(df.groupby("region")["location_id"].count())
    country_counts = dict(df.groupby(["region", "country"])["location_id"].count())
    division_counts = dict(
        df.groupby(["region", "country", "division"])["location_id"].count()
    )
    location_counts = dict(
        df.groupby(["region", "country", "division", "location"])["location_id"].count()
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
                "location_id": i,
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
                "location_id": i,
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
                "location_id": i,
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
                "location_id": i,
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

    with open(geo_select_tree_out, "w") as fp:
        fp.write(json.dumps(select_tree))
