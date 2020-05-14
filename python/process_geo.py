#!/usr/bin/env python3
# coding: utf-8

'''Process clade information
Author: Albert Chen (Deverman Lab, Broad Institute)
'''

import json
import numpy as np
import networkx as nx
import pandas as pd
import re


def load_ns_metadata():
    '''Load NextStrain metadata
    '''

    ns_meta_df = pd.read_csv('nextstrain_metadata.tsv', sep='\t')
    # print(ns_meta_df)
    # print(ns_meta_df.columns)
    return ns_meta_df


def load_geo_data():
    '''Load location data
    '''

    ns_meta_df = load_ns_metadata()

    loc_keys = ['region', 'country', 'division', 'location']
    location_df = (
        ns_meta_df
        .loc[:, loc_keys]
        # Placeholder for missing values, so that it will still 
        # be caught by groupby()
        .fillna(-1) 
        .groupby(loc_keys)
        .count()
        .reset_index() # Unset the groupby keys
        .reset_index() # Produce an index column
    )
    # print(location_df)

    return location_df


def build_select_tree(location_df):
    '''Build tree for ReactDropdownTreeSelect

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
    '''

    # Root node
    select_tree = {
        'label': 'All',
        'value': 'All',
        'children': []
    }

    for i, loc in location_df.iterrows():
        # Add region node
        if loc['region'] == -1:
            continue
        
        region_node = [c for c in select_tree['children'] if c['value'] == loc['region']]
        if region_node:
            region_node = region_node[0]
        else:
            region_node = {
                'label': loc['region'],
                'value': loc['region'],
                'level': 'region',
                'children': []
            }
            select_tree['children'].append(region_node)

        # Add country --> region
        if loc['country'] == -1:
            continue

        country_node = [c for c in region_node['children'] if c['value'] == loc['country']]
        if country_node:
            country_node = country_node[0]
        else:
            country_node = {
                'label': loc['country'],
                'value': loc['country'],
                'region': loc['region'],
                'level': 'country',
                'children': []
            }
            region_node['children'].append(country_node)
        
        # Add division --> country
        if loc['division'] == -1:
            continue

        division_node = [c for c in country_node['children'] if c['value'] == loc['division']]
        if division_node:
            division_node = division_node[0]
        else:
            division_node = {
                'label': loc['division'],
                'value': loc['division'],
                'region': loc['region'],
                'country': loc['country'],
                'level': 'division',
                'children': []
            }
            country_node['children'].append(division_node)

        # Add location --> division
        if loc['location'] == -1:
            continue
        
        location_node = [c for c in division_node['children'] if c['value'] == loc['location']]
        if location_node:
            location_node = location_node[0]
        else:
            location_node = {
                'label': loc['location'],
                'value': loc['location'],
                'region': loc['region'],
                'country': loc['country'],
                'division': loc['division'],
                'level': 'location',
                'children': []
            }
            division_node['children'].append(location_node)

    # Save tree as json file
    with open('processed_data/geo_select_tree.json', 'w') as fp:
        fp.write(json.dumps(select_tree))

    # print(loc_tree.nodes)
    # print(loc_tree.in_edges('New York'))
    return select_tree



def main():
    location_df = load_geo_data()

    # Save location_df
    location_df.to_csv('processed_data/locations.csv', index=False)
    location_df.to_json('processed_data/locations.json', orient='records')

    select_tree = build_select_tree(location_df)
    # print(select_tree)


if __name__ == '__main__':
    main()

