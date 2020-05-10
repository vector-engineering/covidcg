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
        .reset_index()
    )
    # print(location_df)

    return location_df


def build_geo_graph(location_df):
    '''Build Directed Graph of Geographies

    Tree structure:
    Region <-- Country <-- Division <-- Location
    '''

    loc_tree = nx.DiGraph()
    loc_tree.add_node('All') # Root node
    for i, loc in location_df.iterrows():
        # Add region node
        if loc['region'] == -1:
            continue

        loc_tree.add_node(loc['region'])
        loc_tree.add_edge(loc['region'], 'All')

        # Add country --> region
        if loc['country'] == -1:
            continue

        loc_tree.add_node(loc['country'])
        loc_tree.add_edge(loc['country'], loc['region'])
        
        # Add division --> country
        if loc['division'] == -1:
            continue

        loc_tree.add_node(loc['division'])
        loc_tree.add_edge(loc['division'], loc['country'])

        # Add location --> division
        if loc['location'] == -1:
            continue
        
        loc_tree.add_node(loc['location'])
        loc_tree.add_edge(loc['location'], loc['division'])

    # print(loc_tree.nodes)
    # print(loc_tree.in_edges('New York'))
    return loc_tree


def build_select_tree(loc_tree):
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

    # Transform names into more programmatic values
    def geo_label_to_value(label):
        # Lowercase
        label = label.lower()
        # Replace spaces with '_'
        label = re.sub(r'\s', '_', label)
        return label
    # print(geo_label_to_value('New York City'))


    # Traverse tree with recursion
    def build_select_tree_recurse(node):
        # Build self
        node_obj = {
            'label': node,
            'value': node
        }

        # Get all child nodes
        child_nodes = [n[0] for n in loc_tree.in_edges(node) if n[0] != node]
        child_objs = []
        if child_nodes:
            node_obj['children'] = [
                build_select_tree_recurse(child_node)
                for child_node in child_nodes
            ]
            
        
        return node_obj

    select_tree = build_select_tree_recurse('All')
    # print(json.dumps(select_tree, indent=2))

    # Save tree as json file
    with open('processed_data/geography.json', 'w') as fp:
        fp.write(json.dumps(select_tree))

    return select_tree


def main():
    location_df = load_geo_data()
    loc_tree = build_geo_graph(location_df)

    build_select_tree(loc_tree)


if __name__ == '__main__':
    main()

