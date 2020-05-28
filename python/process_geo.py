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

from pathlib import Path

project_root_path = Path(__file__).resolve().parent.parent
data_dir = (project_root_path / 'data').resolve() # Resolve any symlinks --> absolute path

def load_geo_data():
    '''Load location data
    '''

    patient_meta_files = sorted((data_dir / 'patient_meta').glob('*.tsv'))
    print('Collecting {} patient metadata files...'.format(len(patient_meta_files)), end='', flush=True)
    patient_meta_df = pd.DataFrame()
    for f in patient_meta_files:
        _df = pd.read_csv(f, sep='\t', skiprows=2)
        patient_meta_df = pd.concat([patient_meta_df, _df], ignore_index=True)

    # Save dataframe
    patient_meta_df.to_csv('processed_data/patient_meta.csv', index=False)
    print('done', flush=True)

    # Location data is stored in one column, "region / country / division / location"
    location_df = (
        patient_meta_df['Location'].str.split('/', expand=True)
        .iloc[:, :4] # Only take 4 columns
        # Rename columns
        .rename(columns={0: 'region', 1: 'country', 2: 'division', 3: 'location'})
        .applymap(lambda x: x.strip() if x else x)
        # Placeholder for missing values, so that it will still 
        # be caught by groupby() later on
        .fillna(-1)
    )
    # Re-add metadata columns
    location_df['name'] = patient_meta_df['Virus name']
    location_df['gisaid_id'] = patient_meta_df['Accession ID']
    location_df['sample_date'] = patient_meta_df['Collection date']

    # Convert sample_date to datetime
    location_df['sample_date'] = pd.to_datetime(location_df['sample_date'], yearfirst=True)

    # Create complete location column from the separate parts
    # This time with no padding
    location_df['loc_str'] = location_df['region'].str.cat(
        [location_df['country'].astype(str), location_df['division'].astype(str), location_df['location'].astype(str)],
        sep='/'
    )

    unique_location_df = (
        location_df['loc_str']
        .drop_duplicates()
        .sort_values(ignore_index=True)
        .reset_index() # Produce an index column
    )
    # Location data is stored in one column, "region/country/division/location"
    unique_location_df = pd.concat([unique_location_df, (
        unique_location_df['loc_str'].str.split('/', expand=True)
        .iloc[:, :4] # Only take 4 columns
        # Rename columns
        .rename(columns={0: 'region', 1: 'country', 2: 'division', 3: 'location'})
    )], axis=1)

    # Map location IDs back to taxon_locations dataframe
    location_df['location_id'] = location_df['loc_str'].map(pd.Series(
        unique_location_df['index'].values, index=unique_location_df['loc_str'].values
    ))

    # Take subset of columns, re-index
    location_df = (
        location_df
        [['name', 'gisaid_id', 'sample_date', 'location_id']]
        .sort_values('location_id')
        .reset_index(drop=True)
    )

    # print(location_df)

    print('Saving taxon locations')
    location_df.to_csv(data_dir / 'taxon_locations.csv', index=False)
    # location_df.to_json(data_dir / 'taxon_locations.json', orient='records')

    print('Saving unique locations')
    unique_location_df.drop(columns=['loc_str']).to_csv(data_dir / 'location_map.csv', index=False)
    unique_location_df.drop(columns=['loc_str']).to_json(data_dir / 'location_map.json', orient='records')

    return location_df, unique_location_df


def build_select_tree(unique_location_df):
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

    for i, loc in unique_location_df.iterrows():
        # Add region node
        if loc['region'] == '-1':
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
        if loc['country'] == '-1':
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
        if loc['division'] == '-1':
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
        if loc['location'] == '-1':
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
    print('Saving geo select tree')
    select_tree_path = data_dir / 'geo_select_tree.json'
    with select_tree_path.open('w') as fp:
        fp.write(json.dumps(select_tree))

    # print(loc_tree.nodes)
    # print(loc_tree.in_edges('New York'))
    return select_tree


def main():
    location_df, unique_location_df = load_geo_data()
    select_tree = build_select_tree(unique_location_df)
    # print(select_tree)


if __name__ == '__main__':
    main()

