"""
Script to create the dominant sequences data from the raw data for mobile use.
"""

import pandas as pd
from pathlib import Path
import argparse
import json
import os

def dominant_sequences():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-o", "--output", type=str, required=True, help="Path to output directory",
    )
    parser.add_argument(
        "--case-data", type=str, required=True, help="Path to case data JSON file",
    )
    parser.add_argument("--metadata-map", type=str, required=True, 
        help="Path to metadata mapping JSON file.")
    args = parser.parse_args()
    out_path = Path(args.output)
    with open(args.metadata_map, "r") as fp:
            metadata_map = json.loads(fp.read())
    case_df = pd.read_json(args.case_data).set_index('Accession ID')
    # Join locations onto case_data

    loc_levels = ["region", "country", "division", "location"]
    for loc_level in loc_levels:
        case_df.loc[:, loc_level] = case_df[loc_level].map(
            {int(k): v for k, v in metadata_map[loc_level].items()}
        )
        case_df.loc[case_df[loc_level].isna(), loc_level] = None
    case_df = case_df[['country', 'collection_date', 'lineage']]
    # Iterate over columns and apply mapping.
    for c in case_df:
        if c in metadata_map:
            case_df[c] = case_df[c].astype(str).replace(metadata_map[c])
    case_df = case_df.groupby(['country', 'lineage']).count()
    
    # Rename count column.
    case_df = case_df.rename({'collection_date': 'count'}, axis='columns')
    # Reset index so that country is available for every single row.
    case_df = case_df.reset_index()
    case_df.to_json(str(out_path / "dominant_sequences.json"), orient='records')

    # Export data
    #case_df.to_json(out_path + 'dominant_sequences.json')
    

if __name__ == '__main__':
    dominant_sequences()