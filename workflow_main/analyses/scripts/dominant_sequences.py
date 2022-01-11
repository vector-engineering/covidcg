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
    parser.add_argument("--metadata", type=str, required=True, 
        help="Path to metadata mapping JSON file.")
    args = parser.parse_args()
    out_path = Path(args.output)
    with open(args.metadata, "r") as fp:
            metadata_map = json.loads(fp.read())
    case_df = pd.read_json(args.case_data).set_index('Accession ID')
    case_df = case_df[['country', 'collection_date', 'lineage']]
    # Iterate over columns and apply mapping.
    for c in case_df:
        if c in metadata_map:
            case_df[c] = case_df[c].astype(str).replace(metadata_map[c])
    case_df = case_df.groupby(['country', 'lineage']).count()
    case_df.to_csv(str(out_path / "sequences_per_month.json"))

    # Export data
    #case_df.to_json(out_path + 'dominant_sequences.json')
    

if __name__ == '__main__':
    dominant_sequences()