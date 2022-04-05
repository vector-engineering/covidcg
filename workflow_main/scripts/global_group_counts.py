#!/usr/bin/env python3
# coding: utf-8

"""Count sequences per group

Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""

import argparse
import json
import pandas as pd

from itertools import chain
from collections import Counter


def main():
    """Get the number of sequences in each group
    Doing this in the pipeline just saves some work for the browser later
    """

    parser = argparse.ArgumentParser()

    parser.add_argument('--case-data', type=str, required=True, help='Case data JSON file')
    parser.add_argument('--out-global-group-counts', type=str, required=True, help='Global group counts output file')
    parser.add_argument('--group-cols', type=str, nargs='+', default=[], help='Group columns')

    args = parser.parse_args()

    case_df = pd.read_json(args.case_data).set_index("Accession ID")

    global_group_counts = {}
    # Count sequence groupings (e.g., lineages and clades)
    for col in args.group_cols:
        global_group_counts[col] = case_df[col].value_counts().to_dict()

    # Count global mutation frequencies
    # Collapse list of lists into one list, then count individual
    # occurrences, then cast to a regular dict
    global_group_counts["dna_mutation"] = dict(
        Counter(list(chain.from_iterable(case_df["dna_mutation_str"])))
    )
    global_group_counts["gene_aa_mutation"] = dict(
        Counter(list(chain.from_iterable(case_df["gene_aa_mutation_str"])))
    )
    global_group_counts["protein_aa_mutation"] = dict(
        Counter(list(chain.from_iterable(case_df["protein_aa_mutation_str"])))
    )

    with open(args.out_global_group_counts, "w") as fp:
        fp.write(json.dumps(global_group_counts))

if __name__ == '__main__':
    main()