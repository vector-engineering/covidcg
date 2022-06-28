#!/usr/bin/env python3
# coding: utf-8

"""Count sequences per group, plus mutation frequencies per group

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

    parser.add_argument(
        "--isolate-data", type=str, required=True, help="Isolate data JSON file"
    )
    parser.add_argument(
        "--reference",
        type=str,
        required=True,
        help="Path to reference sequences fasta file",
    )
    parser.add_argument(
        "--out-global-group-counts",
        type=str,
        required=True,
        help="Global group counts output file",
    )
    parser.add_argument(
        "--group-cols", type=str, nargs="+", default=[], help="Group columns"
    )

    args = parser.parse_args()

    isolate_df = pd.read_json(args.isolate_data)

    # Get reference sequence names
    with open(args.reference, "rt") as fp:
        references = json.load(fp)

    subtypes = list(sorted(references.keys()))
    reference_names = sorted(
        sum([list(references[subtype].keys()) for subtype in subtypes], [])
    )

    global_group_counts = {}
    for ref_name in reference_names:
        global_group_counts[ref_name] = {}

    # Count sequence groupings (e.g., lineages and clades)
    for ref_name in reference_names:
        for col in args.group_cols:
            global_group_counts[ref_name][col] = (
                isolate_df.drop_duplicates("isolate_id")[col].value_counts().to_dict()
            )

    # Count global mutation frequencies
    # Collapse list of lists into one list, then count individual
    # occurrences, then cast to a regular dict
    for ref_name in reference_names:
        global_group_counts[ref_name]["dna_mutation"] = dict(
            Counter(
                list(
                    chain.from_iterable(
                        isolate_df.loc[
                            isolate_df["reference"] == ref_name, "dna_mutation"
                        ]
                    )
                )
            )
        )
        global_group_counts[ref_name]["gene_aa_mutation"] = dict(
            Counter(
                list(
                    chain.from_iterable(
                        isolate_df.loc[
                            isolate_df["reference"] == ref_name, "gene_aa_mutation"
                        ]
                    )
                )
            )
        )
        global_group_counts[ref_name]["protein_aa_mutation"] = dict(
            Counter(
                list(
                    chain.from_iterable(
                        isolate_df.loc[
                            isolate_df["reference"] == ref_name, "protein_aa_mutation",
                        ]
                    )
                )
            )
        )

    with open(args.out_global_group_counts, "w") as fp:
        fp.write(json.dumps(global_group_counts))


if __name__ == "__main__":
    main()
