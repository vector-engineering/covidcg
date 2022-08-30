#!/usr/bin/env python3
# coding: utf-8

"""Process pangolin lineage/clade information

Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""

import argparse
import json
import pandas as pd

from itertools import chain
from collections import Counter


def get_consensus_mutations(
    isolate_df,
    group_key,
    reference,
    consensus_fraction=0.9,
    min_reporting_fraction=0.05,
):
    """Generalized for lineages/clades

    Parameters
    ----------
    isolate_df: pandas.DataFrame
    group_key: str
        - 'lineage' or 'clade'
    reference: str
        - Reference name
    consensus_fraction: float
        - Fraction of taxons that need to have a mutation for it to be considered 
          a consensus mutation for a lineage/clade
    min_reporting_fraction: float
        - Mutations below this fractional frequency will be ignored

    Returns
    -------
    group_mutation_df: pandas.DataFrame
    """

    isolate_dff = isolate_df.loc[isolate_df["reference"] == reference]

    collapsed_mutations = (
        isolate_dff.groupby(group_key)[
            ["dna_mutation", "gene_aa_mutation", "protein_aa_mutation"]
        ]
        # Instead of trying to flatten this list of lists
        # (and calling like 100K+ mem allocs)
        # Just make an iterator over each element of the nested lists
        # Also, package the number of isolates per lineage in with the
        # iterator so we can use a single aggregate function later
        # to determine which mutations are consensus mutations
        .agg(list).applymap(lambda x: (len(x), chain.from_iterable(x)))
    )

    # For a given number of isolates per group, and the iterator
    # over all mutations in that group, get the consensus mutations
    def count_consensus(pkg):
        # Unbundle the input tuple
        n_isolates, it = pkg
        mutations = list(it)
        # Count unique occurrences
        counts = dict(Counter(mutations))
        # Base the minimum frequency off of the required consensus fraction
        # and the number of genomes for this group
        min_freq = int(min_reporting_fraction * n_isolates)

        # Return all mutations which pass the min_freq threshold, sorted
        # in ascending order by mutation ID
        # Also include the counts and the fraction of counts relative to
        # the number of genomes for this group
        return sorted(
            [(int(k), v, v / n_isolates) for k, v in counts.items() if v >= min_freq]
        )

    # Do this column-by-column because for some reason pandas applymap()
    # misses the first column. I have no idea why...
    for col in ["dna_mutation", "gene_aa_mutation", "protein_aa_mutation"]:
        collapsed_mutations[col] = collapsed_mutations[col].apply(count_consensus)

    def mutations_to_df(field):
        _df = collapsed_mutations[field].explode()
        _df = _df.loc[~pd.isna(_df)]
        _df = pd.DataFrame.from_records(
            _df, columns=["mutation_id", "count", "fraction"], index=_df.index
        )
        return _df

    collapsed_dna_mutations = mutations_to_df("dna_mutation")
    collapsed_gene_aa_mutations = mutations_to_df("gene_aa_mutation")
    collapsed_protein_aa_mutations = mutations_to_df("protein_aa_mutation")

    return (
        # CONSENSUS MUTATIONS
        (
            collapsed_mutations.applymap(
                lambda x: [mut[0] for mut in x if mut[2] > consensus_fraction]
            ).rename(
                columns={
                    "dna_mutation": "dna_mutation_ids",
                    "gene_aa_mutation": "gene_aa_mutation_ids",
                    "protein_aa_mutation": "protein_aa_mutation_ids",
                }
            )
        ).to_dict(orient="index"),
        # MUTATION FREQUENCIES PER GROUP
        {
            "dna": (
                collapsed_dna_mutations.reset_index()
                .rename(columns={group_key: "group"})
                .to_dict(orient="records")
            ),
            "gene_aa": (
                collapsed_gene_aa_mutations.reset_index()
                .rename(columns={group_key: "group"})
                .to_dict(orient="records")
            ),
            "protein_aa": (
                collapsed_protein_aa_mutations.reset_index()
                .rename(columns={group_key: "group"})
                .to_dict(orient="records")
            ),
        },
    )


def main():
    """For each lineage and clade, get the lineage/clade-defining mutations,
    on both the NT and AA level
    Lineage/clade-defining mutations are defined as mutations which occur in
    >= [consensus_fraction] of sequences within that lineage/clade.
    [consensus_fraction] is a parameter which can be adjusted here
    """

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--isolate-data", type=str, required=True, help="Path to isolate data JSON file"
    )
    parser.add_argument(
        "--reference", type=str, required=True, help="Path to Reference JSON file"
    )
    parser.add_argument(
        "--consensus-out", type=str, required=True, help="Consensus mutations JSON file"
    )
    parser.add_argument(
        "--frequencies-out", type=str, required=True, help="Frequencies JSON file"
    )
    parser.add_argument(
        "--group-cols", type=str, nargs="+", default=[], help="Group columns"
    )
    parser.add_argument(
        "--consensus-fraction", type=float, default=0.9, help="Consensus fraction"
    ),
    parser.add_argument(
        "--min-reporting-fraction",
        type=float,
        default=0.05,
        help="Min reporting fraction",
    )

    args = parser.parse_args()

    # If no group columns are defined, then just write an empty JSON object
    if not args.group_cols:
        with open(args.consensus_out, "w") as fp:
            fp.write(json.dumps({}))

    isolate_df = pd.read_json(args.isolate_data).set_index("isolate_id")[
        args.group_cols
        + ["reference", "dna_mutation", "gene_aa_mutation", "protein_aa_mutation"]
    ]

    # Get references
    with open(args.reference, "rt") as fp:
        references = json.load(fp)

    reference_names = sorted(list(references.keys()))

    consensus_dict = {}
    frequencies_dict = {}
    # group_cols is defined under the "group_cols" field
    # in the config.yaml file
    for group in args.group_cols:

        consensus_dict[group] = {}
        frequencies_dict[group] = {}

        for reference_name in reference_names:

            group_consensus, group_frequencies = get_consensus_mutations(
                isolate_df,
                group,
                reference=reference_name,
                consensus_fraction=args.consensus_fraction,
                min_reporting_fraction=args.min_reporting_fraction,
            )
            consensus_dict[group][reference_name] = group_consensus
            frequencies_dict[group][reference_name] = group_frequencies

    with open(args.consensus_out, "w") as fp:
        fp.write(json.dumps(consensus_dict))

    with open(args.frequencies_out, "w") as fp:
        fp.write(json.dumps(frequencies_dict))


if __name__ == "__main__":
    main()
