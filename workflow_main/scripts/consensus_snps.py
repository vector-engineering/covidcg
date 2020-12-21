# coding: utf-8

"""Process pangolin lineage/clade information

Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""

import pandas as pd
import numpy as np

from itertools import chain
from collections import Counter


def get_consensus_snps(case_df, group_key, consensus_fraction=0.9):
    """Generalized for lineages/clades

    Parameters
    ----------
    case_df: pandas.DataFrame
    group_key: str
        - 'lineage' or 'clade'
    consensus_fraction: float
        - Fraction of taxons that need to have a SNP for it to be considered 
          a consensus SNP for a lineage/clade

    Returns
    -------
    group_snp_df: pandas.DataFrame
    """

    group_snp_df = []
    unique_groups = sorted(case_df[group_key].unique())

    collapsed_snvs = (
        case_df.groupby(group_key)[
            ["dna_snp_str", "gene_aa_snp_str", "protein_aa_snp_str"]
        ]
        # Instead of trying to flatten this list of lists
        # (and calling like 100K+ mem allocs)
        # Just make an iterator over each element of the nested lists
        # Also, package the number of isolates per lineage in with the
        # iterator so we can use a single aggregate function later
        # to determine which SNVs are consensus SNVs
        .agg(list).applymap(lambda x: (len(x), chain.from_iterable(x)))
    )

    # For a given number of isolates per group, and the iterator
    # over all SNVs in that group, get the consensus SNVs
    def count_consensus(pkg):
        # Unbundle the input tuple
        n_isolates, it = pkg
        snvs = list(it)
        # Count unique occurrences
        counts = dict(Counter(snvs))
        # Base the minimum frequency off of the required consensus fraction
        # and the number of isolates for this group
        min_freq = int(consensus_fraction * n_isolates)

        # Return all SNVs which pass the min_freq threshold, sorted
        # in ascending order?
        return sorted([int(k) for k, v in counts.items() if v >= min_freq])

    # Do this column-by-column because for some reason pandas applymap()
    # misses the first column. I have no idea why...
    for col in ["dna_snp_str", "gene_aa_snp_str", "protein_aa_snp_str"]:
        collapsed_snvs[col] = collapsed_snvs[col].apply(count_consensus)

    collapsed_snvs = (
        collapsed_snvs
        .reset_index()
        .rename(
            columns={
                "dna_snp_str": "dna_snp_ids",
                "gene_aa_snp_str": "gene_aa_snp_ids",
                "protein_aa_snp_str": "protein_aa_snp_ids",
            }
        )
    )

    return collapsed_snvs


def get_all_consensus_snps(case_data, lineage_out, clade_out, consensus_fraction=0.9):
    """For each lineage and clade, get the lineage/clade-defining SNVs,
    on both the NT and AA level
    Lineage/clade-defining SNVs are defined as SNVs which occur in
    >= [consensus_fraction] of sequences within that lineage/clade.
    [consensus_fraction] is a parameter which can be adjusted here
    """

    case_df = pd.read_csv(case_data, index_col="Accession ID")

    # Serialized list back to list 
    cols = ["dna_snp_str", "gene_aa_snp_str", "protein_aa_snp_str"]
    for col in cols:
        case_df[col] = (
            case_df[col]
            .str.strip("[]")
            .str.split(",")
            .apply(lambda x: [int(_x) for _x in x])
        )

    lineage_snp_df = get_consensus_snps(
        case_df, "lineage",
        consensus_fraction=consensus_fraction
    )
    lineage_snp_df.to_json(lineage_out, orient="records")

    clade_snp_df = get_consensus_snps(
        case_df, "clade",
        consensus_fraction=consensus_fraction
    )
    clade_snp_df.to_json(clade_out, orient="records")
