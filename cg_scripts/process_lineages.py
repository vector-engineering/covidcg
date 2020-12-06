# coding: utf-8

"""Process pangolin lineage/clade information

Author: Albert Chen (Deverman Lab, Broad Institute)
"""

import itertools
import pandas as pd
import numpy as np

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
        .agg(list).applymap(lambda x: (len(x), itertools.chain.from_iterable(x)))
    )

    # For a given number of isolates per group, and the iterator
    # over all SNVs in that group, get the consensus SNVs
    def count_consensus(pkg):
        # Unbundle the input tuple
        n_isolates, it = pkg
        # Count unique occurrences
        counts = dict(Counter(it))
        # Base the minimum frequency off of the required consensus fraction
        # and the number of isolates for this group
        min_freq = int(consensus_fraction * n_isolates)
        # Return all SNVs which pass the min_freq threshold, sorted
        # in ascending order?
        return sorted([int(k) for k, v in counts.items() if v >= min_freq])

    consensus_snvs = (
        collapsed_snvs.applymap(count_consensus)
        .reset_index()
        .rename(
            columns={
                "dna_snp_str": "dna_snp_ids",
                "gene_aa_snp_str": "gene_aa_snp_ids",
                "protein_aa_snp_str": "protein_aa_snp_ids",
            }
        )
    )

    return consensus_snvs
