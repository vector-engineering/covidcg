# coding: utf-8

"""Process pangolin lineage/clade information

Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""

import json
import numpy as np
import pandas as pd

from itertools import chain
from collections import Counter


def get_consensus_mutations(
    case_df, group_key, consensus_fraction=0.9, min_reporting_fraction=0.05
):
    """Generalized for lineages/clades

    Parameters
    ----------
    case_df: pandas.DataFrame
    group_key: str
        - 'lineage' or 'clade'
    consensus_fraction: float
        - Fraction of taxons that need to have a SNP for it to be considered 
          a consensus SNP for a lineage/clade
    min_reporting_fraction: float
        - ...

    Returns
    -------
    group_snp_df: pandas.DataFrame
    """

    group_snp_df = []
    unique_groups = sorted(case_df[group_key].unique())

    collapsed_mutations = (
        case_df.groupby(group_key)[
            ["dna_snp_str", "gene_aa_snp_str", "protein_aa_snp_str"]
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
    for col in ["dna_snp_str", "gene_aa_snp_str", "protein_aa_snp_str"]:
        collapsed_mutations[col] = collapsed_mutations[col].apply(count_consensus)

    def mutations_to_df(field):
        _df = collapsed_mutations[field].explode()
        _df = pd.DataFrame.from_records(
            _df, columns=["mutation_id", "count", "fraction"], index=_df.index
        )
        return _df

    collapsed_dna_mutations = mutations_to_df("dna_snp_str")
    collapsed_gene_aa_mutations = mutations_to_df("gene_aa_snp_str")
    collapsed_protein_aa_mutations = mutations_to_df("protein_aa_snp_str")

    return (
        # CONSENSUS MUTATIONS
        (
            collapsed_mutations.applymap(
                lambda x: [mut[0] for mut in x if mut[2] > consensus_fraction]
            ).rename(
                columns={
                    "dna_snp_str": "dna_snp_ids",
                    "gene_aa_snp_str": "gene_aa_snp_ids",
                    "protein_aa_snp_str": "protein_aa_snp_ids",
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


def get_all_consensus_mutations(
    case_data,
    consensus_out,
    frequencies_out,
    group_cols=[],
    consensus_fraction=0.9,
    min_reporting_fraction=0.05,
):
    """For each lineage and clade, get the lineage/clade-defining mutations,
    on both the NT and AA level
    Lineage/clade-defining mutations are defined as mutations which occur in
    >= [consensus_fraction] of sequences within that lineage/clade.
    [consensus_fraction] is a parameter which can be adjusted here
    """

    # If no group columns are defined, then just write an empty JSON object
    if not group_cols:
        with open(consensus_out, "w") as fp:
            fp.write(json.dumps({}))

    case_df = pd.read_json(case_data).set_index("Accession ID")

    consensus_dict = {}
    frequencies_dict = {}
    # group_cols is defined under the "group_cols" field
    # in the config.yaml file
    for group in group_cols:
        group_consensus, group_frequencies = get_consensus_mutations(
            case_df,
            group,
            consensus_fraction=consensus_fraction,
            min_reporting_fraction=min_reporting_fraction,
        )
        consensus_dict[group] = group_consensus
        frequencies_dict[group] = group_frequencies

    with open(consensus_out, "w") as fp:
        fp.write(json.dumps(consensus_dict))

    with open(frequencies_out, "w") as fp:
        fp.write(json.dumps(frequencies_dict))
