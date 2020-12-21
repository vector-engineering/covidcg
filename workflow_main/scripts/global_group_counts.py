# coding: utf-8

"""Count sequences per group

Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""

import json
import pandas as pd

from itertools import chain
from collections import Counter


def global_group_counts(case_data, out_global_group_counts):
    """Get the number of sequences in each group
    Doing this in the pipeline just saves some work for the browser later
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

    global_group_counts = {}
    # Count lineages and clades
    global_group_counts["lineage"] = case_df["lineage"].value_counts().to_dict()
    global_group_counts["clade"] = case_df["clade"].value_counts().to_dict()

    # Count global SNV frequencies
    # Collapse list of lists into one list, then count individual
    # occurrences, then cast to a regular dict
    global_group_counts["dna_snp"] = dict(
        Counter(list(chain.from_iterable(case_df["dna_snp_str"])))
    )
    global_group_counts["gene_aa_snp"] = dict(
        Counter(list(chain.from_iterable(case_df["gene_aa_snp_str"])))
    )
    global_group_counts["protein_aa_snp"] = dict(
        Counter(list(chain.from_iterable(case_df["protein_aa_snp_str"])))
    )

    with open(out_global_group_counts, "w") as fp:
        fp.write(json.dumps(global_group_counts))
