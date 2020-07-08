# coding: utf-8

"""Process clade information

Author: Albert Chen (Deverman Lab, Broad Institute)
"""

import pandas as pd
import numpy as np

from collections import Counter
from pathlib import Path

from util import translate, data_dir


def get_consensus_snps(case_df):

    # Fraction of taxons that need to have a SNP for it to be considered a consensus
    # SNP for a lineage
    # TODO: make this a CLI arg
    consensus_fraction = 0.9

    lineage_snp_df = []
    unique_lineages = sorted(case_df["lineage"].unique())

    for i, lineage in enumerate(unique_lineages):
        lin_df = case_df.loc[case_df["lineage"] == lineage, :]

        dna_snp_freqs = dict(
            Counter(
                sum(lin_df["dna_snp_str"].apply(lambda x: x.split(";")).values, [])
            ).most_common()
        )

        dna_consensus_snps = sorted(
            [
                int(k)
                for k, v in dna_snp_freqs.items()
                if (v / len(lin_df)) >= consensus_fraction
            ]
        )

        aa_snp_freqs = dict(
            Counter(
                sum(lin_df["aa_snp_str"].apply(lambda x: x.split(";")).values, [])
            ).most_common()
        )

        aa_consensus_snps = sorted(
            [
                int(k)
                for k, v in aa_snp_freqs.items()
                if (v / len(lin_df)) >= consensus_fraction
            ]
        )

        lineage_snp_df.append((lineage, dna_consensus_snps, aa_consensus_snps))

    lineage_snp_df = pd.DataFrame.from_records(
        lineage_snp_df, columns=["lineage", "dna_snp_ids", "aa_snp_ids"]
    )
    lineage_snp_df["dna_snp_ids"] = lineage_snp_df["dna_snp_ids"].apply(
        lambda x: ";".join([str(_x) for _x in x])
    )
    lineage_snp_df["aa_snp_ids"] = lineage_snp_df["aa_snp_ids"].apply(
        lambda x: ";".join([str(_x) for _x in x])
    )

    # Save to disk
    lineage_snp_df.to_csv(data_dir / "lineage_snp.csv", index=False)
    lineage_snp_df.to_json(data_dir / "lineage_snp.json", orient="records")

    return lineage_snp_df
