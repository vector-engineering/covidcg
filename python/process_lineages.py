# coding: utf-8

"""Process clade information

Author: Albert Chen (Deverman Lab, Broad Institute)
"""

import pandas as pd
import numpy as np

from collections import Counter

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
            Counter(sum(lin_df["dna_snp_str"].values, [])).most_common()
        )

        dna_consensus_snps = sorted(
            [
                int(k)
                for k, v in dna_snp_freqs.items()
                if (v / len(lin_df)) >= consensus_fraction
            ]
        )

        gene_aa_snp_freqs = dict(
            Counter(sum(lin_df["gene_aa_snp_str"].values, [])).most_common()
        )

        gene_aa_consensus_snps = sorted(
            [
                int(k)
                for k, v in gene_aa_snp_freqs.items()
                if (v / len(lin_df)) >= consensus_fraction
            ]
        )

        protein_aa_snp_freqs = dict(
            Counter(sum(lin_df["protein_aa_snp_str"].values, [])).most_common()
        )

        protein_aa_consensus_snps = sorted(
            [
                int(k)
                for k, v in protein_aa_snp_freqs.items()
                if (v / len(lin_df)) >= consensus_fraction
            ]
        )

        lineage_snp_df.append(
            (
                lineage,
                dna_consensus_snps,
                gene_aa_consensus_snps,
                protein_aa_consensus_snps,
            )
        )

    lineage_snp_df = pd.DataFrame.from_records(
        lineage_snp_df,
        columns=["lineage", "dna_snp_ids", "gene_aa_snp_ids", "protein_aa_snp_ids"],
    )

    # Save to disk
    lineage_snp_df.to_csv(data_dir / "lineage_snp.csv", index=False)
    lineage_snp_df.to_json(data_dir / "lineage_snp.json", orient="records")

    return lineage_snp_df
