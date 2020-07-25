# coding: utf-8

"""Process pangolin lineage/clade information

Author: Albert Chen (Deverman Lab, Broad Institute)
"""

import pandas as pd
import numpy as np

from collections import Counter

from util import translate, data_dir


def get_all_consensus_snps(case_df):
    get_consensus_snps(case_df, "lineage")
    get_consensus_snps(case_df, "clade")


def get_consensus_snps(case_df, group_key):
    """Generalized for lineages/clades

    Parameters
    ----------
    case_df: pandas.DataFrame
    group_key: str
        - 'lineage' or 'clade'

    Returns
    -------
    group_snp_df: pandas.DataFrame
    """

    # Fraction of taxons that need to have a SNP for it to be considered a consensus
    # SNP for a lineage/clade
    # TODO: make this a CLI arg
    consensus_fraction = 0.9

    group_snp_df = []
    unique_groups = sorted(case_df[group_key].unique())

    for i, group in enumerate(unique_groups):
        lin_df = case_df.loc[case_df[group_key] == group, :]

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

        group_snp_df.append(
            (
                group,
                dna_consensus_snps,
                gene_aa_consensus_snps,
                protein_aa_consensus_snps,
            )
        )

    group_snp_df = pd.DataFrame.from_records(
        group_snp_df,
        columns=[group_key, "dna_snp_ids", "gene_aa_snp_ids", "protein_aa_snp_ids"],
    )

    # Save to disk
    group_snp_df.to_csv(data_dir / (group_key + "_snp.csv"), index=False)
    group_snp_df.to_json(data_dir / (group_key + "_snp.json"), orient="records")

    return group_snp_df


if __name__ == "__main__":
    case_df = pd.read_csv(data_dir / "case_data2.csv")
    case_df.set_index("Accession ID")

    case_df["dna_snp_str"] = (
        case_df["dna_snp_str"]
        .str.strip("[]")
        .str.split(",")
        .apply(lambda x: [int(_x) for _x in x])
    )
    case_df["gene_aa_snp_str"] = (
        case_df["gene_aa_snp_str"]
        .str.strip("[]")
        .str.split(",")
        .apply(lambda x: [int(_x) for _x in x])
    )
    case_df["protein_aa_snp_str"] = (
        case_df["protein_aa_snp_str"]
        .str.strip("[]")
        .str.split(",")
        .apply(lambda x: [int(_x) for _x in x])
    )

    get_all_consensus_snps(case_df)
