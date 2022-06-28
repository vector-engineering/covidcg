#!/usr/bin/env python3
# coding: utf-8

"""Table building for AZ

Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""

import argparse
import pandas as pd
import numpy as np
import json

from collections import Counter
from itertools import chain
from pathlib import Path


def dna_mutation_to_name(mut_str):
    split = mut_str.split("|")
    # REF POS ALT
    return split[2] + split[1] + split[3]


def aa_mutation_to_name(mut_str):
    split = mut_str.split("|")
    # REF POS ALT
    return split[2] + split[1] + split[3]


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--isolate-data", type=str, required=True, help="Isolate data JSON file"
    )
    parser.add_argument(
        "--metadata-map", type=str, required=True, help="Metadata map JSON file"
    )
    parser.add_argument(
        "--report-gene",
        type=str,
        required=True,
        help="Gene name for gene-focused report files",
    )
    parser.add_argument(
        "--report-group-col",
        type=str,
        required=True,
        help="Grouping column for the grouping-focused report files",
    )
    parser.add_argument(
        "--report-out", type=str, required=True, help="Output directory"
    )
    args = parser.parse_args()

    report_out = Path(args.report_out)

    # ------------------------
    #        LOAD DATA
    # ------------------------

    isolate_df = pd.read_json(args.isolate_data)

    # Load DNA mutation ID map
    with open(args.metadata_map, "r") as fp:
        metadata_map = json.loads(fp.read())

    id_to_dna_mutation = {v: k for k, v in metadata_map["dna_mutation"].items()}
    id_to_gene_aa_mutation = {v: k for k, v in metadata_map["gene_aa_mutation"].items()}
    id_to_protein_aa_mutation = {
        v: k for k, v in metadata_map["protein_aa_mutation"].items()
    }

    # Join locations onto isolate_df
    loc_levels = ["region", "country", "division", "location"]
    for loc_level in loc_levels:
        isolate_df.loc[:, loc_level] = isolate_df[loc_level].map(
            {int(k): v for k, v in metadata_map[loc_level].items()}
        )
        isolate_df.loc[isolate_df[loc_level].isna(), loc_level] = None

    # -----------------------------
    #     MUTATION DEFINITIONS
    # -----------------------------

    dna_mutation_def = (
        pd.Series(id_to_dna_mutation)
        .rename("mutation_str")
        .rename_axis("mutation_id")
        .to_frame()
    )
    dna_mutation_def.insert(
        0, "segment", dna_mutation_def["mutation_str"].apply(lambda x: x.split("|")[0])
    )
    dna_mutation_def.insert(
        1, "pos", dna_mutation_def["mutation_str"].apply(lambda x: int(x.split("|")[1]))
    )
    dna_mutation_def.insert(
        2, "ref", dna_mutation_def["mutation_str"].apply(lambda x: x.split("|")[2])
    )
    dna_mutation_def.insert(
        3, "alt", dna_mutation_def["mutation_str"].apply(lambda x: x.split("|")[3])
    )
    dna_mutation_def.insert(
        0,
        "mutation_name",
        dna_mutation_def["mutation_str"].apply(dna_mutation_to_name),
    )
    dna_mutation_def.drop(columns=["mutation_str"], inplace=True)

    gene_aa_mutation_def = (
        pd.Series(id_to_gene_aa_mutation)
        .rename("mutation_str")
        .rename_axis("mutation_id")
        .to_frame()
    )
    gene_aa_mutation_def.insert(
        1, "gene", gene_aa_mutation_def["mutation_str"].apply(lambda x: x.split("|")[0])
    )
    gene_aa_mutation_def.insert(
        2,
        "pos",
        gene_aa_mutation_def["mutation_str"].apply(lambda x: int(x.split("|")[1])),
    )
    gene_aa_mutation_def.insert(
        3, "ref", gene_aa_mutation_def["mutation_str"].apply(lambda x: x.split("|")[2])
    )
    gene_aa_mutation_def.insert(
        4, "alt", gene_aa_mutation_def["mutation_str"].apply(lambda x: x.split("|")[3])
    )
    gene_aa_mutation_def.insert(
        0,
        "mutation_name",
        gene_aa_mutation_def["mutation_str"].apply(aa_mutation_to_name),
    )
    gene_aa_mutation_def.drop(columns=["mutation_str"], inplace=True)

    protein_aa_mutation_def = (
        pd.Series(id_to_protein_aa_mutation)
        .rename("mutation_str")
        .rename_axis("mutation_id")
        .to_frame()
    )
    protein_aa_mutation_def.insert(
        1,
        "protein",
        protein_aa_mutation_def["mutation_str"].apply(lambda x: x.split("|")[0]),
    )
    protein_aa_mutation_def.insert(
        2,
        "pos",
        protein_aa_mutation_def["mutation_str"].apply(lambda x: int(x.split("|")[1])),
    )
    protein_aa_mutation_def.insert(
        3,
        "ref",
        protein_aa_mutation_def["mutation_str"].apply(lambda x: x.split("|")[2]),
    )
    protein_aa_mutation_def.insert(
        4,
        "alt",
        protein_aa_mutation_def["mutation_str"].apply(lambda x: x.split("|")[3]),
    )
    protein_aa_mutation_def.insert(
        0,
        "mutation_name",
        protein_aa_mutation_def["mutation_str"].apply(aa_mutation_to_name),
    )
    protein_aa_mutation_def.drop(columns=["mutation_str"], inplace=True)

    # -----------------------------------
    #    GLOBAL REPORT GENE MUTATIONS
    # -----------------------------------

    collapsed_mutations = (
        isolate_df.groupby(lambda x: True)[
            ["dna_mutation", "gene_aa_mutation", "protein_aa_mutation"]
        ]
        # Instead of trying to flatten this list of lists
        # (and calling like 100K+ mem allocs)
        # Just make an iterator over each element of the nested lists
        # Also, package the number of isolates per lineage in with the
        # iterator so we can use a single aggregate function later
        # to determine which mutations are consensus mutations
        .agg(list).applymap(lambda x: chain.from_iterable(x))
    )
    dna_mutations = Counter(collapsed_mutations.iat[0, 0])
    gene_aa_mutations = Counter(collapsed_mutations.iat[0, 1])
    protein_aa_mutations = Counter(collapsed_mutations.iat[0, 2])

    # DNA
    dna_mutation_df = dna_mutation_def.join(
        pd.Series(dict(dna_mutations))
        .rename("count")
        .rename_axis("dna_mutation_id")
        .to_frame(),
        how="inner",
    )
    dna_mutation_df = dna_mutation_df.assign(
        freq=(dna_mutation_df["count"] / len(isolate_df)) * 100
    )
    dna_mutation_df = dna_mutation_df.sort_values("count", ascending=False)
    dna_mutation_df.to_csv(report_out / "dna_mutation_global.csv", index=False)

    # GENE AA
    gene_aa_mutation_df = gene_aa_mutation_def.join(
        pd.Series(dict(gene_aa_mutations))
        .rename("count")
        .rename_axis("gene_aa_mutation_id")
        .to_frame(),
        how="inner",
    )
    gene_aa_mutation_df = gene_aa_mutation_df.assign(
        freq=(gene_aa_mutation_df["count"] / len(isolate_df)) * 100
    )
    gene_aa_mutation_df = gene_aa_mutation_df.sort_values("count", ascending=False)
    gene_aa_mutation_df.to_csv(report_out / "gene_aa_mutation_global.csv", index=False)

    (
        gene_aa_mutation_df.loc[gene_aa_mutation_df["gene"] == args.report_gene]
        .drop(columns=["gene"])
        .to_csv(report_out / f"{args.report_gene}_mutation_global.csv", index=False)
    )

    # PROTEIN AA
    protein_aa_mutation_df = protein_aa_mutation_def.join(
        pd.Series(dict(protein_aa_mutations))
        .rename("count")
        .rename_axis("protein_aa_mutation_id")
        .to_frame(),
        how="inner",
    )
    protein_aa_mutation_df = protein_aa_mutation_df.assign(
        freq=(protein_aa_mutation_df["count"] / len(isolate_df)) * 100
    )
    protein_aa_mutation_df = protein_aa_mutation_df.sort_values(
        "count", ascending=False
    )
    protein_aa_mutation_df.to_csv(
        report_out / "protein_aa_mutation_global.csv", index=False
    )

    # ----------------------------------
    #   REGIONAL REPORT GENE MUTATIONS
    # ----------------------------------

    collapsed_region_mutations = (
        isolate_df.groupby("region")[
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

        # Return all mutations which pass the min_freq threshold, sorted
        # in ascending order by mutation ID
        # Also include the counts and the fraction of counts relative to
        # the number of genomes for this group
        return sorted([(int(k), v, (v / n_isolates) * 100) for k, v in counts.items()])

    # Do this column-by-column because for some reason pandas applymap()
    # misses the first column. I have no idea why...
    for col in ["dna_mutation", "gene_aa_mutation", "protein_aa_mutation"]:
        collapsed_region_mutations[col] = collapsed_region_mutations[col].apply(
            count_consensus
        )

    def mutations_to_df(field):
        _df = collapsed_region_mutations[field].explode()
        _df = pd.DataFrame.from_records(
            _df, columns=["mutation_id", "count", "freq"], index=_df.index
        )
        return _df

    collapsed_dna_mutations = mutations_to_df("dna_mutation")
    collapsed_gene_aa_mutations = mutations_to_df("gene_aa_mutation")
    collapsed_protein_aa_mutations = mutations_to_df("protein_aa_mutation")

    # DNA
    dna_mutation_region = pd.pivot_table(
        collapsed_dna_mutations.reset_index(),
        index="mutation_id",
        columns="region",
        values=["count", "freq"],
        aggfunc="first",  # should never trigger, all pairs should be unique
    ).fillna(0)
    dna_mutation_region.columns = [
        col[1] + " " + col[0] for col in dna_mutation_region.columns
    ]

    count_cols = [col for col in dna_mutation_region.columns if "count" in col]
    dna_mutation_region[count_cols] = dna_mutation_region[count_cols].astype(int)
    dna_mutation_region.insert(
        0, "sum_counts", dna_mutation_region[count_cols].apply(np.sum, axis=1)
    )

    dna_mutation_region = dna_mutation_def.join(dna_mutation_region, how="inner")
    dna_mutation_region = dna_mutation_region.sort_values("sum_counts", ascending=False)
    dna_mutation_region.to_csv(report_out / "dna_mutation_region.csv", index=False)

    # GENE AA
    gene_aa_mutation_region = pd.pivot_table(
        collapsed_gene_aa_mutations.reset_index(),
        index="mutation_id",
        columns="region",
        values=["count", "freq"],
        aggfunc="first",  # should never trigger, all pairs should be unique
    ).fillna(0)
    gene_aa_mutation_region.columns = [
        col[1] + " " + col[0] for col in gene_aa_mutation_region.columns
    ]

    count_cols = [col for col in gene_aa_mutation_region.columns if "count" in col]
    gene_aa_mutation_region[count_cols] = gene_aa_mutation_region[count_cols].astype(
        int
    )
    gene_aa_mutation_region.insert(
        0, "sum_counts", gene_aa_mutation_region[count_cols].apply(np.sum, axis=1)
    )

    gene_aa_mutation_region = gene_aa_mutation_def.join(
        gene_aa_mutation_region, how="inner"
    )
    gene_aa_mutation_region = gene_aa_mutation_region.sort_values(
        "sum_counts", ascending=False
    )
    gene_aa_mutation_region.to_csv(
        report_out / "gene_aa_mutation_region.csv", index=False
    )

    gene_aa_mutation_region.loc[
        gene_aa_mutation_region["gene"] == args.report_gene
    ].to_csv(report_out / f"{args.report_gene}_mutation_region.csv", index=False)

    # PROTEIN AA
    protein_aa_mutation_region = pd.pivot_table(
        collapsed_protein_aa_mutations.reset_index(),
        index="mutation_id",
        columns="region",
        values=["count", "freq"],
        aggfunc="first",  # should never trigger, all pairs should be unique
    ).fillna(0)
    protein_aa_mutation_region.columns = [
        col[1] + " " + col[0] for col in protein_aa_mutation_region.columns
    ]

    count_cols = [col for col in protein_aa_mutation_region.columns if "count" in col]
    protein_aa_mutation_region[count_cols] = protein_aa_mutation_region[
        count_cols
    ].astype(int)
    protein_aa_mutation_region.insert(
        0, "sum_counts", protein_aa_mutation_region[count_cols].apply(np.sum, axis=1)
    )

    protein_aa_mutation_region = protein_aa_mutation_def.join(
        protein_aa_mutation_region, how="inner"
    )
    protein_aa_mutation_region = protein_aa_mutation_region.sort_values(
        "sum_counts", ascending=False
    )
    protein_aa_mutation_region.to_csv(
        report_out / "protein_aa_mutation_region.csv", index=False
    )

    # -------------------------------------
    #  CO-OCCURRING REPORT GENE MUTATIONS
    # -------------------------------------

    isolate_df["dna_mutation"] = (
        isolate_df["dna_mutation"].apply(np.unique).apply(tuple)
    )
    isolate_df["gene_aa_mutation"] = (
        isolate_df["gene_aa_mutation"].apply(np.unique).apply(tuple)
    )
    isolate_df["protein_aa_mutation"] = (
        isolate_df["protein_aa_mutation"].apply(np.unique).apply(tuple)
    )

    isolate_df_gene_mutation = isolate_df[
        ["isolate_id", "gene_aa_mutation", args.report_group_col, "region"]
    ].explode("gene_aa_mutation")

    # Filter for only report gene mutations
    isolate_df_report_gene_mutation = (
        isolate_df_gene_mutation.reset_index(drop=True).join(
            (
                gene_aa_mutation_def.loc[
                    gene_aa_mutation_def["gene"] == args.report_gene
                ]
            ),
            on="gene_aa_mutation",
            how="inner",
        )
        # Sort this dataframe by position first so that
        # the final co-occuring mutation string is in the right order
        .sort_values("pos")
    )

    report_gene_cooc = isolate_df_report_gene_mutation.groupby("isolate_id").agg(
        report_gene_mutation=("mutation_name", lambda x: ":".join(x)),
        region=("region", "first"),
        **{args.report_group_col: (args.report_group_col, "first")},
    )

    # GLOBAL
    report_gene_cooc_global = (
        report_gene_cooc.reset_index()
        .groupby("report_gene_mutation", as_index=False)
        .agg(count=("isolate_id", "count"))
        .assign(freq=lambda x: (x["count"] / len(isolate_df)) * 100)
        .sort_values("count", ascending=False)
        .reset_index(drop=True)
    )
    report_gene_cooc_global_group = (
        report_gene_cooc.reset_index()
        .groupby(["report_gene_mutation", args.report_group_col], as_index=False)
        .agg(count=("isolate_id", "count"))
        .assign(freq=lambda x: (x["count"] / len(isolate_df)) * 100)
        .sort_values("count", ascending=False)
        .reset_index(drop=True)
    )
    report_gene_cooc_global.to_csv(
        report_out / f"{args.report_gene}_cooc_global.csv", index=False
    )
    report_gene_cooc_global_group.to_csv(
        report_out / f"{args.report_gene}_cooc_{args.report_group_col}_global.csv",
        index=False,
    )

    # REGIONAL

    # Get total region counts for normalization
    region_counts = (
        isolate_df.groupby("region")
        .agg(count=("isolate_id", "count"))
        .squeeze()
        .to_dict()
    )

    report_gene_cooc_region = (
        report_gene_cooc.reset_index()
        .groupby(["report_gene_mutation", "region"], as_index=False)
        .agg(count=("isolate_id", "count"))
        .assign(
            region_counts=lambda x: x["region"].map(region_counts),
            freq=lambda x: (x["count"] / x["region_counts"]) * 100,
        )
        .drop(columns=["region_counts"])
        .sort_values("count", ascending=False)
        .reset_index(drop=True)
    )
    report_gene_cooc_region = pd.pivot_table(
        report_gene_cooc_region,
        index="report_gene_mutation",
        columns="region",
        values=["count", "freq"],
    ).fillna(0)
    report_gene_cooc_region.columns = [
        col[1] + " " + col[0] for col in report_gene_cooc_region.columns
    ]
    count_cols = [col for col in report_gene_cooc_region.columns if "count" in col]
    report_gene_cooc_region[count_cols] = report_gene_cooc_region[count_cols].astype(
        int
    )
    report_gene_cooc_region.insert(
        0, "sum_counts", report_gene_cooc_region[count_cols].apply(np.sum, axis=1)
    )
    report_gene_cooc_region = report_gene_cooc_region.sort_values(
        "sum_counts", ascending=False
    )
    report_gene_cooc_region.reset_index(inplace=True)

    report_gene_cooc_region_group = (
        report_gene_cooc.reset_index()
        .groupby(
            ["report_gene_mutation", args.report_group_col, "region"], as_index=False
        )
        .agg(count=("isolate_id", "count"))
        .assign(
            region_counts=lambda x: x["region"].map(region_counts),
            freq=lambda x: (x["count"] / x["region_counts"]) * 100,
        )
        .drop(columns=["region_counts"])
        .sort_values("count", ascending=False)
        .reset_index(drop=True)
    )
    report_gene_cooc_region_group = pd.pivot_table(
        report_gene_cooc_region_group,
        index=["report_gene_mutation", args.report_group_col],
        columns="region",
        values=["count", "freq"],
    ).fillna(0)
    report_gene_cooc_region_group.columns = [
        col[1] + " " + col[0] for col in report_gene_cooc_region_group.columns
    ]
    count_cols = [
        col for col in report_gene_cooc_region_group.columns if "count" in col
    ]
    report_gene_cooc_region_group[count_cols] = report_gene_cooc_region_group[
        count_cols
    ].astype(int)
    report_gene_cooc_region_group.insert(
        0, "sum_counts", report_gene_cooc_region_group[count_cols].apply(np.sum, axis=1)
    )
    report_gene_cooc_region_group = report_gene_cooc_region_group.sort_values(
        "sum_counts", ascending=False
    )
    report_gene_cooc_region_group.reset_index(inplace=True)

    report_gene_cooc_region.to_csv(
        report_out / f"{args.report_gene}_cooc_region.csv", index=False
    )
    report_gene_cooc_region_group.to_csv(
        report_out / f"{args.report_gene}_cooc_{args.report_group_col}_region.csv",
        index=False,
    )

    # ------------------------
    #     REPORT GROUPING
    # ------------------------

    group_global = (
        isolate_df.groupby(args.report_group_col)
        .agg(count=("isolate_id", "count"))
        .assign(freq=lambda x: (x["count"] / len(isolate_df)) * 100)
        .sort_values("count", ascending=False)
        .reset_index()
    )
    group_global.to_csv(report_out / f"{args.report_group_col}_global.csv", index=False)

    group_region = (
        isolate_df.groupby([args.report_group_col, "region"], as_index=False)
        .agg(count=("isolate_id", "count"))
        .assign(
            region_counts=lambda x: x["region"].map(region_counts),
            freq=lambda x: (x["count"] / x["region_counts"]) * 100,
        )
    )

    group_region = pd.pivot_table(
        group_region,
        index=[args.report_group_col],
        columns="region",
        values=["count", "freq"],
    ).fillna(0)
    group_region.columns = [col[1] + " " + col[0] for col in group_region.columns]
    count_cols = [col for col in group_region.columns if "count" in col]
    group_region[count_cols] = group_region[count_cols].astype(int)
    group_region.insert(0, "sum_counts", group_region[count_cols].apply(np.sum, axis=1))
    group_region = group_region.sort_values("sum_counts", ascending=False)
    group_region.reset_index(inplace=True)
    group_region.to_csv(report_out / f"{args.report_group_col}_region.csv", index=False)


if __name__ == "__main__":
    main()
