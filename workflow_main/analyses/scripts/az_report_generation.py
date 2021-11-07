# coding: utf-8

"""Table building for AZ

Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""

import pandas as pd
import numpy as np
import json

from collections import Counter
from itertools import chain
from pathlib import Path


def dna_mutation_str_to_name(mut_str):
    split = mut_str.split("|")
    # REF POS ALT
    return split[1] + split[0] + split[2]


def aa_mutation_str_to_name(mut_str):
    split = mut_str.split("|")
    # REF POS ALT
    return split[2] + split[1] + split[3]


def generate_az_reports(case_data_path, metadata_map_path, report_out):
    report_out = Path(report_out)

    # ------------------------
    #        LOAD DATA
    # ------------------------

    case_data = pd.read_json(case_data_path)

    # Load DNA mutation ID map
    with open(metadata_map_path, "r") as fp:
        metadata_map = json.loads(fp.read())

    id_to_dna_mutation = {v: k for k, v in metadata_map["dna_mutation"].items()}
    id_to_gene_aa_mutation = {v: k for k, v in metadata_map["gene_aa_mutation"].items()}
    id_to_protein_aa_mutation = {
        v: k for k, v in metadata_map["protein_aa_mutation"].items()
    }

    # Join locations onto case_data
    loc_levels = ["region", "country", "division", "location"]
    for loc_level in loc_levels:
        case_data.loc[:, loc_level] = case_data[loc_level].map(
            {int(k): v for k, v in metadata_map[loc_level].items()}
        )
        case_data.loc[case_data[loc_level].isna(), loc_level] = None

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
        1, "pos", dna_mutation_def["mutation_str"].apply(lambda x: int(x.split("|")[0]))
    )
    dna_mutation_def.insert(
        2, "ref", dna_mutation_def["mutation_str"].apply(lambda x: x.split("|")[1])
    )
    dna_mutation_def.insert(
        3, "alt", dna_mutation_def["mutation_str"].apply(lambda x: x.split("|")[2])
    )
    dna_mutation_def.insert(
        0,
        "mutation_name",
        dna_mutation_def["mutation_str"].apply(dna_mutation_str_to_name),
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
        gene_aa_mutation_def["mutation_str"].apply(aa_mutation_str_to_name),
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
        protein_aa_mutation_def["mutation_str"].apply(aa_mutation_str_to_name),
    )
    protein_aa_mutation_def.drop(columns=["mutation_str"], inplace=True)

    # -----------------------------
    #    GLOBAL SPIKE MUTATIONS
    # -----------------------------

    collapsed_mutations = (
        case_data.groupby(lambda x: True)[
            ["dna_mutation_str", "gene_aa_mutation_str", "protein_aa_mutation_str"]
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
        freq=(dna_mutation_df["count"] / len(case_data)) * 100
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
        freq=(gene_aa_mutation_df["count"] / len(case_data)) * 100
    )
    gene_aa_mutation_df = gene_aa_mutation_df.sort_values("count", ascending=False)
    gene_aa_mutation_df.to_csv(report_out / "gene_aa_mutation_global.csv", index=False)

    (
        gene_aa_mutation_df.loc[gene_aa_mutation_df["gene"] == "S"]
        .drop(columns=["gene"])
        .to_csv(report_out / "spike_mutation_global.csv", index=False)
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
        freq=(protein_aa_mutation_df["count"] / len(case_data)) * 100
    )
    protein_aa_mutation_df = protein_aa_mutation_df.sort_values(
        "count", ascending=False
    )
    protein_aa_mutation_df.to_csv(
        report_out / "protein_aa_mutation_global.csv", index=False
    )

    # ------------------------------
    #   REGIONAL SPIKE MUTATIONS
    # ------------------------------

    collapsed_region_mutations = (
        case_data.groupby("region")[
            ["dna_mutation_str", "gene_aa_mutation_str", "protein_aa_mutation_str"]
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
    for col in ["dna_mutation_str", "gene_aa_mutation_str", "protein_aa_mutation_str"]:
        collapsed_region_mutations[col] = collapsed_region_mutations[col].apply(
            count_consensus
        )

    def mutations_to_df(field):
        _df = collapsed_region_mutations[field].explode()
        _df = pd.DataFrame.from_records(
            _df, columns=["mutation_id", "count", "freq"], index=_df.index
        )
        return _df

    collapsed_dna_mutations = mutations_to_df("dna_mutation_str")
    collapsed_gene_aa_mutations = mutations_to_df("gene_aa_mutation_str")
    collapsed_protein_aa_mutations = mutations_to_df("protein_aa_mutation_str")

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

    gene_aa_mutation_region.loc[gene_aa_mutation_region["gene"] == "S"].to_csv(
        report_out / "spike_mutation_region.csv", index=False
    )

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

    # -------------------------------
    #  CO-OCCURRING SPIKE MUTATIONS
    # -------------------------------

    case_data["dna_mutation_str"] = (
        case_data["dna_mutation_str"].apply(np.unique).apply(tuple)
    )
    case_data["gene_aa_mutation_str"] = (
        case_data["gene_aa_mutation_str"].apply(np.unique).apply(tuple)
    )
    case_data["protein_aa_mutation_str"] = (
        case_data["protein_aa_mutation_str"].apply(np.unique).apply(tuple)
    )

    case_data_gene_mutation = case_data[
        ["Accession ID", "gene_aa_mutation_str", "lineage", "region"]
    ].explode("gene_aa_mutation_str")

    # Filter for only spike mutations
    case_data_spike_mutation = (
        case_data_gene_mutation.reset_index(drop=True).join(
            (gene_aa_mutation_def.loc[gene_aa_mutation_def["gene"] == "S"]),
            on="gene_aa_mutation_str",
            how="inner",
        )
        # Sort this dataframe by position first so that
        # the final co-occuring mutation string is in the right order
        .sort_values("pos")
    )

    spike_cooc = case_data_spike_mutation.groupby("Accession ID").agg(
        spike_mutation=("mutation_name", lambda x: ":".join(x)),
        region=("region", "first"),
        lineage=("lineage", "first"),
    )

    # GLOBAL
    spike_cooc_global = (
        spike_cooc.reset_index()
        .groupby("spike_mutation", as_index=False)
        .agg(count=("Accession ID", "count"))
        .assign(freq=lambda x: (x["count"] / len(case_data)) * 100)
        .sort_values("count", ascending=False)
        .reset_index(drop=True)
    )
    spike_cooc_global_lineage = (
        spike_cooc.reset_index()
        .groupby(["spike_mutation", "lineage"], as_index=False)
        .agg(count=("Accession ID", "count"))
        .assign(freq=lambda x: (x["count"] / len(case_data)) * 100)
        .sort_values("count", ascending=False)
        .reset_index(drop=True)
    )
    spike_cooc_global.to_csv(report_out / "spike_cooc_global.csv", index=False)
    spike_cooc_global_lineage.to_csv(
        report_out / "spike_cooc_lineage_global.csv", index=False
    )

    # REGIONAL

    # Get total region counts for normalization
    region_counts = (
        case_data.groupby("region")
        .agg(count=("Accession ID", "count"))
        .squeeze()
        .to_dict()
    )

    spike_cooc_region = (
        spike_cooc.reset_index()
        .groupby(["spike_mutation", "region"], as_index=False)
        .agg(count=("Accession ID", "count"))
        .assign(
            region_counts=lambda x: x["region"].map(region_counts),
            freq=lambda x: (x["count"] / x["region_counts"]) * 100,
        )
        .drop(columns=["region_counts"])
        .sort_values("count", ascending=False)
        .reset_index(drop=True)
    )
    spike_cooc_region = pd.pivot_table(
        spike_cooc_region,
        index="spike_mutation",
        columns="region",
        values=["count", "freq"],
    ).fillna(0)
    spike_cooc_region.columns = [
        col[1] + " " + col[0] for col in spike_cooc_region.columns
    ]
    count_cols = [col for col in spike_cooc_region.columns if "count" in col]
    spike_cooc_region[count_cols] = spike_cooc_region[count_cols].astype(int)
    spike_cooc_region.insert(
        0, "sum_counts", spike_cooc_region[count_cols].apply(np.sum, axis=1)
    )
    spike_cooc_region = spike_cooc_region.sort_values("sum_counts", ascending=False)
    spike_cooc_region.reset_index(inplace=True)

    spike_cooc_region_lineage = (
        spike_cooc.reset_index()
        .groupby(["spike_mutation", "lineage", "region"], as_index=False)
        .agg(count=("Accession ID", "count"))
        .assign(
            region_counts=lambda x: x["region"].map(region_counts),
            freq=lambda x: (x["count"] / x["region_counts"]) * 100,
        )
        .drop(columns=["region_counts"])
        .sort_values("count", ascending=False)
        .reset_index(drop=True)
    )
    spike_cooc_region_lineage = pd.pivot_table(
        spike_cooc_region_lineage,
        index=["spike_mutation", "lineage"],
        columns="region",
        values=["count", "freq"],
    ).fillna(0)
    spike_cooc_region_lineage.columns = [
        col[1] + " " + col[0] for col in spike_cooc_region_lineage.columns
    ]
    count_cols = [col for col in spike_cooc_region_lineage.columns if "count" in col]
    spike_cooc_region_lineage[count_cols] = spike_cooc_region_lineage[
        count_cols
    ].astype(int)
    spike_cooc_region_lineage.insert(
        0, "sum_counts", spike_cooc_region_lineage[count_cols].apply(np.sum, axis=1)
    )
    spike_cooc_region_lineage = spike_cooc_region_lineage.sort_values(
        "sum_counts", ascending=False
    )
    spike_cooc_region_lineage.reset_index(inplace=True)

    spike_cooc_region.to_csv(report_out / "spike_cooc_region.csv", index=False)
    spike_cooc_region_lineage.to_csv(
        report_out / "spike_cooc_lineage_region.csv", index=False
    )

    # ------------------------
    #        LINEAGES
    # ------------------------

    lineage_global = (
        case_data.groupby("lineage")
        .agg(count=("Accession ID", "count"))
        .assign(freq=lambda x: (x["count"] / len(case_data)) * 100)
        .sort_values("count", ascending=False)
        .reset_index()
    )
    lineage_global.to_csv(report_out / "lineage_global.csv", index=False)

    lineage_region = (
        case_data.groupby(["lineage", "region"], as_index=False)
        .agg(count=("Accession ID", "count"))
        .assign(
            region_counts=lambda x: x["region"].map(region_counts),
            freq=lambda x: (x["count"] / x["region_counts"]) * 100,
        )
    )

    lineage_region = pd.pivot_table(
        lineage_region, index=["lineage"], columns="region", values=["count", "freq"],
    ).fillna(0)
    lineage_region.columns = [col[1] + " " + col[0] for col in lineage_region.columns]
    count_cols = [col for col in lineage_region.columns if "count" in col]
    lineage_region[count_cols] = lineage_region[count_cols].astype(int)
    lineage_region.insert(
        0, "sum_counts", lineage_region[count_cols].apply(np.sum, axis=1)
    )
    lineage_region = lineage_region.sort_values("sum_counts", ascending=False)
    lineage_region.reset_index(inplace=True)
    lineage_region.to_csv(report_out / "lineage_region.csv", index=False)

