#!/usr/bin/env python3
# coding: utf-8

"""Combine metadata and mutation data, create metadata maps

Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""

import argparse
import json
from pathlib import Path

import pandas as pd

from process_mutations import process_mutations


def main():
    parser = argparse.ArgumentParser()

    # Input
    parser.add_argument("--metadata", type=str, required=True, help="Metadata file")

    parser.add_argument(
        "--manifest", type=str, required=True, help="Path to manifest CSV file",
    )
    parser.add_argument(
        "--dna-mutation-dir", type=str, required=True, help="DNA mutation directory",
    )
    parser.add_argument(
        "--gene-aa-mutation-dir",
        type=str,
        required=True,
        help="Gene AA mutation directory",
    )
    parser.add_argument(
        "--protein-aa-mutation-dir",
        type=str,
        required=True,
        help="Protein AA mutation directory",
    )
    parser.add_argument(
        "--dna-coverage", type=str, required=True, help="DNA coverage CSV file"
    )
    parser.add_argument(
        "--gene-aa-coverage", type=str, required=True, help="Gene AA coverage CSV file"
    )
    parser.add_argument(
        "--protein-aa-coverage",
        type=str,
        required=True,
        help="Protein AA coverage CSV file",
    )

    # Output
    parser.add_argument(
        "--metadata-map", type=str, required=True, help="Metadata map file output"
    )
    parser.add_argument(
        "--sequence-data",
        type=str,
        required=True,
        help="Sequence data JSON file output",
    )
    parser.add_argument(
        "--sequence-data-csv",
        type=str,
        required=True,
        help="Sequence data CSV file output",
    )

    # Parameters
    parser.add_argument(
        "--count-threshold", type=int, default=3, help="Global mutation count threshold"
    )
    parser.add_argument(
        "--metadata-cols", type=str, nargs="+", default=[], help="Metadata columns"
    )

    args = parser.parse_args()

    manifest = pd.read_csv(args.manifest)

    # Count mutations
    dna_mutation_group_df, dna_mutation_map = process_mutations(
        manifest,
        sorted(Path(args.dna_mutation_dir).glob("*.csv")),
        mode="dna",
        count_threshold=args.count_threshold,
    )
    gene_aa_mutation_group_df, gene_aa_mutation_map = process_mutations(
        manifest,
        sorted(Path(args.gene_aa_mutation_dir).glob("*.csv")),
        mode="gene_aa",
        count_threshold=args.count_threshold,
    )
    protein_aa_mutation_group_df, protein_aa_mutation_map = process_mutations(
        manifest,
        sorted(Path(args.protein_aa_mutation_dir).glob("*.csv")),
        mode="protein_aa",
        count_threshold=args.count_threshold,
    )

    # Load metadata
    df = pd.read_csv(args.metadata)
    # If metadata contains subtype column, remove it to avoid duplicate column
    # during merge with manifest
    if "subtype" in df.columns:
        df.drop(columns=["subtype"], inplace=True)

    # Join sequence ID and reference from manifest
    df = manifest.merge(
        df.drop(columns=["segment"]), on="Accession ID", how="left"
    ).drop(columns=["file_name", "date"])

    # Exclude sequences without a group assignment
    # (i.e., lineage or clade assignment)
    # "group_cols" is defined in the "group_cols" field in the
    # config.yaml file
    # for col in group_cols:
    #     df = df.loc[~pd.isnull(df[col]), :]

    # Join mutations to main dataframe
    # inner join to exclude filtered out sequences
    df = (
        df.merge(
            dna_mutation_group_df.drop(columns=["sequence_id"]),
            left_on=["Accession ID", "reference"],
            right_on=["Accession ID", "reference"],
            how="inner",
        )
        .rename(columns={"mutation_id": "dna_mutation"})
        .merge(
            gene_aa_mutation_group_df.drop(columns=["sequence_id"]),
            left_on=["Accession ID", "reference"],
            right_on=["Accession ID", "reference"],
            how="inner",
        )
        .rename(columns={"mutation_id": "gene_aa_mutation"})
        .merge(
            protein_aa_mutation_group_df.drop(columns=["sequence_id"]),
            left_on=["Accession ID", "reference"],
            right_on=["Accession ID", "reference"],
            how="inner",
        )
        .rename(columns={"mutation_id": "protein_aa_mutation"})
    )

    # Join coverage data to main dataframe
    coverage_dna = pd.read_csv(args.dna_coverage, dtype={"segment": str},)
    # coverage_dna["range"] = coverage_dna[["start", "end"]].apply(
    #     lambda x: tuple(x), axis=1
    # )
    # Sometimes with multiple alignments to the same reference
    # we get multiple DNA coverage entries for the same
    # Accession ID / reference pair.
    # Collapse these duplicate entries, by assuming that the
    # sequence is contiguous -- i.e., take the max of the ends
    # and the min of the starts
    coverage_dna_group = (
        coverage_dna.groupby(["Accession ID", "reference"], as_index=False)
        .agg(start=("start", min), end=("end", max))
        .assign(range=lambda x: x[["start", "end"]].apply(lambda _x: tuple(_x), axis=1))
    )

    coverage_gene_aa = pd.read_csv(
        args.gene_aa_coverage,
        dtype={"segment": str, "feature": str},
        keep_default_na=False,
        na_values=[""],
    )
    coverage_gene_aa["range"] = coverage_gene_aa[["feature", "start", "end"]].apply(
        lambda x: tuple(x), axis=1
    )
    coverage_gene_aa_group = (
        coverage_gene_aa.drop(columns=["feature", "start", "end"])
        .groupby(["Accession ID", "reference"], as_index=False)
        .agg(range=("range", list))
    )

    coverage_protein_aa = pd.read_csv(
        args.protein_aa_coverage,
        dtype={"segment": str, "feature": str},
        keep_default_na=False,
        na_values=[""],
    )
    coverage_protein_aa["range"] = coverage_protein_aa[
        ["feature", "start", "end"]
    ].apply(lambda x: tuple(x), axis=1)
    coverage_protein_aa_group = (
        coverage_protein_aa.drop(columns=["feature", "start", "end"])
        .groupby(["Accession ID", "reference"], as_index=False)
        .agg(range=("range", list))
    )

    df = (
        df.merge(
            coverage_dna_group[["Accession ID", "reference", "range"]],
            left_on=["Accession ID", "reference"],
            right_on=["Accession ID", "reference"],
            how="left",
        )
        .rename(columns={"range": "dna_range"})
        .merge(
            coverage_gene_aa_group[["Accession ID", "reference", "range"]],
            left_on=["Accession ID", "reference"],
            right_on=["Accession ID", "reference"],
            how="left",
        )
        .rename(columns={"range": "gene_aa_range"})
        .merge(
            coverage_protein_aa_group[["Accession ID", "reference", "range"]],
            left_on=["Accession ID", "reference"],
            right_on=["Accession ID", "reference"],
            how="left",
        )
        .rename(columns={"range": "protein_aa_range"})
    )

    # Remove rows with null dna_range
    # These are failed alignments with no mutation information
    df.drop(df.index[df["dna_range"].isna()], inplace=True)

    # Factorize some more metadata columns
    metadata_maps = {}

    # Metadata cols passed in as kwarg, and defined
    # in config.yaml as "metadata_cols"
    for col in args.metadata_cols:
        factor, labels = pd.factorize(df[col])
        df.loc[:, col] = factor
        metadata_maps[col] = pd.Series(labels).to_dict()

    # Special processing for locations - leave missing data as -1
    for col in ["region", "country", "division", "location"]:
        missing_inds = df[col] == "-1"
        factor, labels = pd.factorize(df.loc[~missing_inds, col])
        df.loc[~missing_inds, col] = factor
        metadata_maps[col] = pd.Series(labels).to_dict()
        df.loc[:, col] = df[col].astype(int)

    # Add mutation maps into the metadata map
    metadata_maps["dna_mutation"] = dna_mutation_map.to_dict()
    metadata_maps["gene_aa_mutation"] = gene_aa_mutation_map.to_dict()
    metadata_maps["protein_aa_mutation"] = protein_aa_mutation_map.to_dict()

    # Write the metadata map to a JSON file
    with open(args.metadata_map, "w") as fp:
        fp.write(json.dumps(metadata_maps))

    # Write final dataframe
    df.to_csv(args.sequence_data_csv, index=False)
    df.to_json(args.sequence_data, orient="records")


if __name__ == "__main__":
    main()
