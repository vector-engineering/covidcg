#!/usr/bin/env python3
# coding: utf-8

"""Combine coverage data -- PIVOT VERSION -- UNUSED

python3 scripts/combine_coverage.py \
    --manifest ../example_data_genbank/sars2/sequence_manifest.csv \
    --reference ../static_data/sars2/reference.json \
    --mode dna \
    --coverage-dir ../example_data_genbank/sars2/coverage_dna \
    --out test.csv

python3 scripts/combine_coverage.py \
    --manifest ../example_data_genbank/sars2/sequence_manifest.csv \
    --reference ../static_data/sars2/reference.json \
    --gene-protein-def ../static_data/sars2/genes_processed.json \
    --mode gene_aa \
    --coverage-dir ../example_data_genbank/sars2/coverage_gene_aa \
    --out test.csv

python3 scripts/combine_coverage.py \
    --manifest ../data_gisaid_flu/sequence_manifest.csv \
    --reference ../static_data/flu/reference.json \
    --mode dna \
    --coverage-dir ../data_gisaid_flu/coverage_dna \
    --out test.csv

python3 scripts/combine_coverage.py \
    --manifest ../data_gisaid_flu/sequence_manifest.csv \
    --reference ../static_data/flu/reference.json \
    --gene-protein-def ../static_data/flu/genes_processed.json \
    --mode gene_aa \
    --coverage-dir ../data_gisaid_flu/coverage_gene_aa \
    --out test.csv


Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""

import argparse
import io
import json
from pathlib import Path

import pandas as pd


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--manifest", type=str, required=True, help="Path to manifest CSV file"
    )
    parser.add_argument(
        "--reference", type=str, required=True, help="Path to reference file"
    )
    parser.add_argument(
        "--gene-protein-def", type=str, help="Path to gene/protein definition JSON",
    )
    parser.add_argument(
        "--coverage-dir", type=str, required=True, help="Coverage files",
    )
    parser.add_argument(
        "--mode", type=str, required=True, help="Mode (dna, gene_aa, protein_aa)"
    )

    # Output
    parser.add_argument(
        "--out", type=str, required=True, help="Combined coverage file CSV output"
    )

    args = parser.parse_args()

    # Load the reference sequences
    with open(args.reference, "r") as fp:
        references = json.loads(fp.read())

    # Build list of all feature columns
    # For DNA mode - columns are just the possible segments
    # For AA mode - all possible ORFs across all references
    feature_cols = set()
    if args.mode == "dna":
        for reference in references.keys():
            for segment in references[reference]["segments"].keys():
                feature_cols.add(f"start_{segment}")
                feature_cols.add(f"end_{segment}")
    else:
        # Load gene/protein defs
        # JSON to dataframe
        with open(args.gene_protein_def, "r") as fp:
            feature_dicts = json.loads(fp.read())

        for reference in feature_dicts.keys():
            for feature in feature_dicts[reference]:
                if feature["protein_coding"]:
                    feature_cols.add(f"start_{feature['name']}")
                    feature_cols.add(f"end_{feature['name']}")
    feature_cols = list(sorted(feature_cols))

    # Read in chunks by date, and combine segment files
    # This is necessary since we want all segments (or features, for gene/protein mode)
    # on the same row. Therefore, we need all segment files per sequence
    # to be in memory before we can do the pivot operation
    chunk_paths = sorted(Path(args.coverage_dir).glob("*.csv"))
    chunks_by_segment = pd.DataFrame(chunk_paths, columns=["path"])
    chunks_by_segment["name"] = chunks_by_segment["path"].apply(
        lambda x: x.name.replace("_coverage_" + args.mode + ".csv", "")
    )
    chunks_by_segment = pd.concat(
        [
            chunks_by_segment,
            (
                chunks_by_segment["name"]
                .str.extract(r"^([0-9]+)_([A-Za-z0-9-]+)_([0-9-]+)$")
                .rename(columns={0: "segment", 1: "subtype", 2: "date"})
            ),
        ],
        axis=1,
    )
    chunks_by_segment = chunks_by_segment.groupby(
        ["subtype", "date"], as_index=False
    ).agg(paths=("path", list), names=("name", list), segments=("segment", list))

    # Dump all mutation chunks into a text buffer
    # This should be more memory efficient than storing a list of dataframes
    # then doing a huge pd.concat() operation
    coverage_df_io = io.StringIO()
    for _, row in chunks_by_segment.iterrows():
        paths = row["paths"]
        names = row["names"]
        segments = row["segments"]
        # Combine segment csvs into one dataframe
        chunk_df = pd.concat(
            [
                pd.read_csv(path).assign(
                    file_name=lambda _: name, segment=lambda _: segment
                )
                for path, name, segment in zip(paths, names, segments)
            ],
            axis=0,
            ignore_index=True,
            sort=False,
            copy=False,
        )

        # For DNA mode, the segment column is the feature column
        if args.mode == "dna":
            chunk_df["subtype"] = row["subtype"]
            chunk_df.rename(columns={"segment": "feature"}, inplace=True)
        else:
            chunk_df.drop(columns=["segment"], inplace=True)

        # Pivot for one row per sequence
        chunk_df = pd.pivot_table(
            chunk_df,
            index=["Accession ID", "subtype", "reference", "file_name"],
            columns=["feature"],
            values=["start", "end"],
            sort=False,
        )

        # Unset multi-index columns
        chunk_df.columns = chunk_df.columns.map("_".join)
        # Add missing feature columns, then sort columns
        chunk_df = chunk_df.reindex(columns=feature_cols, fill_value=None)
        chunk_df.reset_index(inplace=True)
        chunk_df.to_csv(coverage_df_io, index=False, header=False)

    coverage_df_io.seek(0)

    # For gene/protein mode, read segment as string
    dtype_opts = {}
    if args.mode == "gene_aa" or args.mode == "protein_aa":
        # dtype_opts["segment"] = str
        dtype_opts["feature"] = str

    # Specify only '' for null values
    # Otherwise gene name of 'NA' will be interpreted as missing
    coverage_df = pd.read_csv(
        coverage_df_io,
        header=None,
        dtype=dtype_opts,
        keep_default_na=False,
        na_values=[""],
        names=["Accession ID", "subtype", "reference", "file_name"] + feature_cols,
    )
    coverage_df_io.close()

    # --------------------------
    # Remove duplicate sequences
    # --------------------------

    # Get sequence manifest
    manifest = pd.read_csv(args.manifest)

    # The segment/subtype cols are provided by the coverage files themselves
    # Remove now to prevent duplicate columns during the join
    manifest.drop(columns=["segment", "subtype"], inplace=True)

    # Also has the effect of adding rows for sequences without mutations
    # (pos, ref, alt, etc filled with NaNs)

    coverage_df = manifest.merge(
        coverage_df,
        how="left",
        left_on=["Accession ID", "reference", "file_name"],
        right_on=["Accession ID", "reference", "file_name"],
    )

    # Remove rows without coverage data
    coverage_df.drop(
        coverage_df.index[
            coverage_df[feature_cols].isna().sum(axis=1) == len(feature_cols)
        ],
        inplace=True,
    )
    coverage_df.reset_index(drop=True, inplace=True)

    # Drop unnecessary metadata columns
    coverage_df.drop(columns=["file_name", "date"], inplace=True)

    coverage_df.to_csv(args.out, index=False)


if __name__ == "__main__":
    main()

