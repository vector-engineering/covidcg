#!/usr/bin/env python3
# coding: utf-8

"""Generate main case data file for the UI

Author: Albert Chen (Deverman Lab, Broad Institute)
"""

import hashlib
import json
import numpy as np
import networkx as nx
import pandas as pd
import re

from functools import reduce
from pathlib import Path

from fasta import read_fasta_file
from clean_patient_metadata import clean_patient_metadata
from clean_seq_metadata import clean_seq_metadata
from process_ack import process_ack
from process_lineages import get_consensus_snps
from process_locations import process_location_data
from process_snps import load_dna_snps, load_aa_snps, process_snp_data
from reference import ref_seq, genes, gene_aa
from util import translate, data_dir


def load_patient_metadata():
    """Load patient metadata files
    """

    patient_meta_files = sorted((data_dir / "patient_meta").glob("*.tsv"))
    print(
        "Collecting {} patient metadata files...".format(len(patient_meta_files)),
        end="",
        flush=True,
    )
    patient_meta_df = pd.DataFrame()
    for f in patient_meta_files:
        _df = pd.read_csv(f, sep="\t", skiprows=2)
        patient_meta_df = pd.concat([patient_meta_df, _df], ignore_index=True)

    # Drop columns we don't need
    patient_meta_df = patient_meta_df.drop(
        columns=[
            "Virus name",
            "Host",
            "Additional location information",
            "Additional host information",
        ]
    )

    # Set Accession ID as the index
    patient_meta_df = patient_meta_df.set_index("Accession ID")

    # Save dataframe
    # patient_meta_df.to_csv(data_dir / 'patient_meta.csv', index=False)
    print("done", flush=True)

    return patient_meta_df


def load_seq_metadata():
    """Load sequencing technology metadata files
    """

    seq_meta_files = sorted((data_dir / "seq_meta").glob("*.tsv"))
    print(
        "Collecting {} sequencing technology metadata files...".format(
            len(seq_meta_files)
        ),
        end="",
        flush=True,
    )

    seq_meta_df = pd.DataFrame()
    for f in seq_meta_files:
        _df = pd.read_csv(f, sep="\t", skiprows=2)
        seq_meta_df = pd.concat([seq_meta_df, _df], ignore_index=True)

    # Drop columns we don't need
    # A lot of these we'll get from the patient metadata files instead
    seq_meta_df = seq_meta_df.drop(
        columns=[
            "Virus name",
            "Collection date",
            "Location",
            "Host",
            "Passage",
            "Specimen",
            "Additional host information",
            "Lineage",
            "Clade",
            "Comment",
        ]
    )

    # Set Accession ID as index
    seq_meta_df = seq_meta_df.set_index("Accession ID")

    print("done", flush=True)

    return seq_meta_df


def hash_accession_id(accession_id):
    m = hashlib.sha256()
    m.update(str(accession_id).encode("utf-8"))
    return m.hexdigest()


def write_reference_files():
    print("Writing reference sequence files...", end="", flush=True)
    # Write the reference fasta file to json
    ref_json_path = data_dir / "reference.json"

    ref_obj = {"ref_seq": ref_seq, "gene_aa": gene_aa}

    with ref_json_path.open("w") as fp:
        fp.write(json.dumps(ref_obj))

    print("done")


def main():

    # Load patient metadata, clean up
    patient_meta_df = load_patient_metadata()
    patient_meta_df = clean_patient_metadata(patient_meta_df)

    # Load sequencing metadata, clean up
    seq_meta_df = load_seq_metadata()
    seq_meta_df = clean_seq_metadata(seq_meta_df)
    # print(seq_meta_df)

    # Join patient and sequencing metadata on Accession ID
    case_df = patient_meta_df.join(
        seq_meta_df, on="Accession ID", how="left", sort=True
    )

    # Filter out "None" lineages
    case_df = case_df.loc[case_df["lineage"] != "None", :]

    # Process location data
    # Unset the Accession ID as the index, to speed up filtering
    location_df, unique_location_df = process_location_data(case_df.reset_index())

    # Join location IDs onto main metadata dataframe
    case_df = case_df.join(
        location_df[["location_id"]], on="Accession ID", how="inner", sort=False
    )
    # Drop original Location column
    case_df = case_df.drop(columns=["Location"])

    # Join acknowledgement IDs onto main metadata dataframe
    ack_df = process_ack()
    print("Joining acknowledgements to main dataframe...", end="", flush=True)
    case_df = case_df.join(ack_df, on="Accession ID", how="left", sort=False)
    # Replace missing acknowledgement IDs with -1
    case_df["ack_id"].fillna(-1, inplace=True)
    # Cast ack_id to integer
    case_df["ack_id"] = case_df["ack_id"].astype(int)
    print("done")

    dna_snp_df = load_dna_snps()
    aa_snp_df = load_aa_snps()
    dna_snp_group_df, aa_snp_group_df = process_snp_data(dna_snp_df, aa_snp_df)

    print("Joining SNP data to main dataframe...", end="", flush=True)
    dna_snp_group_df = dna_snp_group_df.rename(
        columns={"gisaid_id": "Accession ID"}
    ).set_index("Accession ID")
    aa_snp_group_df = aa_snp_group_df.rename(
        columns={"gisaid_id": "Accession ID"}
    ).set_index("Accession ID")

    # Inner join to main dataframe, to exclude filtered out sequences
    case_df = case_df.join(
        dna_snp_group_df[["snp_str", "snp_sig"]],
        on="Accession ID",
        how="inner",
        sort=False,
    ).rename(columns={"snp_str": "dna_snp_str", "snp_sig": "dna_snp_sig"})
    case_df = case_df.join(
        aa_snp_group_df[["snp_str", "snp_sig"]],
        on="Accession ID",
        how="inner",
        sort=False,
    ).rename(columns={"snp_str": "aa_snp_str", "snp_sig": "aa_snp_sig"})
    print("done")

    # Get consensus SNPs for each lineage
    print("Getting consensus SNPs for each lineage...", end="", flush=True)
    get_consensus_snps(case_df)
    print("done")

    # Hash Accession IDs
    print("Anonymizing/hashing accession IDs...", end="", flush=True)
    case_df["hashed_id"] = case_df.index.to_series().apply(hash_accession_id)

    # Create map of hash -> Accession ID
    accession_hash_df = case_df[["hashed_id"]]
    accession_hash_df.to_csv(
        data_dir / "accession_hashmap.csv", index_label="Accession ID"
    )

    # Delete old accession ID column, reassign to hashed ID
    case_df = (
        case_df.reset_index()
        .drop(columns=["Accession ID"])
        .rename(columns={"hashed_id": "Accession ID"})
        .set_index("Accession ID")
    )
    print("done")

    # Factorize some more metadata columns
    print("Mapping more metadata columns...", end="", flush=True)
    map_cols = [
        "gender",
        "patient_status",
        "passage",
        "specimen",
        "sequencing_tech",
        "assembly_method",
        "comment_type",
    ]
    metadata_maps = {}

    for i, col in enumerate(map_cols):
        factor = pd.factorize(case_df[col])

        id_col = col + "_id"
        case_df[id_col] = factor[0]

        metadata_maps[col] = pd.Series(factor[1]).to_dict()

    # Drop the original metadata columns
    case_df = case_df.drop(columns=map_cols)

    # Write the metadata map to a JSON file
    metadata_map_path = data_dir / "metadata_map.json"
    with metadata_map_path.open("w") as fp:
        fp.write(json.dumps(metadata_maps))
    print("done")

    print("Writing final case data...", end="", flush=True)
    case_df.to_csv(data_dir / "case_data2.csv", index_label="Accession ID")
    case_df.reset_index().to_json(data_dir / "case_data2.json", orient="records")
    print("done")

    write_reference_files()


if __name__ == "__main__":
    main()

