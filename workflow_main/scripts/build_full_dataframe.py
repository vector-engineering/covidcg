# coding: utf-8

import pandas as pd
import numpy as np
import json


def build_full_dataframe(case_data, metadata_map, df_out):

    # Load data
    df = pd.read_json(case_data).set_index("Accession ID")
    with open(metadata_map, "r") as fp:
        metadata_map = json.loads(fp.read())

    # Join SNV information
    dna_snp_map = {v: k for k, v in metadata_map["dna_snp"].items()}
    gene_aa_snp_map = {v: k for k, v in metadata_map["gene_aa_snp"].items()}
    protein_aa_snp_map = {v: k for k, v in metadata_map["protein_aa_snp"].items()}

    df.loc[:, "dna_snp_str"] = df["dna_snp_str"].apply(
        lambda x: ";".join([dna_snp_map[i] for i in x])
    )
    df.loc[:, "gene_aa_snp_str"] = df["gene_aa_snp_str"].apply(
        lambda x: ";".join([gene_aa_snp_map[i] for i in x])
    )
    df.loc[:, "protein_aa_snp_str"] = df["protein_aa_snp_str"].apply(
        lambda x: ";".join([protein_aa_snp_map[i] for i in x])
    )

    df = df.rename(
        columns={
            "dna_snp_str": "dna_snp",
            "gene_aa_snp_str": "gene_aa_snp",
            "protein_aa_snp_str": "protein_aa_snp",
        }
    )

    # Metadata
    for col in metadata_map.keys():
        # We already did these
        if col in ["dna_snp", "gene_aa_snp", "protein_aa_snp"]:
            continue

        mmap = {int(k): v for k, v in metadata_map[col].items()}
        df[col] = df[col].map(mmap)

    df.to_csv(df_out)

