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
<<<<<<< HEAD
    dna_mutation_map = {v: k for k, v in metadata_map["dna_mutation"].items()}
    gene_aa_mutation_map = {v: k for k, v in metadata_map["gene_aa_mutation"].items()}
    protein_aa_mutation_map = {v: k for k, v in metadata_map["protein_aa_mutation"].items()}

    df.loc[:, "dna_mutation_str"] = df["dna_mutation_str"].apply(
        lambda x: ";".join([dna_mutation_map[i] for i in x])
    )
    df.loc[:, "gene_aa_mutation_str"] = df["gene_aa_mutation_str"].apply(
        lambda x: ";".join([gene_aa_mutation_map[i] for i in x])
    )
    df.loc[:, "protein_aa_mutation_str"] = df["protein_aa_mutation_str"].apply(
        lambda x: ";".join([protein_aa_mutation_map[i] for i in x])
=======
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
>>>>>>> e6dd8312 (Rsvg workflow main (#420))
    )

    df = df.rename(
        columns={
<<<<<<< HEAD
            "dna_mutation_str": "dna_mutation",
            "gene_aa_mutation_str": "gene_aa_mutation",
            "protein_aa_mutation_str": "protein_aa_mutation",
=======
            "dna_snp_str": "dna_snp",
            "gene_aa_snp_str": "gene_aa_snp",
            "protein_aa_snp_str": "protein_aa_snp",
>>>>>>> e6dd8312 (Rsvg workflow main (#420))
        }
    )

    # Metadata
    for col in metadata_map.keys():
        # We already did these
<<<<<<< HEAD
        if col in ["dna_mutation", "gene_aa_mutation", "protein_aa_mutation"]:
=======
        if col in ["dna_snp", "gene_aa_snp", "protein_aa_snp"]:
>>>>>>> e6dd8312 (Rsvg workflow main (#420))
            continue

        mmap = {int(k): v for k, v in metadata_map[col].items()}
        df[col] = df[col].map(mmap)

    df.to_csv(df_out)
