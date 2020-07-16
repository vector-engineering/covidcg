#!/usr/bin/env python3
# coding: utf-8

"""Load SNP data, create SNP signatures

Author: Albert Chen - Deverman Lab, Broad Institute
"""

import numpy as np
import pandas as pd

from reference import ref_seq
from util import data_dir, static_data_dir, translate


def load_dna_snps():
    # Load all DNA SNP files
    dna_snp_files = sorted((data_dir / "dna_snp").glob("*.csv"))
    print("Loading {} DNA SNP files...".format(len(dna_snp_files)), end="", flush=True)
    # Load into dataframe
    dna_snp_df = pd.DataFrame()
    for f in dna_snp_files:
        dna_snp_df = pd.concat([dna_snp_df, pd.read_csv(f)], ignore_index=True)
    # Extract the GISAID ID
    dna_snp_df["gisaid_id"] = dna_snp_df["taxon"].str.split("|", expand=True)[1]

    # Fill NaN values
    dna_snp_df["ref"].fillna("", inplace=True)
    dna_snp_df["alt"].fillna("", inplace=True)

    # Drop duplicate entries
    dna_snp_df.drop_duplicates(["taxon", "pos"], inplace=True)

    dna_snp_df = dna_snp_df.reset_index(drop=True)
    print("done. Loaded {} DNA SNPs".format(len(dna_snp_df)), flush=True)
    # dna_snp_df.to_csv('dna_snp.csv', index=False)

    return dna_snp_df


def process_gene_aa_snps(dna_snp_df):
    genes_df = pd.read_csv(static_data_dir / "genes.csv")

    # Only take protein-coding genes
    genes_df = (
        genes_df.loc[genes_df["protein_coding"] == 1, :]
        # set the gene as the index
        .set_index("gene")
    )

    gene_aa_snps = []
    gene_aa_seqs = {}
    for gene_name, gene in genes_df.iterrows():
        print(gene_name)

        gene_start = int(gene["start"])
        gene_end = int(gene["end"])

        gene_aa_seqs[gene_name] = translate(ref_seq[gene_start - 1 : gene_end])

        gene_snp_df = dna_snp_df.loc[
            (dna_snp_df["pos"] >= gene_start) & (dna_snp_df["pos"] <= gene_end), :
        ].copy()
        gene_snp_df["codon_ind"] = (gene_snp_df["pos"] - gene_start) // 3
        gene_snp_df["codon_start"] = gene_start + (gene_snp_df["codon_ind"] * 3)
        # spike_snps['codon_end'] = spike_snps['codon_start'] + 3
        gene_snp_df["codon_pos"] = gene_snp_df["pos"] - gene_snp_df["codon_start"]

        gene_snp_df = gene_snp_df.groupby(["taxon", "codon_ind"], as_index=False)[
            ["codon_pos", "ref", "alt"]
        ].agg(list)

        gene_snp_df["codon"] = gene_snp_df["codon_ind"].apply(
            lambda x: list(
                ref_seq[(gene_start + (x * 3) - 1) : (gene_start + (x * 3) + 2)]
            )
        )
        gene_snp_df["ref_aa"] = gene_snp_df["codon"].apply(
            lambda x: translate("".join(x))
        )

        for i, row in gene_snp_df.iterrows():
            codon = row["codon"]
            for codon_pos, ref, alt in zip(row["codon_pos"], row["ref"], row["alt"]):
                # assert codon[codon_pos] == ref
                if codon[codon_pos] != ref:
                    print(row)
                codon[codon_pos] = alt

            alt_aa = translate("".join(codon))

            if row["ref_aa"] != alt_aa:
                gene_aa_snps.append(
                    (
                        row["taxon"],
                        gene_name,
                        row["codon_ind"] + 1,
                        row["ref_aa"],
                        alt_aa,
                    )
                )

    gene_aa_snp_df = pd.DataFrame.from_records(
        gene_aa_snps, columns=["taxon", "gene", "pos", "ref", "alt"]
    )
    return gene_aa_snp_df


def process_protein_aa_snps(dna_snp_df):
    # Load proteins data
    proteins_df = pd.read_csv(static_data_dir / "proteins.csv", comment="#")
    proteins_df = proteins_df.set_index("protein")

    protein_aa_seq = {}
    protein_aa_snps = []

    for protein_name, protein in proteins_df.iterrows():
        print(protein_name)

        segments = protein["segments"].split(";")

        resi_counter = 0
        protein_aa_seq[protein_name] = []
        for segment in segments:
            segment_start = int(segment.split("..")[0])
            segment_end = int(segment.split("..")[1])
            resi_counter += (segment_end - segment_start + 1) // 3

            protein_aa_seq[protein_name] += list(
                translate(ref_seq[segment_start - 1 : segment_end])
            )

            segment_snp_df = dna_snp_df.loc[
                (dna_snp_df["pos"] >= segment_start)
                & (dna_snp_df["pos"] <= segment_end),
                :,
            ].copy()
            segment_snp_df["codon_ind"] = (segment_snp_df["pos"] - segment_start) // 3
            segment_snp_df["codon_start"] = segment_start + (
                segment_snp_df["codon_ind"] * 3
            )
            # spike_snps['codon_end'] = spike_snps['codon_start'] + 3
            segment_snp_df["codon_pos"] = (
                segment_snp_df["pos"] - segment_snp_df["codon_start"]
            )

            segment_snp_df = segment_snp_df.groupby(
                ["taxon", "codon_ind"], as_index=False
            )[["codon_pos", "ref", "alt"]].agg(list)

            segment_snp_df["codon"] = segment_snp_df["codon_ind"].apply(
                lambda x: list(
                    ref_seq[
                        (segment_start + (x * 3) - 1) : (segment_start + (x * 3) + 2)
                    ]
                )
            )
            segment_snp_df["ref_aa"] = segment_snp_df["codon"].apply(
                lambda x: translate("".join(x))
            )

            for i, row in segment_snp_df.iterrows():
                codon = row["codon"]
                for codon_pos, ref, alt in zip(
                    row["codon_pos"], row["ref"], row["alt"]
                ):
                    assert codon[codon_pos] == ref
                    codon[codon_pos] = alt

                alt_aa = translate("".join(codon))

                if row["ref_aa"] != alt_aa:
                    protein_aa_snps.append(
                        (
                            row["taxon"],
                            protein_name,
                            row["codon_ind"] + 1,
                            row["ref_aa"],
                            alt_aa,
                        )
                    )

        # print("".join(protein_aa_seq[protein_name]))

    protein_aa_snp_df = pd.DataFrame.from_records(
        protein_aa_snps, columns=["taxon", "protein", "pos", "ref", "alt"]
    )

    return protein_aa_snp_df


def generate_snp_signatures(dna_snp_df, gene_aa_snp_df, protein_aa_snp_df):
    # SNPs must occur at least this many times to pass filters
    count_threshold = 3
    # SNP signatures must be present in at least this many sequences to pass filters
    sig_count_threshold = 1

    # Collapse by taxon and count occurrences
    print("Collapsing by taxon and counting occurrences...", end="", flush=True)
    gene_aa_snp_count_df = (
        gene_aa_snp_df.groupby(["gene", "pos", "ref", "alt"], as_index=False)
        .count()
        .rename(columns={"taxon": "count"})
    )
    protein_aa_snp_count_df = (
        protein_aa_snp_df.groupby(["protein", "pos", "ref", "alt"], as_index=False)
        .count()
        .rename(columns={"taxon": "count"})
    )
    dna_snp_count_df = (
        dna_snp_df.groupby(["pos", "ref", "alt"], as_index=False)
        .count()
        .rename(columns={"taxon": "count"})
    )
    print("done")

    # Filter out SNPs
    print("Filtering low-occurrence SNPs...", end="", flush=True)
    valid_gene_aa_snps = (
        gene_aa_snp_count_df.loc[gene_aa_snp_count_df["count"] >= count_threshold, :]
        .reset_index(drop=True)
        .sort_values("count", ascending=False)
    )
    valid_protein_aa_snps = (
        protein_aa_snp_count_df.loc[
            protein_aa_snp_count_df["count"] >= count_threshold, :
        ]
        .reset_index(drop=True)
        .sort_values("count", ascending=False)
    )
    valid_dna_snps = (
        dna_snp_count_df.loc[dna_snp_count_df["count"] >= count_threshold, :]
        .reset_index(drop=True)
        .sort_values("count", ascending=False)
    )
    print("done")

    # Create unique SNP string
    print("Creating unique SNP strings...", end="", flush=True)
    valid_gene_aa_snps["snp_str"] = (
        valid_gene_aa_snps[["gene", "pos", "ref", "alt"]]
        .applymap(str)
        .agg("|".join, axis=1)
    )
    valid_protein_aa_snps["snp_str"] = (
        valid_protein_aa_snps[["protein", "pos", "ref", "alt"]]
        .applymap(str)
        .agg("|".join, axis=1)
    )
    valid_dna_snps["snp_str"] = (
        valid_dna_snps[["pos", "ref", "alt"]].applymap(str).agg("|".join, axis=1)
    )
    print("done")

    # Save to disk
    # valid_dna_snps.to_csv(data_dir / 'valid_dna_snps.csv', index=False)
    # valid_gene_aa_snps.to_csv(data_dir / 'valid_gene_aa_snps.csv', index=False)
    # valid_protein_aa_snps.to_csv(data_dir / 'valid_protein_aa_snps.csv', index=False)

    # Generate SNP strings for initial dataframes
    print("Generating SNP strings for initial dataframes...", end="", flush=True)
    gene_aa_snp_df["snp_str"] = gene_aa_snp_df["gene"].str.cat(
        [
            gene_aa_snp_df["pos"].astype(str),
            gene_aa_snp_df["ref"],
            gene_aa_snp_df["alt"],
        ],
        sep="|",
    )
    protein_aa_snp_df["snp_str"] = protein_aa_snp_df["protein"].str.cat(
        [
            protein_aa_snp_df["pos"].astype(str),
            protein_aa_snp_df["ref"],
            protein_aa_snp_df["alt"],
        ],
        sep="|",
    )
    dna_snp_df["snp_str"] = (
        dna_snp_df["pos"]
        .astype(str)
        .str.cat([dna_snp_df["ref"], dna_snp_df["alt"]], sep="|")
    )
    print("done")

    # Filter SNPs by valid SNPs
    print("Generating SNP signatures...", end="", flush=True)
    gene_aa_snp_df = gene_aa_snp_df.loc[
        gene_aa_snp_df["snp_str"].isin(valid_gene_aa_snps["snp_str"]), :
    ].reset_index(drop=True)
    protein_aa_snp_df = protein_aa_snp_df.loc[
        protein_aa_snp_df["snp_str"].isin(valid_protein_aa_snps["snp_str"]), :
    ].reset_index(drop=True)
    dna_snp_df = dna_snp_df.loc[
        dna_snp_df["snp_str"].isin(valid_dna_snps["snp_str"]), :
    ].reset_index(drop=True)

    # Group by taxon and make a ';' delimited list of snp_strs
    gene_aa_snp_group_df = (
        gene_aa_snp_df.groupby("taxon")["snp_str"].agg(";".join).reset_index()
    )
    protein_aa_snp_group_df = (
        protein_aa_snp_df.groupby("taxon")["snp_str"].agg(";".join).reset_index()
    )
    dna_snp_group_df = (
        dna_snp_df.groupby("taxon")["snp_str"].agg(";".join).reset_index()
    )

    # Count occurences of SNP signatures
    # gene_aa_snp_sig_count_df = (
    #     gene_aa_snp_group_df.groupby("snp_str").count().sort_values("taxon", ascending=False)
    # )
    # protein_aa_snp_sig_count_df = (
    #     protein_aa_snp_group_df.groupby("snp_str").count().sort_values("taxon", ascending=False)
    # )
    # dna_snp_sig_count_df = (
    #     dna_snp_group_df.groupby("snp_str")
    #     .count()
    #     .sort_values("taxon", ascending=False)
    # )

    # Only take SNP signatures with >= N sequences
    # valid_gene_aa_snp_sigs = gene_aa_snp_sig_count_df.index.values[
    #     gene_aa_snp_sig_count_df["taxon"] >= sig_count_threshold
    # ]
    # valid_protein_aa_snp_sigs = protein_aa_snp_sig_count_df.index.values[
    #     protein_aa_snp_sig_count_df["taxon"] >= sig_count_threshold
    # ]
    # valid_dna_snp_sigs = dna_snp_sig_count_df.index.values[
    #     dna_snp_sig_count_df["taxon"] >= sig_count_threshold
    # ]

    # Expand string form for easier searching
    # valid_gene_aa_snp_sigs_explode = [set(s.split(";")) for s in valid_gene_aa_snp_sigs]
    # valid_protein_aa_snp_sigs_explode = [set(s.split(";")) for s in valid_protein_aa_snp_sigs]
    # valid_dna_snp_sigs_explode = [set(s.split(";")) for s in valid_dna_snp_sigs]
    print("done")

    # Assign taxons that don't already fit perfectly into a group, to a group

    # gene_aa_snp_group_df["snp_sig"] = None
    # protein_aa_snp_group_df["snp_sig"] = None
    # dna_snp_group_df["snp_sig"] = None
    # gene_aa_snp_group_df["snp_sig"] = gene_aa_snp_group_df["snp_str"]
    # protein_aa_snp_group_df["snp_sig"] = protein_aa_snp_group_df["snp_str"]
    # dna_snp_group_df["snp_sig"] = dna_snp_group_df["snp_str"]

    # Save to disk
    dna_snp_group_df.to_csv(data_dir / "dna_snp_group.csv", index=False)
    gene_aa_snp_group_df.to_csv(data_dir / "gene_aa_snp_group.csv", index=False)
    protein_aa_snp_group_df.to_csv(data_dir / "protein_aa_snp_group.csv", index=False)

    return dna_snp_group_df, gene_aa_snp_group_df, protein_aa_snp_group_df


def process_snp_data():

    dna_snp_df = load_dna_snps()

    gene_aa_snp_df = pd.read_csv(data_dir / "gene_aa_snp.csv")
    protein_aa_snp_df = pd.read_csv(data_dir / "protein_aa_snp.csv")

    # Group SNPs by taxon, find and assign SNP signatures
    (
        dna_snp_group_df,
        gene_aa_snp_group_df,
        protein_aa_snp_group_df,
    ) = generate_snp_signatures(dna_snp_df, gene_aa_snp_df, protein_aa_snp_df)

    # Extract the GISAID id from the taxon column
    dna_snp_group_df["Accession ID"] = (
        dna_snp_group_df["taxon"].str.split("|").apply(lambda x: x[1])
    )
    gene_aa_snp_group_df["Accession ID"] = (
        gene_aa_snp_group_df["taxon"].str.split("|").apply(lambda x: x[1])
    )
    protein_aa_snp_group_df["Accession ID"] = (
        protein_aa_snp_group_df["taxon"].str.split("|").apply(lambda x: x[1])
    )

    # Map SNPs to integer IDs
    print("Mapping SNPs to integers...", end="", flush=True)
    dna_snp_map = pd.Series(
        np.unique(
            np.concatenate(dna_snp_group_df["snp_str"].str.split(";").values).ravel()
        )
    )
    gene_aa_snp_map = pd.Series(
        np.unique(
            np.concatenate(
                gene_aa_snp_group_df["snp_str"].str.split(";").values
            ).ravel()
        )
    )
    protein_aa_snp_map = pd.Series(
        np.unique(
            np.concatenate(
                protein_aa_snp_group_df["snp_str"].str.split(";").values
            ).ravel()
        )
    )

    # Flip index and values
    dna_snp_map = pd.Series(dna_snp_map.index.values, index=dna_snp_map.values)
    gene_aa_snp_map = pd.Series(
        gene_aa_snp_map.index.values, index=gene_aa_snp_map.values
    )
    protein_aa_snp_map = pd.Series(
        protein_aa_snp_map.index.values, index=protein_aa_snp_map.values
    )

    # Save maps
    dna_snp_map.to_csv(data_dir / "dna_snp_map.csv", index_label="snp", header=["id"])
    dna_snp_map.to_json(data_dir / "dna_snp_map.json", orient="index")
    gene_aa_snp_map.to_csv(
        data_dir / "gene_aa_snp_map.csv", index_label="snp", header=["id"]
    )
    gene_aa_snp_map.to_json(data_dir / "gene_aa_snp_map.json", orient="index")
    protein_aa_snp_map.to_csv(
        data_dir / "protein_aa_snp_map.csv", index_label="snp", header=["id"]
    )
    protein_aa_snp_map.to_json(data_dir / "protein_aa_snp_map.json", orient="index")

    # Convert SNP strings to integer lists
    dna_snp_group_df["snp_str"] = (
        dna_snp_group_df["snp_str"]
        .str.split(";")
        .apply(lambda x: ";".join([str(dna_snp_map[a]) for a in x] if x else None))
    )
    # dna_snp_group_df["snp_sig"] = (
    #     dna_snp_group_df["snp_sig"]
    #     .str.split(";")
    #     .apply(lambda x: ";".join([str(dna_snp_map[a]) for a in x] if x else None))
    # )
    gene_aa_snp_group_df["snp_str"] = (
        gene_aa_snp_group_df["snp_str"]
        .str.split(";")
        .apply(lambda x: ";".join([str(gene_aa_snp_map[a]) for a in x] if x else None))
    )
    protein_aa_snp_group_df["snp_str"] = (
        protein_aa_snp_group_df["snp_str"]
        .str.split(";")
        .apply(
            lambda x: ";".join([str(protein_aa_snp_map[a]) for a in x] if x else None)
        )
    )
    # gene_aa_snp_group_df["snp_sig"] = (
    #     gene_aa_snp_group_df["snp_sig"]
    #     .str.split(";")
    #     .apply(lambda x: ";".join([str(gene_aa_snp_map[a]) for a in x] if x else None))
    # )
    # protein_aa_snp_group_df["snp_sig"] = (
    #     protein_aa_snp_group_df["snp_sig"]
    #     .str.split(";")
    #     .apply(lambda x: ";".join([str(protein_aa_snp_map[a]) for a in x] if x else None))
    # )
    print("done. Saved SNP -> integer maps")

    return dna_snp_group_df, gene_aa_snp_group_df, protein_aa_snp_group_df


def main():
    
    dna_snp_df = load_dna_snps()
    # print(dna_snp_df)

    # Filter out indels and SNPs with length > 1
    # Need to figure out what to do with those...
    dna_snp_df = (
        dna_snp_df.fillna("")
        .loc[(dna_snp_df["ref"].str.len() == 1) & (dna_snp_df["alt"].str.len() == 1), :]
        .reset_index(drop=True)
    )

    gene_aa_snp_df = process_gene_aa_snps(dna_snp_df)
    # print(gene_aa_snp_df)
    protein_aa_snp_df = process_protein_aa_snps(dna_snp_df)
    # print(protein_aa_snp_df)

    gene_aa_snp_df.to_csv(data_dir / "gene_aa_snp.csv", index=False)
    protein_aa_snp_df.to_csv(data_dir / "protein_aa_snp.csv", index=False)
    
    dna_snp_group_df, gene_aa_snp_group_df, protein_aa_snp_group_df = process_snp_data()
    print(dna_snp_group_df)
    print(gene_aa_snp_group_df)
    print(protein_aa_snp_group_df)

if __name__ == "__main__":
    main()
