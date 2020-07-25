import pandas as pd
import numpy as np


def map_snps(dna_snp_file, gene_aa_snp_file, protein_aa_snp_file):

    dna_snp_df = pd.read_csv(dna_snp_file)
    gene_aa_snp_df = pd.read_csv(gene_aa_snp_file)
    protein_aa_snp_df = pd.read_csv(protein_aa_snp_file)

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

