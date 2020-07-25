import pandas as pd
import numpy as np

from cg_scripts.get_dna_snps import load_dna_snp_file
from cg_scripts.reference import ref_seq
from cg_scripts.util import translate


def process_gene_aa_snps(dna_snp_df, genes_file):
    genes_df = pd.read_csv(genes_file)

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


def process_protein_aa_snps(dna_snp_df, proteins_file):
    # Load proteins data
    proteins_df = pd.read_csv(proteins_file)
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


def get_aa_snps(
    dna_snp_file, genes_file, proteins_file, gene_aa_snp_file, protein_aa_snp_file
):
    dna_snp_df = load_dna_snp_file(dna_snp_file)

    # Filter out indels and SNPs with length > 1
    # Need to figure out what to do with those...
    # dna_snp_df = (
    #     dna_snp_df.fillna("")
    #     .loc[(dna_snp_df["ref"].str.len() == 1) & (dna_snp_df["alt"].str.len() == 1), :]
    #     .reset_index(drop=True)
    # )

    gene_aa_snp_df = process_gene_aa_snps(dna_snp_df, genes_file)
    # print(gene_aa_snp_df)
    protein_aa_snp_df = process_protein_aa_snps(dna_snp_df, proteins_file)
    # print(protein_aa_snp_df)

    gene_aa_snp_df.to_csv(gene_aa_snp_file, index=False)
    protein_aa_snp_df.to_csv(protein_aa_snp_file, index=False)
