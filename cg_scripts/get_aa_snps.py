import json
import numpy as np
import pandas as pd

from cg_scripts.fasta import read_fasta_file
from cg_scripts.util import translate


def get_aa_snps(dna_snp_file, gene_or_protein_file, reference_file, mode="gene"):
    # Load the reference sequence
    with open(reference_file, "r") as fp:
        lines = fp.readlines()
        ref = read_fasta_file(lines)
        ref_seq = list(ref.values())[0]

    # JSON to dataframe
    with open(gene_or_protein_file) as fp:
        gene_or_protein_df = json.loads(fp.read())
        gene_or_protein_df = pd.DataFrame(gene_or_protein_df)

    if mode == "gene":
        # Only take protein-coding genes
        gene_or_protein_df = (
            gene_or_protein_df.loc[gene_or_protein_df["protein_coding"] == 1, :]
            # set the gene as the index
            .set_index("gene")
        )
    else:
        gene_or_protein_df = gene_or_protein_df.set_index("protein")

    dna_snp_df = pd.read_csv(dna_snp_file).fillna("")
    # Filter out any big SNPs in the 5' or 3' UTR
    dna_snp_df = dna_snp_df.loc[
        (dna_snp_df["pos"] < 29675) & (dna_snp_df["pos"] > 265), :
    ].reset_index(drop=True)
    # Filter out any frameshifting indels
    dna_snp_df = dna_snp_df.loc[
        ((dna_snp_df["ref"].str.len() == 1) & (dna_snp_df["alt"].str.len() == 1))
        | (
            (dna_snp_df["ref"].str.len() > 1)
            & (dna_snp_df["alt"].str.len() == 0)
            & (dna_snp_df["ref"].str.len() % 3 == 0)
        )
        | (
            (dna_snp_df["alt"].str.len() > 1)
            & (dna_snp_df["ref"].str.len() == 0)
            & (dna_snp_df["alt"].str.len() % 3 == 0)
        )
    ].reset_index(drop=True)

    aa_snps = []
    aa_seqs = {}

    for ref_name, ref_row in gene_or_protein_df.iterrows():
        # print(ref_name)

        segments = ref_row["segments"].split(";")

        resi_counter = 0
        aa_seqs[ref_name] = []

        for segment in segments:
            # Get the region in coordinates to translate/look for SNPs in
            segment_start = int(segment.split("..")[0])
            segment_end = int(segment.split("..")[1])

            # Translate the sequence and store it for later
            aa_seqs[ref_name] += list(
                translate(ref_seq[segment_start - 1 : segment_end])
            )

            # Get all NT SNPs in this segment
            segment_snp_df = dna_snp_df.loc[
                (dna_snp_df["pos"] >= segment_start)
                & (dna_snp_df["pos"] <= segment_end),
                :,
            ].copy()

            # For each NT SNP in this segment:
            for i, snp in segment_snp_df.iterrows():

                # Get the affected region, in codon-indexes
                # (Relative to the segment start)
                codon_ind_start = (snp["pos"] - segment_start) // 3
                codon_ind_end = (
                    snp["pos"]
                    + (0 if len(snp["ref"]) == 0 else len(snp["ref"]) - 1)
                    - segment_start
                ) // 3
                # print(codon_ind_start, codon_ind_end)

                # Get region start/end, 0-indexed
                region_start = segment_start + (codon_ind_start * 3) - 1
                region_end = segment_start + (codon_ind_end * 3) + 2
                # Position of the SNP inside the region (0-indexed)
                pos_inside_region = snp["pos"] - region_start - 1

                # Get the reference sequence of the region
                region_seq = list(ref_seq[region_start:region_end])
                # Translate the reference region sequence
                ref_aa = list(translate("".join(region_seq)))
                # print(region_seq, ref_aa)

                # Make sure the reference matches
                # print(region_seq[pos_inside_region:(pos_inside_region + len(snp['ref']))])
                if len(snp["ref"]) > 0:
                    ref_snp_seq = "".join(
                        region_seq[
                            pos_inside_region : (pos_inside_region + len(snp["ref"]))
                        ]
                    )
                    if not ref_snp_seq == snp["ref"]:
                        print(
                            "REF MISMATCH:\n\tReference sequence:\t{}\n\tSNP sequence\t\t{}\n".format(
                                ref_snp_seq, snp["ref"],
                            )
                        )
                        # I guess just move on.....

                # Remove the reference base(s)
                if len(snp["ref"]) > 0:
                    region_seq = (
                        region_seq[:pos_inside_region]
                        + region_seq[(pos_inside_region + len(snp["ref"])) :]
                    )
                # Add the alt base(s)
                if len(snp["alt"]) > 0:
                    for base in list(snp["alt"])[::-1]:
                        region_seq.insert(pos_inside_region, base)

                # Translate the new region
                alt_aa = list(translate("".join(region_seq)))
                # print(region_seq, alt_aa)

                # Remove matching AAs from the start of ref_aa
                # i.e., if ref = 'FF', and alt = 'F', then:
                #       ref = 'F' and alt = ''
                remove_inds = []
                for b in range(max(len(ref_aa), len(alt_aa))):
                    if b >= len(ref_aa) or b >= len(alt_aa):
                        break

                    if ref_aa[b] == alt_aa[b]:
                        remove_inds.append(b)
                    else:
                        break

                for ind in remove_inds[::-1]:
                    ref_aa.pop(ind)
                    alt_aa.pop(ind)
                # print(ref_aa, alt_aa)

                # If there's no mutation, or synonymous mutation,
                # then move on
                if not ref_aa and not alt_aa:
                    continue

                aa_snps.append(
                    (
                        snp["taxon"],
                        ref_name,
                        resi_counter + codon_ind_start + 1,
                        "".join(ref_aa),
                        "".join(alt_aa),
                    )
                )
            # END FOR TAXON

            resi_counter += (segment_end - segment_start + 1) // 3
        # END FOR SEGMENT
    # END FOR GENE/PROTEIN

    if mode == "gene":
        aa_snp_df = pd.DataFrame.from_records(
            aa_snps, columns=["taxon", "gene", "pos", "ref", "alt"]
        )
    else:
        aa_snp_df = pd.DataFrame.from_records(
            aa_snps, columns=["taxon", "protein", "pos", "ref", "alt"]
        )

    return aa_snp_df
