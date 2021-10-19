# coding: utf-8

"""Extract AA SNPs from NT SNPs

Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""


import json
import numpy as np
import pandas as pd

from scripts.fasta import read_fasta_file
from scripts.util import translate


def extract_aa_snps(dna_snp_file, gene_or_protein_file, reference_file, mode="gene"):
    # Load the reference sequence
    with open(reference_file, "r") as fp:
        lines = fp.readlines()
        ref = read_fasta_file(lines)
        ref_seq = list(ref.values())[0]

    # JSON to dataframe
    gene_or_protein_df = pd.read_json(gene_or_protein_file)

    if mode == "gene":
        # Only take protein-coding genes
        gene_or_protein_df = (
            gene_or_protein_df.loc[gene_or_protein_df["protein_coding"] == 1, :]
            # set the gene as the index
            .set_index("name")
        )
    else:
        gene_or_protein_df = gene_or_protein_df.set_index("name")

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

        resi_counter = 0
        aa_seqs[ref_name] = []

        for segment in ref_row["segments"]:
            # Get the region in coordinates to translate/look for SNPs in
            segment_start = segment[0]
            segment_end = segment[1]
            segment_len = ((segment_end - segment_start) // 3) + 1

            # Translate the sequence and store it for later
            aa_seqs[ref_name] += list(
                translate(ref_seq[segment_start - 1 : segment_end])
            )

            # Get all NT SNPs in this segment
            segment_snp_df = (
                dna_snp_df.loc[
                    (dna_snp_df["pos"] >= segment_start)
                    & (dna_snp_df["pos"] <= segment_end),
                    :,
                ]
                .copy()
                .sort_values(["Accession ID", "pos"])
            )

            segment_snp_df["codon_ind_start"] = (
                segment_snp_df["pos"] - segment_start
            ) // 3
            segment_snp_df["codon_ind_end"] = (
                segment_snp_df["pos"]
                + segment_snp_df["ref"].apply(
                    lambda x: 0 if len(x) == 0 else len(x) - 1
                )
                - segment_start
            ) // 3

            # For each NT SNP in this segment:
            i = 0
            while i < len(segment_snp_df):
                cur_snp = segment_snp_df.iloc[i, :]

                # Look ahead in the SNV list (ordered by position)
                # for any other SNVs within this codon
                j = i + 1
                snps = [cur_snp]
                codon_ind_start = cur_snp["codon_ind_start"]
                codon_ind_end = cur_snp["codon_ind_end"]

                if j < len(segment_snp_df):
                    new_snp = segment_snp_df.iloc[j]
                    while (
                        new_snp["codon_ind_start"] <= codon_ind_end
                        and new_snp["Accession ID"] == cur_snp["Accession ID"]
                    ):
                        snps.append(new_snp)
                        # Update end position
                        codon_ind_end = new_snp["codon_ind_end"]
                        j += 1
                        if j < len(segment_snp_df):
                            new_snp = segment_snp_df.iloc[j]
                        else:
                            break

                # Increment our counter so we don't reprocess any grouped SNVs
                i = j

                # Get region start/end, 0-indexed
                region_start = segment_start + (codon_ind_start * 3) - 1
                region_end = segment_start + (codon_ind_end * 3) + 2

                # Position of the SNP inside the region (0-indexed)
                pos_inside_region = []
                for snp in snps:
                    pos_inside_region.append(snp["pos"] - region_start - 1)

                # Get the reference sequence of the region
                region_seq = list(ref_seq[region_start:region_end])
                initial_region_seq = [s for s in region_seq]
                # Translate the reference region sequence
                ref_aa = list(translate("".join(region_seq)))
                # print(region_seq, ref_aa)

                # print([snp['pos'] for snp in snps])
                # print(region_seq, ref_aa)

                for k, snp in enumerate(snps):

                    # Make sure the reference matches
                    if len(snp["ref"]) > 0:
                        ref_snp_seq = "".join(
                            initial_region_seq[
                                pos_inside_region[k] : (
                                    pos_inside_region[k] + len(snp["ref"])
                                )
                            ]
                        )
                        if not ref_snp_seq == snp["ref"]:
                            print(
                                "REF MISMATCH:\n\tReference sequence:\t{}\n\tSNP sequence\t\t{}\n".format(
                                    ref_snp_seq, snp["ref"],
                                )
                            )
                            continue

                    # Remove the reference base(s)
                    if len(snp["ref"]) > 0:
                        region_seq = (
                            region_seq[: pos_inside_region[k]]
                            + region_seq[(pos_inside_region[k] + len(snp["ref"])) :]
                        )
                    # Add the alt base(s)
                    if len(snp["alt"]) > 0:
                        for base in list(snp["alt"])[::-1]:
                            region_seq.insert(pos_inside_region[k], base)

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
                        # Increment the AA index start so that
                        # we end up on the correct position
                        # (This should only affect deletions)
                        codon_ind_start += 1
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

                pos = resi_counter + codon_ind_start + 1

                # If the positiion is outside of the segment, then skip
                # (This happens sometimes for long deletions)
                if pos > (resi_counter + segment_len):
                    continue

                # If the ref AA sequence overruns the segment
                # (this happens for long deletions near the end of a gene)
                # then truncate the ref AAs
                overrun = (pos + len(ref_aa) - 1) - (resi_counter + segment_len)
                if overrun > 0:
                    ref_aa = ref_aa[:-overrun]

                aa_snps.append(
                    (
                        cur_snp["Accession ID"],
                        ref_name,
                        pos,
                        "".join(ref_aa),
                        "".join(alt_aa),
                    )
                )
            # END FOR ACCESSION ID

            resi_counter += segment_len
        # END FOR SEGMENT
    # END FOR GENE/PROTEIN

    if mode == "gene":
        aa_snp_df = pd.DataFrame.from_records(
            aa_snps, columns=["Accession ID", "gene", "pos", "ref", "alt"]
        )
    else:
        aa_snp_df = pd.DataFrame.from_records(
            aa_snps, columns=["Accession ID", "protein", "pos", "ref", "alt"]
        )

    return aa_snp_df
