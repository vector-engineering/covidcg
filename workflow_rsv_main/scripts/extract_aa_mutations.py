# coding: utf-8

"""Extract AA mutations from NT mutations

Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""


import pandas as pd
import json

from scripts.fasta import read_fasta_file
from scripts.util import translate


def extract_aa_mutations(
    dna_mutation_file, gene_or_protein_file, reference_file, mode="gene"
):
    # Load the reference sequence
    with open(reference_file, "r") as fp:
        lines = fp.readlines()
        ref = read_fasta_file(lines)

    # JSON to dataframe
    with open(gene_or_protein_file, "r") as f:
        data = json.load(f)

    dna_mutation_df = pd.read_csv(dna_mutation_file).fillna("")
    # Filter out any frameshifting indels
    dna_mutation_df = dna_mutation_df.loc[
        (
            (dna_mutation_df["ref"].str.len() == 1)
            & (dna_mutation_df["alt"].str.len() == 1)
        )
        | (
            (dna_mutation_df["ref"].str.len() > 1)
            & (dna_mutation_df["alt"].str.len() == 0)
            & (dna_mutation_df["ref"].str.len() % 3 == 0)
        )
        | (
            (dna_mutation_df["alt"].str.len() > 1)
            & (dna_mutation_df["ref"].str.len() == 0)
            & (dna_mutation_df["alt"].str.len() % 3 == 0)
        )
    ].reset_index(drop=True)

    aa_mutations = []
    aa_seqs = {}
    for reference in ref:
        ref_seq = ref[reference]
        gene_or_protein_df = pd.read_json(json.dumps(data[reference]))
        if mode == "gene":
            # Only take protein-coding genes
            gene_or_protein_df = (
                gene_or_protein_df.loc[gene_or_protein_df["protein_coding"] == 1, :]
                # set the gene as the index
                .set_index("name")
            )
        else:
            gene_or_protein_df = gene_or_protein_df.set_index("name")

        if reference == "NC_038235.1":
            subtype = "A"
        else:
            subtype = "B"

        for ref_name, ref_row in gene_or_protein_df.iterrows():
            # print(ref_name)

            resi_counter = 0
            aa_seqs[ref_name] = []

            for segment in ref_row["segments"]:
                # Get the region in coordinates to translate/look for mutations in
                segment_start = segment[0]
                segment_end = segment[1]
                segment_len = ((segment_end - segment_start) // 3) + 1

                # Get all NT mutations in this segment
                segment_mutation_df = (
                    dna_mutation_df.loc[
                        (dna_mutation_df["pos"] >= segment_start)
                        & (dna_mutation_df["pos"] <= segment_end)
                        & (dna_mutation_df["subtype"] == subtype),
                        :,
                    ]
                    .copy()
                    .sort_values(["Accession ID", "pos"])
                )

                segment_mutation_df["codon_ind_start"] = (
                    segment_mutation_df["pos"] - segment_start
                ) // 3
                segment_mutation_df["codon_ind_end"] = (
                    segment_mutation_df["pos"]
                    + segment_mutation_df["ref"].apply(
                        lambda x: 0 if len(x) == 0 else len(x) - 1
                    )
                    - segment_start
                ) // 3

                # For each NT mutation in this segment:
                i = 0
                while i < len(segment_mutation_df):
                    cur_mutation = segment_mutation_df.iloc[i, :]

                    # Translate the sequence and store it for later
                    aa_seqs[ref_name] += list(
                        translate(ref_seq[(segment_start - 1):segment_end])
                    )

                    # Look ahead in the mutation list (ordered by position)
                    # for any other mutations within this codon
                    j = i + 1
                    mutations = [cur_mutation]
                    codon_ind_start = cur_mutation["codon_ind_start"]
                    codon_ind_end = cur_mutation["codon_ind_end"]

                    if j < len(segment_mutation_df):
                        new_mutation = segment_mutation_df.iloc[j]
                        while (
                            new_mutation["codon_ind_start"] <= codon_ind_end
                            and new_mutation["Accession ID"] == cur_mutation["Accession ID"]
                        ):
                            mutations.append(new_mutation)
                            # Update end position
                            codon_ind_end = new_mutation["codon_ind_end"]
                            j += 1
                            if j < len(segment_mutation_df):
                                new_mutation = segment_mutation_df.iloc[j]
                            else:
                                break

                    # Increment our counter so we don't reprocess any grouped mutations
                    i = j

                    # Get region start/end, 0-indexed
                    region_start = segment_start + (codon_ind_start * 3) - 1
                    region_end = segment_start + (codon_ind_end * 3) + 2

                    # Position of the mutation inside the region (0-indexed)
                    pos_inside_region = []
                    for mut in mutations:
                        pos_inside_region.append(mut["pos"] - region_start - 1)

                    # Get the reference sequence of the region
                    region_seq = list(ref_seq[region_start:region_end])
                    initial_region_seq = [s for s in region_seq]
                    # Translate the reference region sequence
                    ref_aa = list(translate("".join(region_seq)))
                    # print(region_seq, ref_aa)

                    # print([mut['pos'] for mut in mutations])
                    # print(region_seq, ref_aa)

                    for k, mut in enumerate(mutations):

                        # Make sure the reference matches
                        if len(mut["ref"]) > 0:
                            ref_mutation_seq = "".join(
                                initial_region_seq[
                                    pos_inside_region[k] : (
                                        pos_inside_region[k] + len(mut["ref"])
                                    )
                                ]
                            )
                            if not ref_mutation_seq == mut["ref"]:
                                print(
                                    "REF MISMATCH:\n\tReference sequence:\t{}\n\tMutation sequence\t\t{}\n".format(
                                        ref_mutation_seq, mut["ref"],
                                    )
                                )
                                continue

                        # Remove the reference base(s)
                        if len(mut["ref"]) > 0:
                            region_seq = (
                                region_seq[: pos_inside_region[k]]
                                + region_seq[(pos_inside_region[k] + len(mut["ref"])) :]
                            )
                        # Add the alt base(s)
                        if len(mut["alt"]) > 0:
                            for base in list(mut["alt"])[::-1]:
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

                    aa_mutations.append(
                        (
                            cur_mutation["Accession ID"],
                            ref_name,
                            pos,
                            "".join(ref_aa),
                            "".join(alt_aa),
                            subtype
                        )
                    )
                # END FOR ACCESSION ID

                resi_counter += segment_len
            # END FOR SEGMENT
        # END FOR GENE/PROTEIN
    # END FOR REFERENCE

    if mode == "gene":
        aa_mutation_df = pd.DataFrame.from_records(
            aa_mutations, columns=["Accession ID", "gene", "pos", "ref", "alt", "subtype"]
        )
    else:
        aa_mutation_df = pd.DataFrame.from_records(
            aa_mutations, columns=["Accession ID", "protein", "pos", "ref", "alt", "subtype"]
        )
        
    return aa_mutation_df
