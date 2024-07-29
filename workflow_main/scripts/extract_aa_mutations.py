#!/usr/bin/env python3
# coding: utf-8

"""Extract AA mutations from NT mutations

Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""

import argparse
import json
import pandas as pd

from util import translate

# Number of codons to look ahead for resolving frameshifts
FS_LOOKAHEAD = 10


def extract_aa_mutations(
    dna_mutation_file,
    gene_or_protein_file,
    reference_file,
    subtype,
    active_segment,
    mode="gene",
):
    """
    Extract AA mutations from NT mutations

    Parameters
    ----------
    dna_mutation_file: str
        Path to DNA mutation CSV
    gene_or_protein_file: str
        Path to gene/protein definition JSON
    reference_file: str
        Path to reference JSON file
    subtype: str
        Subtype
    active_segment: str
        Segment/chromosome
    mode: str
        Mode ('gene' or 'protein')
    """

    # Load the reference sequences
    with open(reference_file, "r") as fp:
        references = json.loads(fp.read())

    # Get references for this subtype
    references = {k: v for k, v in references.items() if v["subtype"] == subtype}

    # Get sequences for all references under this subtype,
    # but only for the active segment
    ref_seqs = {
        ref["name"]: ref["segments"][active_segment]["sequence"]
        for ref in references.values()
        if active_segment in ref["segments"]
    }

    # Load gene/protein defs
    # JSON to dataframe
    with open(gene_or_protein_file, "r") as fp:
        feature_dicts = json.loads(fp.read())

    # Get only features for the above references
    feature_dicts = {k: v for k, v in feature_dicts.items() if k in references.keys()}

    feature_dfs = {}
    for k, v in feature_dicts.items():
        v = pd.DataFrame.from_records(v)

        if mode == "gene":
            # Only take protein-coding genes
            v = (
                v.loc[v["protein_coding"] == 1, :]
                # set the gene as the index
                .set_index("name")
            )
        else:
            v = v.set_index("name")

        feature_dfs[k] = v

    dna_mutation_df = pd.read_csv(dna_mutation_file, index_col=False).fillna("")

    # Filter out any big mutations in the 5' or 3' UTR
    # dna_mutation_df = dna_mutation_df.loc[
    #     (dna_mutation_df["pos"] < 29675) & (dna_mutation_df["pos"] > 265), :
    # ].reset_index(drop=True)

    # Filter out any frameshifting indels
    # dna_mutation_df = dna_mutation_df.loc[
    #     (
    #         (dna_mutation_df["ref"].str.len() == 1)
    #         & (dna_mutation_df["alt"].str.len() == 1)
    #     )
    #     | (
    #         (dna_mutation_df["ref"].str.len() > 1)
    #         & (dna_mutation_df["alt"].str.len() == 0)
    #         & (dna_mutation_df["ref"].str.len() % 3 == 0)
    #     )
    #     | (
    #         (dna_mutation_df["alt"].str.len() > 1)
    #         & (dna_mutation_df["ref"].str.len() == 0)
    #         & (dna_mutation_df["alt"].str.len() % 3 == 0)
    #     )
    # ].reset_index(drop=True)

    aa_mutations = []
    aa_seqs = {}

    # For each reference
    for reference, feature_df in feature_dfs.items():

        # For each feature (gene/protein)
        for feature_name, feature_row in feature_df.iterrows():

            # Only process genes in this segment/chromosome
            if feature_row["segment"] != active_segment:
                continue

            resi_counter = 0
            aa_seqs[feature_name] = []

            # For each segment (ORF) in this feature:
            for segment in feature_row["segments"]:
                # Get the region in coordinates to translate/look for mutations in
                segment_start = segment[0]
                segment_end = segment[1]
                segment_len = ((segment_end - segment_start) // 3) + 1

                # Translate the sequence and store it for later
                aa_seqs[feature_name] += list(
                    translate(ref_seqs[reference][segment_start - 1 : segment_end])
                )

                # Get all NT mutations in this segment, for this reference
                segment_mutation_df = (
                    dna_mutation_df.loc[
                        (dna_mutation_df["reference"] == reference)
                        & (dna_mutation_df["pos"] >= segment_start)
                        & (dna_mutation_df["pos"] <= segment_end),
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

                # GROUP DNA MUTATIONS
                # Group together individual DNA mutations to process them
                # as a single DNA mutation
                # Criteria:
                # - Same Accession ID
                # - Mutations are within codon range

                # For each NT mutation in this segment:
                segment_mutations = []
                i = 0
                while i < len(segment_mutation_df):
                    cur_mutation = segment_mutation_df.iloc[i, :]

                    # If the current mutation results in a frameshift,
                    # then seek to resolve the frameshift by looking ahead FS_LOOKAHEAD codons
                    # i.e., if this mutation is a single deletion, then
                    # look ahead for a single insertion to resolve the frame

                    # Number of bases off from the frame
                    # Use a modulo of 3 later to get the # of bases off
                    frameshift = 0

                    if len(cur_mutation["alt"]) - len(cur_mutation["ref"]) % 3 != 0:
                        frameshift += len(cur_mutation["alt"]) - len(
                            cur_mutation["ref"]
                        )

                    # If we have a frameshift, then look ahead FS_LOOKAHEAD codons
                    # If no frameshift resolution is possible then toss this mutation and move on
                    # Keep track of when the indel resolves
                    resolve_codon_ind = None
                    if frameshift % 3 != 0 and i < len(segment_mutation_df) - 1:
                        j = i + 1
                        new_mutation = segment_mutation_df.iloc[j]
                        while (
                            new_mutation["codon_ind_start"]
                            <= cur_mutation["codon_ind_start"] + FS_LOOKAHEAD
                            and new_mutation["Accession ID"]
                            == cur_mutation["Accession ID"]
                        ):
                            # Adjust frameshift if the new mutation is itself a frameshifting indel
                            if (
                                len(new_mutation["alt"]) - len(new_mutation["ref"]) % 3
                                != 0
                            ):
                                frameshift += len(new_mutation["alt"]) - len(
                                    new_mutation["ref"]
                                )

                            # If the frameshift is resolved, then flag as resolved and break
                            if frameshift % 3 == 0:
                                resolve_codon_ind = new_mutation["codon_ind_end"]
                                break

                            j += 1
                            if j < len(segment_mutation_df):
                                new_mutation = segment_mutation_df.iloc[j]
                            else:
                                break

                        # If there's no resolution to the frameshift within FS_LOOKAHEAD codons,
                        # then ignore this mutation and move onto the next mutation in the list
                        if resolve_codon_ind is None:
                            i += 1
                            continue

                    # If this was a singular frameshift mutation without a resolution,
                    # Then no possibility to resolve the frameshift
                    # Just ignore this mutation
                    # TODO: process nonsense frameshifts if near the end of the ORF?
                    if frameshift % 3 != 0 and resolve_codon_ind is None:
                        i += 1
                        continue

                    # Look ahead in the mutation list (ordered by position)
                    # for any other mutations within the codon range
                    j = i + 1
                    mutations = [cur_mutation.to_dict()]
                    codon_ind_start = cur_mutation["codon_ind_start"]
                    codon_ind_end = cur_mutation["codon_ind_end"]

                    # If we have a resolvable frameshift, then set the mutation combination window
                    # From here til where the indel resolves
                    if resolve_codon_ind is not None:
                        codon_ind_end = resolve_codon_ind

                    if j < len(segment_mutation_df):
                        new_mutation = segment_mutation_df.iloc[j]
                        while (
                            new_mutation["codon_ind_start"] <= codon_ind_end
                            and new_mutation["Accession ID"]
                            == cur_mutation["Accession ID"]
                        ):
                            # If we weren't expecting a frameshift and this mutation
                            # produces a frameshift, then ignore this mutation
                            if resolve_codon_ind is None and (
                                len(new_mutation["alt"]) - len(new_mutation["ref"]) % 3
                                != 0
                            ):
                                pass
                            else:
                                mutations.append(new_mutation.to_dict())

                                # Update end position if beyond current end position
                                if new_mutation["codon_ind_end"] > codon_ind_end:
                                    codon_ind_end = new_mutation["codon_ind_end"]

                            j += 1
                            if j < len(segment_mutation_df):
                                new_mutation = segment_mutation_df.iloc[j]
                            else:
                                break

                    # Increment our counter so we don't reprocess any grouped mutations
                    i = j

                    if mutations:
                        segment_mutations.append(mutations)

                # PROCESS GROUPED MUTATIONS
                for mutations in segment_mutations:

                    codon_ind_start = mutations[0]["codon_ind_start"]
                    codon_ind_end = mutations[-1]["codon_ind_end"]
                    cur_mutation = mutations[0]

                    # Get region start/end, 0-indexed
                    region_start = segment_start + (codon_ind_start * 3) - 1
                    region_end = segment_start + (codon_ind_end * 3) + 2

                    # Position of the mutation inside the region (0-indexed)
                    pos_inside_region = []
                    for mut in mutations:
                        pos_inside_region.append(mut["pos"] - region_start - 1)

                    # Get the reference sequence of the region
                    region_seq = list(ref_seqs[reference][region_start:region_end])
                    initial_region_seq = [s for s in region_seq]
                    # Translate the reference region sequence
                    ref_aa = list(translate("".join(region_seq)))

                    # Process mutations in reverse order to preserve indexing
                    for k, mut in zip(
                        range(len(mutations) - 1, -1, -1), mutations[::-1]
                    ):
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
                                    "REF MISMATCH:\n\tReference sequence:\t{}\n\tMutation sequence:\t{}\n".format(
                                        ref_mutation_seq,
                                        mut["ref"],
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

                    # If synonymous mutation, then move on:
                    if ref_aa == alt_aa:
                        continue

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

                    # Remove matching AAs from the end of ref_aa
                    # i.e., if ref = 'SRG', and alt = 'KG', then:
                    #          ref = 'SR' and alt = 'K'
                    remove_inds = []
                    for b in range(max(len(ref_aa), len(alt_aa))):
                        if b >= len(ref_aa) or b >= len(alt_aa):
                            break

                        if ref_aa[-1 - b] == alt_aa[-1 - b]:
                            remove_inds.append(b)
                            # Increment the AA index end so that
                            # we end up on the correct position
                            # (This should only affect deletions)
                            codon_ind_end -= 1
                        else:
                            break

                    for ind in remove_inds[::-1]:
                        ref_aa.pop(-1 - ind)
                        alt_aa.pop(-1 - ind)

                    # If there's no mutation, or synonymous mutation,
                    # then move on
                    if not ref_aa and not alt_aa:
                        continue

                    pos = resi_counter + codon_ind_start + 1

                    # If the position is outside of the segment, then skip
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
                            cur_mutation["reference"],
                            cur_mutation["Accession ID"],
                            feature_name,
                            # Apply feature residue offset here
                            # for renumbering residues off of, e.g., signal peptide
                            pos - feature_row["residue_offset"],
                            "".join(ref_aa),
                            "".join(alt_aa),
                        )
                    )
                # END FOR ACCESSION ID

                resi_counter += segment_len
            # END FOR SEGMENT
        # END FOR GENE/PROTEIN

    # Combine AA mutations
    # Merge AA mutations only if:
    # - Reference is the same
    # - Accession ID is the same
    # - gene/protein is the same
    # - position is the same
    # - one mutation is an indel, other is a substitution

    """
    if len(aa_mutations) > 0:
        cur_mutation = aa_mutations[0]
        i = 1
        while i < len(aa_mutations):

            next_mutation = aa_mutations[i]

            # (reference, Accession ID, gene/protein, pos, ref, alt)
            if (
                cur_mutation[0] == next_mutation[0]
                and cur_mutation[1] == next_mutation[1]
                and cur_mutation[2] == next_mutation[2]
                and cur_mutation[3] == next_mutation[3]
                and (
                    (cur_mutation[4] == "" or cur_mutation[5] == "")
                    or (next_mutation[4] == "" or next_mutation[5] == "")
                )
            ):

                # Merge mutations

                # Combine refs and alts in a new mutation
                new_mutation = (
                    cur_mutation[0],
                    cur_mutation[1],
                    cur_mutation[2],
                    cur_mutation[3],
                    cur_mutation[4] + next_mutation[4],
                    cur_mutation[5] + next_mutation[5],
                )

                # Pop both cur_mutation and next_mutation from the list
                del aa_mutations[i - 1 : i + 1]

                # Insert new mutation at existing index (so it can be merged further if necessary)
                aa_mutations.insert(i - 1, new_mutation)

                # print('MERGE AA', cur_mutation, next_mutation, '-->', new_mutation)

            else:
                # No merging, move on
                cur_mutation = aa_mutations[i]
                i += 1
    """

    # SPLIT AA MUTATIONS
    # ------------------
    # Split AA mutations:
    # 1. Split deletions of multiple residues into individual deletions
    #    - e.g., ∆HV69 --> ∆H69, ∆V70
    #    - Only applicable to pure deletions, i.e., no substitutions
    # 2. Mutations will be “ungrouped” if possible
    #    - e.g., FR157SG --> F157S, R158G
    #    - Only applicable for substitutions of the same residue length

    # 1. SPLIT DELETIONS
    old_aa_mutation_inds = []
    new_aa_mutations = []  # (new_mutation, insertion_index)
    for i, aa_mutation in enumerate(aa_mutations):
        # (reference, Accession ID, gene/protein, pos, ref, alt)
        ref = aa_mutation[4]
        alt = aa_mutation[5]
        # Skip if the mutation is not a pure deletion of more than 1 residue
        if len(alt) > 0 or len(ref) == 1:
            continue

        # Split the deletion into individual deletions
        for j, r in enumerate(ref):
            new_aa_mutations.append(
                (
                    (
                        aa_mutation[0],
                        aa_mutation[1],
                        aa_mutation[2],
                        aa_mutation[3] + j,
                        r,
                        "",
                    ),
                    i + j,
                )
            )

        old_aa_mutation_inds.append(i)

    # Remove old deletion mutations
    for i in old_aa_mutation_inds[::-1]:
        del aa_mutations[i]
    # Add new mutations
    for new_mutation, insertion_index in new_aa_mutations:
        aa_mutations.insert(insertion_index, new_mutation)

    # 2. SPLIT SUBSTITUTIONS
    old_aa_mutation_inds = []
    new_aa_mutations = [] # (new_mutation, insertion_index)
    for i, aa_mutation in enumerate(aa_mutations):
        # (reference, Accession ID, gene/protein, pos, ref, alt)
        ref = aa_mutation[4]
        alt = aa_mutation[5]
        # Skip if the mutation is not a pure substitution of more than 1 residue
        if len(alt) == 0 or len(ref) == 0 or len(ref) != len(alt) or len(ref) == 1:
            continue

        # Split the substitution into individual substitutions
        for j, (a, b) in enumerate(zip(ref, alt)):
            new_aa_mutations.append((
                (
                    aa_mutation[0],
                    aa_mutation[1],
                    aa_mutation[2],
                    aa_mutation[3] + j,
                    a,
                    b,
                ), i + j
            ))

        old_aa_mutation_inds.append(i)

    # Remove old deletion mutations
    for i in old_aa_mutation_inds[::-1]:
        del aa_mutations[i]
    # Add new mutations
    for new_mutation, insertion_index in new_aa_mutations:
        aa_mutations.insert(insertion_index, new_mutation)

    aa_mutation_df = pd.DataFrame.from_records(
        aa_mutations,
        columns=["reference", "Accession ID", "feature", "pos", "ref", "alt"],
    )

    # Isolate accession ID, if its bound with the strain/isolate name
    aa_mutation_df.loc[:, "Accession ID"] = (
        aa_mutation_df["Accession ID"].str.split("|").apply(lambda x: x[0])
    )

    aa_mutation_df.insert(1, "subtype", subtype)
    aa_mutation_df.insert(2, "segment", active_segment)

    return aa_mutation_df


def main():
    """Entry point"""

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--dna-mutation",
        type=str,
        required=True,
        help="Path to DNA mutation CSV",
    )
    parser.add_argument(
        "--gene-protein-def",
        type=str,
        required=True,
        help="Path to gene/protein definition JSON",
    )
    parser.add_argument(
        "--reference", type=str, required=True, help="Path to reference file"
    )
    parser.add_argument("--subtype", type=str, required=True, help="Subtype")
    parser.add_argument("--segment", type=str, required=True, help="Segment/chromosome")
    parser.add_argument(
        "--mode", type=str, required=True, help="Mode ('gene' or 'protein')"
    )
    parser.add_argument("--out", type=str, required=True, help="Path to output")
    args = parser.parse_args()

    aa_mutation_df = extract_aa_mutations(
        args.dna_mutation,
        args.gene_protein_def,
        args.reference,
        args.subtype,
        args.segment,
        args.mode,
    )
    aa_mutation_df.to_csv(args.out, index=False)


if __name__ == "__main__":
    main()
