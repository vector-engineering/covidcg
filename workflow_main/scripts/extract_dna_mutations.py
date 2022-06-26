#!/usr/bin/env python3
# coding: utf-8

"""Extract DNA mutations from bowtie2 alignments

Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""

import argparse
import json
import pandas as pd
import pysam

from read_extractor_lite import ReadExtractor


def extract_dna_mutations(
    sam_file, reference_file, subtype, active_segment,
):
    """
    Extract DNA mutations from bowtie2 alignments

    Parameters
    ----------
    sam_file: str
        Path to SAM file
    reference_file: str
        Path to reference JSON file
    active_segment: str
        Segment/chromosome
    subtype: str
        Subtype

    Returns
    -------
    dna_mutation_df: pandas.DataFrame
    """

    # Load the reference sequences
    with open(reference_file, "r") as fp:
        references = json.loads(fp.read())

    # Get references for this subtype
    references = references[subtype]

    # Get sequences for all references under this subtype,
    # but only for the active segment
    ref_seqs = {
        ref["segments"][active_segment]["name"]: ref["segments"][active_segment][
            "sequence"
        ]
        for ref in references.values()
    }

    ReadExtractor.RefSeq = ref_seqs

    samfile = pysam.AlignmentFile(sam_file, "r")  # pylint: disable=no-member

    all_dna_mutations = []
    for read in samfile.fetch(until_eof=True):
        # Skip if unmapped
        if read.is_unmapped:
            continue

        read_extractor = ReadExtractor(read)

        # print(read.query_name)
        dna_mutations = read_extractor.process_all()
        all_dna_mutations.extend(dna_mutations)

    samfile.close()

    dna_mutation_df = pd.DataFrame.from_records(
        all_dna_mutations, columns=["reference", "Accession ID", "pos", "ref", "alt"]
    )

    # Fill NaN values
    dna_mutation_df["ref"].fillna("", inplace=True)
    dna_mutation_df["alt"].fillna("", inplace=True)

    # Map segment reference names back to genome reference names
    ref_name_map = {
        ref["segments"][active_segment]["name"]: ref["name"]
        for ref in references.values()
    }
    dna_mutation_df["reference"] = dna_mutation_df["reference"].map(ref_name_map)

    return dna_mutation_df


def main():
    """Entry point"""

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--bam", type=str, required=True, help="Path to BAM file",
    )
    parser.add_argument(
        "--reference", type=str, required=True, help="Path to reference JSON file"
    )
    parser.add_argument("--subtype", type=str, required=True, help="Subtype")
    parser.add_argument("--segment", type=str, required=True, help="Segment/chromosome")
    parser.add_argument("--out", type=str, required=True, help="Path to output")
    args = parser.parse_args()

    dna_mutation_df = extract_dna_mutations(
        args.bam, args.reference, args.subtype, args.segment
    )
    dna_mutation_df.to_csv(args.out, index=False)


if __name__ == "__main__":
    main()
