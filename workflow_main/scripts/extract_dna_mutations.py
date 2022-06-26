#!/usr/bin/env python3
# coding: utf-8

"""Extract DNA mutations from bowtie2 alignments

Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""

import argparse
import pandas as pd
import pysam

from read_extractor_lite import ReadExtractor
from fasta import read_fasta_file


def extract_dna_mutations(sam_file, reference_file):
    """
    Extract DNA mutations from bowtie2 alignments

    Parameters
    ----------
    sam_file: str
        Path to SAM file
    reference_file: str
        Path to reference JSON file

    Returns
    -------
    dna_mutation_df: pandas.DataFrame
    """

    # Load the reference sequences
    with open(reference_file, "r") as fp:
        lines = fp.readlines()
        ref_seqs = read_fasta_file(lines)

    ref_seqs = {ref["name"]: ref["sequence"] for ref in ref_seqs.values()}

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

    return dna_mutation_df


def main():
    """Entry point"""

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--bam", type=str, required=True, help="Path to BAM file",
    )
    parser.add_argument(
        "--reference", type=str, required=True, help="Path to reference FASTA file"
    )
    parser.add_argument("--out", type=str, required=True, help="Path to output")
    args = parser.parse_args()

    dna_mutation_df = extract_dna_mutations(args.bam, args.reference)
    dna_mutation_df.to_csv(args.out, index=False)


if __name__ == "__main__":
    main()
