#!/usr/bin/env python3
# coding: utf-8

"""Extract DNA mutations from bowtie2 alignments

Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""

import argparse
import pandas as pd
import pysam
<<<<<<< HEAD
=======
from pathlib import Path
>>>>>>> 5cdb66c3 (Updated to work with new server layout)

from scripts.fasta import read_fasta_file
from scripts.read_extractor_lite import ReadExtractor


def extract_dna_mutations(sam_file, reference_file):
    # Load the reference sequence
    with open(reference_file, "r") as fp:
        lines = fp.readlines()
        ref_seq = read_fasta_file(lines)

    ReadExtractor.RefSeq = ref_seq

    samfile = pysam.AlignmentFile(sam_file, "r")  # pylint: disable=no-member

    all_dna_mutations = []
    for read in samfile.fetch(until_eof=True):
        # Skip if unmapped
        if read.is_unmapped:
            continue

        read_extractor = ReadExtractor(read)

        # print(read.query_name)
        dna_mutations = read_extractor.process_all()
        if dna_mutations:
            all_dna_mutations.extend(dna_mutations)

    samfile.close()

    dna_mutation_df = pd.DataFrame.from_records(
        all_dna_mutations, columns=["Accession ID", "pos", "ref", "alt", "ref_seq_name"]
    )

    # Fill NaN values
    dna_mutation_df["ref"].fillna("", inplace=True)
    dna_mutation_df["alt"].fillna("", inplace=True)

    return dna_mutation_df


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--bam", type=str, required=True, help="Path to BAM file",
    )
    parser.add_argument(
        "--reference", type=str, required=True, help="Path to reference file"
    )
    parser.add_argument("--out", type=str, required=True, help="Path to output")
    args = parser.parse_args()

    dna_mutation_df = extract_dna_mutations(args.bam, args.reference)
    dna_mutation_df.to_csv(args.out, index=False)


if __name__ == "__main__":
    main()
