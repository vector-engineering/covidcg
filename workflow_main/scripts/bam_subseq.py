# coding: utf-8

"""
Generate subsequence FASTA files from BAM files
TODO: Restrict to one reference or allow multiple start/ends per reference
      Not really necessary for now since this is only used in the SARS2 workflows
      where there is only one reference (WIV04) ...for now

Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""

import argparse
import gzip
import json
import pysam

from read_extractor_lite import ReadExtractor
from util import translate


def extract_subseqs(bam_file, reference_file, subtype, active_segment, start, end):
    """
    Extract nucleotide subsequence of queries from a BAM file

    Parameters
    ----------
    bam_file: str
        Path to BAM file
    reference_file: str
        Path to reference JSON file
    active_segment: str
        Segment/chromosome
    subtype: str
        Subtype
    start: int
        1-indexed start position
    end: int
        1-indexed end position

    Returns
    -------
    None
    """

    # Load the reference sequences
    with open(reference_file, "r") as fp:
        references = json.loads(fp.read())

    # Get references for this subtype
    references = {k: v for k, v in references.items() if v["subtype"] == subtype}

    # Get sequences for all references under this subtype,
    # but only for the active segment
    ref_seqs = {
        ref["segments"][active_segment]["name"]: ref["segments"][active_segment][
            "sequence"
        ]
        for ref in references.values()
    }

    ReadExtractor.RefSeq = ref_seqs

    bamfile = pysam.AlignmentFile(bam_file, "r")  # pylint: disable=no-member

    subseqs = []
    for read in bamfile.fetch(until_eof=True):
        # Skip if unmapped
        if read.is_unmapped:
            continue

        read_extractor = ReadExtractor(read)
        if not read_extractor.valid:
            continue

        # Crawl to region of interest, then start extracting until the end of the region of interest
        read_extractor.crawl_to(start - 1)  # Start is 1-indexed
        read_extractor.crawl_to(end, store_query_sequence=True)
        subseqs.append((read.query_name, read.reference_name, read_extractor.query_seq))

    bamfile.close()

    return subseqs


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
    parser.add_argument(
        "--start",
        type=int,
        required=True,
        help="Substring start coordinate (1-indexed)",
    )
    parser.add_argument(
        "--end", type=int, required=True, help="Substring end coordinate (1-indexed)"
    )
    parser.add_argument(
        "--output-dna", type=str, default=None, help="Path to output DNA FASTA file"
    )
    parser.add_argument(
        "--output-aa",
        type=str,
        default=None,
        help="Path to output amino acid FASTA file",
    )

    args = parser.parse_args()

    subseqs = extract_subseqs(
        args.bam, args.reference, args.subtype, args.segment, args.start, args.end
    )

    if args.output_dna:
        with gzip.open(args.output_dna, "wt") as fp:
            for subseq in subseqs:
                fp.write(f">{subseq[0]}|{subseq[1]}\n{subseq[2]}\n")

    if args.output_aa:
        with gzip.open(args.output_aa, "wt") as fp:
            for subseq in subseqs:
                fp.write(f">{subseq[0]}|{subseq[1]}\n{translate(subseq[2])}\n")


if __name__ == "__main__":
    main()
