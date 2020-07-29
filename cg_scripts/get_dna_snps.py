import pandas as pd
import pysam
from pathlib import Path

from cg_scripts.fasta import read_fasta_file
from cg_scripts.read_extractor_lite import ReadExtractor


def get_dna_snps(sam_file, reference_file):
    # Load the reference sequence
    with open(reference_file, "r") as fp:
        lines = fp.readlines()
        ref = read_fasta_file(lines)
        ref_seq = list(ref.values())[0]

    ReadExtractor.RefSeq = ref_seq

    samfile = pysam.AlignmentFile(sam_file, "r")  # pylint: disable=no-member

    all_dna_snps = []
    for read in samfile.fetch(until_eof=True):
        # Skip if unmapped
        if read.is_unmapped:
            continue

        read_extractor = ReadExtractor(read)

        # print(read.query_name)
        dna_snps = read_extractor.process_all()
        all_dna_snps.extend(dna_snps)

    samfile.close()

    dna_snp_df = pd.DataFrame.from_records(
        all_dna_snps, columns=["taxon", "pos", "ref", "alt"]
    )

    # Extract the GISAID ID
    dna_snp_df["gisaid_id"] = dna_snp_df["taxon"].str.split("|", expand=True)[1]

    # Fill NaN values
    dna_snp_df["ref"].fillna("", inplace=True)
    dna_snp_df["alt"].fillna("", inplace=True)

    # Drop duplicate entries
    # dna_snp_df.drop_duplicates(["taxon", "pos"], inplace=True)
    # dna_snp_df = dna_snp_df.reset_index(drop=True)

    return dna_snp_df

