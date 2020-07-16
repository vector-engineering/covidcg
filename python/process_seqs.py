#!/usr/bin/env python3
# coding: utf-8

"""Process raw sequences from GISAID
"""

import math
import multiprocessing as mp
import os
import pysam
import re
import subprocess

from functools import partial

from fasta import read_fasta_file
from process_snps import load_dna_snps, process_gene_aa_snps, process_protein_aa_snps
from read_extractor_lite import ReadExtractor
from util import data_dir, static_data_dir


def preprocess_sequences():
    """Filter out sequences (adapted from van Dorp et al, 2020)
    1. Filter against nextstrain exclusion list
    2. Remove animal/environmental isolates (bat, pangolin, mink, tiger, cat, canine, env)
    3. Can't be less than 29700NT
	4. Can't have more than 5% ambiguous NT
    """

    print("\nPreprocessing sequences")

    # Get latest nextstrain exclusion file
    ns_exclusion_list = sorted(static_data_dir.glob("nextstrain_exclude*"))[-1]
    print("Using nextstrain exclusion list: {}".format(ns_exclusion_list.name))
    # Load lines, ignoring comments and empty lines
    exclude_taxons = []
    with ns_exclusion_list.open("r") as fp:
        for line in fp.readlines():
            # Exclude comments
            if line[0] == "#":
                continue

            # Strip whitespace
            line = re.sub(r"\s+", "", line).strip()

            # Exclude empty lines
            if not line:
                continue

            exclude_taxons.append(line)

    animal_env_tags = ["bat", "pangolin", "mink", "tiger", "cat", "canine", "env"]
    # Surround tags with '/' for easier and more explicit matching to taxons
    animal_env_tags = ["/{}/".format(t) for t in animal_env_tags]

    fasta_files = [
        f
        for f in sorted((data_dir / "fasta_raw").iterdir())
        if f.suffix == ".fasta" or f.suffix == ".fa"
    ]

    # Make output directory
    (data_dir / "fasta_processed").mkdir(exist_ok=True)

    for ff in fasta_files:
        # Skip if the sam file already exists
        processed_file_path = data_dir / "fasta_processed" / ff.name
        if processed_file_path.exists():
            print("{} already exists, skipping".format(processed_file_path.name))
            continue

        print("Preprocessing {}...".format(ff.name), end="", flush=True)

        num_excluded = 0
        fp_in = ff.open("r")
        fp_out = processed_file_path.open("w")

        cur_entry = ""
        cur_seq = ""
        while True:
            line = fp_in.readline()

            # Beginning of a new entry, or EOF = end of current entry
            if not line or line[0] == ">":

                if cur_entry:
                    num_ambiguous = 0
                    for char in cur_seq:
                        if char == "N":
                            num_ambiguous += 1

                    if (
                        # 1: Check against nextstrain exclusion list
                        (cur_entry in exclude_taxons)
                        or
                        # 2: Remove animal/environmental isolates
                        any([tag in cur_entry for tag in animal_env_tags])
                        or
                        # 3: Can't be less than 29700 NT
                        len(cur_seq) < 29700
                        or
                        # 4: Can't have more than 5% ambiguous (N) NT
                        num_ambiguous > math.floor(len(cur_seq) * 0.05)
                    ):
                        num_excluded += 1
                    else:
                        # It passed, write to output
                        fp_out.write(">" + cur_entry + "\n")
                        fp_out.write(cur_seq + "\n")

                # If it's the end, then break out
                if not line:
                    break

                # Reset sequence and name
                cur_seq = ""
                # Extract the name (up to the first whitespace)
                # [1:] excludes the first '>'
                # .split() breaks up the line into chunks separated by whitespace
                # [0] gets the first chunk
                # cur_entry = line[1:].split()[0]
                # Nevermind, the fasta entries sometimes have spaces.....
                # Just rstrip to remove the newline, that should work good enough
                cur_entry = line[1:].rstrip()

            # Otherwise add sequence to the current entry
            elif cur_entry:
                cur_seq += re.sub(r"\s+", "", line).strip()

        fp_in.close()
        fp_out.close()

        print("done. Removed {:,} sequences".format(num_excluded), flush=True)
    # END FOR FASTA FILE


def align_sequences():

    print("\nAligning sequences")

    # Build bowtie2 index, but skip if already built
    bt2_index_dir = data_dir / "reference_index"
    bt2_index_name_path = bt2_index_dir / "reference"
    ref_fasta_path = static_data_dir / "reference.fasta"

    if not (
        bt2_index_dir.exists()
        and bt2_index_dir.is_dir()
        and bt2_index_dir.glob("*.bt2")
    ):
        bt2_index_dir.mkdir(exist_ok=True)
        subprocess.run(
            ["bowtie2-build", str(ref_fasta_path), str(bt2_index_name_path)],
            stdout=subprocess.DEVNULL,
        )  # hide stdout

    fasta_files = [
        f
        for f in sorted((data_dir / "fasta_processed").iterdir())
        if f.suffix == ".fasta" or f.suffix == ".fa"
    ]
    # sam_files = [f for f in sorted((data_dir / 'sam').iterdir()) if f.suffix == '.sam']

    # Make output directory
    (data_dir / "sam").mkdir(exist_ok=True)

    # Align sequences to reference with bowtie2
    for i, ff in enumerate(fasta_files):

        # Testing
        # if i > 0:
        #     break

        # Skip if the sam file already exists
        sam_file_path = data_dir / "sam" / (ff.stem + ".sam")
        if sam_file_path.exists():
            print("{} already exists, skipping".format(sam_file_path.name))
            continue

        print(
            "Aligning {}. This will take a while, and will benefit from multiprocessing with multiple CPUs/threads. Note that each thread will require ~7.5GB of RAM to run, so please do not specify too many CPUs or else bowtie2 will crash.".format(
                ff.name
            )
        )
        subprocess.run(
            [
                "bowtie2",
                "--end-to-end",
                "--very-fast",  # Global alignment, no overhangs
                "--xeq",  # Track SNPs
                "--reorder",  # Same order as input fasta file
                "--sam-no-qname-trunc",  # Don't truncate entry names
                "-x",
                str(bt2_index_name_path),  # Path to index we just built
                "-f",  # FASTA input
                "-U",
                str(ff),
                "-S",
                str(sam_file_path),  # Output path
                # '--threads', str(os.cpu_count()) # Blast it
                "--threads",
                "1",  # TODO: add this to CLI args
            ]
        )


def process_snps_from_reads(sam_file_path, start, end):

    samfile = pysam.AlignmentFile(sam_file_path, "r")  # pylint: disable=no-member

    all_dna_snps = []
    counter = 0
    for read in samfile.fetch(until_eof=True):

        if counter < start or counter >= end:
            counter += 1
            continue

        # Skip if unmapped
        if read.is_unmapped:
            counter += 1
            continue

        read_extractor = ReadExtractor(read)

        # print(read.query_name)
        dna_snps = read_extractor.process_all()

        all_dna_snps.extend(dna_snps)

        counter += 1

    samfile.close()

    return all_dna_snps


def process_snps(num_processes=0):
    print("\nGetting SNPs/indels")
    sam_files = sorted((data_dir / "sam").glob("*.sam"))

    dna_snp_folder_path = data_dir / "dna_snp"
    dna_snp_folder_path.mkdir(exist_ok=True)

    # Get SNPs from each sam file
    for i, sam_file_path in enumerate(sam_files):
        # Output paths
        dna_snp_path = dna_snp_folder_path / (sam_file_path.stem + "_dna_snp.csv")
        # Skip if already done
        if dna_snp_path.exists():
            print("{} already exist, skipping".format(dna_snp_path.name))
            continue

        # Testing
        # if i > 0:
        #     break

        print(
            "Getting SNPs/indels from {}...".format(sam_file_path.name),
            end="",
            flush=True,
        )

        samfile = pysam.AlignmentFile(
            str(sam_file_path), "r"
        )  # pylint: disable=no-member
        num_reads = samfile.count(until_eof=True)
        samfile.close()

        # TODO: make this a CLI arg
        if num_processes == 0:
            num_processes = os.cpu_count()
        pool = mp.Pool(processes=num_processes)

        all_dna_snps = []

        def process_callback(ret):
            all_dna_snps.extend(ret)

        def handle_error(e):
            raise e

        chunk_size = math.ceil(num_reads / num_processes)
        for chunk_start in range(0, num_reads, chunk_size):
            pool.apply_async(
                partial(
                    process_snps_from_reads,
                    str(sam_file_path),
                    chunk_start,
                    chunk_start + chunk_size,
                ),
                callback=process_callback,
                error_callback=handle_error,
            )

        pool.close()
        pool.join()

        # Write to disk
        # Open output files for writing
        fp_dna = dna_snp_path.open("w")
        # Write headers
        fp_dna.write("taxon,pos,ref,alt\n")

        # Write DNA SNPs/indels
        for snp in all_dna_snps:
            fp_dna.write(
                "{},{},{},{}\n".format(snp[0], str(snp[1]), str(snp[2]), str(snp[3]))
            )

        fp_dna.close()

        print("done")


def main():
    preprocess_sequences()
    align_sequences()
    process_snps()

    dna_snp_df = load_dna_snps()
    # print(dna_snp_df)

    # Filter out indels and SNPs with length > 1
    # Need to figure out what to do with those...
    dna_snp_df = (
        dna_snp_df.fillna("")
        .loc[(dna_snp_df["ref"].str.len() == 1) & (dna_snp_df["alt"].str.len() == 1), :]
        .reset_index(drop=True)
    )

    gene_aa_snp_df = process_gene_aa_snps(dna_snp_df)
    # print(gene_aa_snp_df)
    protein_aa_snp_df = process_protein_aa_snps(dna_snp_df)
    # print(protein_aa_snp_df)

    gene_aa_snp_df.to_csv(data_dir / "gene_aa_snp.csv", index=False)
    protein_aa_snp_df.to_csv(data_dir / "protein_aa_snp.csv", index=False)


if __name__ == "__main__":
    main()
